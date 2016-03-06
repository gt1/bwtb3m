/**
    feramanzgen
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
#include <libmaus2/bitio/CompactArrayWriterFile.hpp>
#include <libmaus2/fastx/CharBuffer.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/lz/BufferedGzipStream.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/aio/InputStreamFactoryContainer.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/aio/CircularWrapper.hpp>
#include <libmaus2/fastx/DNAIndexMetaDataBigBandBiDir.hpp>
#include <libmaus2/bitio/CompactArray.hpp>

/*
 * produce 2 bit representation plus meta data from (compressed) fasta file
 */

std::string formatBytes(uint64_t n)
{
	char const * units[] = { "", "k", "m", "g", "t", "p", "e", "z", "y" };

	std::vector < std::string > parts;
	unsigned int uindex = 0;

	for ( ; n; ++uindex )
	{
		std::ostringstream ostr;
		ostr << (n%1024) << units[uindex];
		parts.push_back(ostr.str());

		n /= 1024;
	}

	std::reverse(parts.begin(),parts.end());

	std::ostringstream ostr;
	for ( uint64_t i = 0; i < parts.size(); ++i )
	{
		ostr << parts[i];
		if ( i+1 < parts.size() )
			ostr << " ";
	}

	return ostr.str();
}

std::string basename(std::string const s)
{
	uint64_t l = s.size();
	for ( uint64_t i = 0; i < s.size(); ++i )
		if ( s[i] == '/' )
			l = i;

	if ( l != s.size() )
		return s.substr(l+1);
	else
		return s;
}

std::string stripAfterDot(std::string const s)
{
	for ( uint64_t i = 0; i < s.size(); ++i )
		if ( s[i] == '.' )
			return s.substr(0,i);

	return s;
}

int fagzToCompact4BigBandBiDir(libmaus2::util::ArgInfo const & arginfo)
{
	bool const rc = arginfo.getValue<unsigned int>("rc",1);
	bool const gz = arginfo.getValue<unsigned int>("gz",1);

	std::vector<std::string> inputfilenames;
	inputfilenames = arginfo.restargs;

	if ( arginfo.hasArg("inputfilenames") )
	{
		std::string const inf = arginfo.getUnparsedValue("inputfilenames",std::string());
		libmaus2::aio::InputStream::unique_ptr_type Pinf(libmaus2::aio::InputStreamFactoryContainer::constructUnique(inf));
		while ( *Pinf )
		{
			std::string line;
			std::getline(*Pinf,line);
			if ( line.size() )
				inputfilenames.push_back(line);
		}
	}

	std::string const inlcp = libmaus2::util::OutputFileNameTools::lcp(inputfilenames);
	std::string defout = inlcp;
	defout = libmaus2::util::OutputFileNameTools::clipOff(defout,".gz");
	defout = libmaus2::util::OutputFileNameTools::clipOff(defout,".fasta");
	defout = libmaus2::util::OutputFileNameTools::clipOff(defout,".fa");

	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename",defout + ".compact");
	std::string const outputfilenameforw = outputfilename + ".forw";
	std::string const outputfilenamereco = outputfilename + ".reco";
	std::string const finalmetaoutputfilename = outputfilename + ".meta";
	std::string const metaoutputfilename = outputfilename + ".meta.tmp";
	std::string const metametaoutputfilename = outputfilename + ".meta.meta";
	std::string const repfastaoutputfilename = outputfilename + ".repl.fasta";
	std::string const seqtmpfilename = outputfilename + ".repl.fasta.tmp";
	int const verbose = arginfo.getValue<int>("verbose",1);

	libmaus2::aio::OutputStreamInstance::unique_ptr_type replOSI(new libmaus2::aio::OutputStreamInstance(repfastaoutputfilename));
	
	libmaus2::util::TempFileRemovalContainer::addTempFile(seqtmpfilename);
	libmaus2::util::TempFileRemovalContainer::addTempFile(metametaoutputfilename);

	if ( ! rc )
		std::cerr << "[V] not storing reverse complements" << std::endl;

	// forward mapping table
	libmaus2::autoarray::AutoArray<uint8_t> ftable(256,false);
	// rc mapping for mapped symbols
	libmaus2::autoarray::AutoArray<uint8_t> ctable(256,false);
	// backward mapping table
	libmaus2::autoarray::AutoArray<uint8_t> rtable(256,false);

	std::fill(ftable.begin(),ftable.end(),4);
	std::fill(ctable.begin(),ctable.end(),4);
	std::fill(rtable.begin(),rtable.end(),'N');
	ftable['a'] = ftable['A'] = 0;
	ftable['c'] = ftable['C'] = 1;
	ftable['g'] = ftable['G'] = 2;
	ftable['t'] = ftable['T'] = 3;

	ctable[0] = 3; // A->T
	ctable[1] = 2; // C->G
	ctable[2] = 1; // G->C
	ctable[3] = 0; // T->A

	rtable[0] = 'A';
	rtable[1] = 'C';
	rtable[2] = 'G';
	rtable[3] = 'T';

	libmaus2::aio::OutputStreamInstance::unique_ptr_type metaOut(new libmaus2::aio::OutputStreamInstance(metaoutputfilename));
	libmaus2::aio::OutputStreamInstance::unique_ptr_type metametaOut(new libmaus2::aio::OutputStreamInstance(metametaoutputfilename));

	libmaus2::aio::OutputStreamInstance::unique_ptr_type seqtmpOut(new libmaus2::aio::OutputStreamInstance(seqtmpfilename));
	libmaus2::aio::OutputStreamInstance::unique_ptr_type repfastaOut(new libmaus2::aio::OutputStreamInstance(repfastaoutputfilename));

	// number of sequences
	uint64_t metapos = 0;
	libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,0);
	metapos += sizeof(uint64_t);

	std::vector<uint64_t> lvec;
	uint64_t nseq = 0;
	uint64_t seqoff = 0;
	for ( uint64_t i = 0; i < inputfilenames.size(); ++i )
	{
		std::string const fn = inputfilenames[i];
		libmaus2::aio::InputStreamInstance CIS(fn);
		libmaus2::lz::BufferedGzipStream::unique_ptr_type BGS;
		std::istream * istr = 0;
		if ( gz )
		{
			libmaus2::lz::BufferedGzipStream::unique_ptr_type tBGS(
				new libmaus2::lz::BufferedGzipStream(CIS));
			BGS = UNIQUE_PTR_MOVE(tBGS);
			istr = BGS.get();
		}
		else
		{
			istr = &CIS;
		}
		libmaus2::fastx::StreamFastAReaderWrapper fain(*istr);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

		while ( fain.getNextPatternUnlocked(pattern) )
		{
			if ( verbose )
				std::cerr << (i+1) << " " << stripAfterDot(basename(fn)) << " " << pattern.sid << "...";
				
			std::string const & id = pattern.sid;
			std::string & spattern = pattern.spattern;
			uint64_t const seqlen = spattern.size();

			assert ( static_cast<int64_t>(metapos) == metaOut->tellp() );
			libmaus2::util::NumberSerialisation::serialiseNumber(*metametaOut,metapos);
			libmaus2::util::NumberSerialisation::serialiseNumber(*metametaOut,seqoff);
			
			// length of sequence
			libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,seqlen);
			lvec.push_back(seqlen);
			
			for ( uint64_t i = 0; i < seqlen; ++i )
				spattern[i] = ftable[spattern[i]];

			std::vector < std::pair<uint64_t,uint64_t > > Vreplace;

			uint64_t low = 0;
			while ( low < seqlen )
			{
				while ( low < seqlen && spattern[low] < 4 )
					++low;
				
				if ( low < seqlen )
				{
					assert ( spattern[low] >= 4 );
					uint64_t high = low+1;
					while ( high < seqlen && spattern[high] >= 4 )
						++high;
					for ( uint64_t i = low; i < high; ++i )
						assert ( spattern[i] >= 4 );
					if ( low )
						assert ( spattern[low-1] < 4 );
					if ( high < seqlen )
						assert ( spattern[high] < 4 );
						
					for ( uint64_t i = low; i < high; ++i )
						pattern.spattern[i] = (libmaus2::random::Random::rand8() & 3);
						
					Vreplace.push_back(std::pair<uint64_t,uint64_t>(low,high));

					low = high;
				}
			}
			
			assert ( seqtmpOut->tellp() == static_cast<int64_t>(seqoff) );
			seqtmpOut->write(spattern.c_str(),seqlen);
			seqoff += seqlen;

			// map back to clear text
			for ( uint64_t i = 0; i < seqlen; ++i )
				spattern[i] = rtable[spattern[i]];

			// write recoded sequence
			char const * rout = spattern.c_str();
			char const * route = rout + seqlen;
			(*repfastaOut) << '>' << id << '\n';
			while ( rout != route )
			{
				uint64_t const rest = route-rout;
				uint64_t const cols = 80;
				uint64_t const towrite = std::min(rest,cols);
				
				repfastaOut->write(rout,towrite);
				repfastaOut->put('\n');
				
				rout += towrite;
			}
			
			// number of replaced intervals
			libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,Vreplace.size());
			// replaced intervals
			for ( uint64_t i = 0; i < Vreplace.size(); ++i )
			{
				libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,Vreplace[i].first);
				libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,Vreplace[i].second);
			}

			metapos += 2*sizeof(uint64_t) + 2 * Vreplace.size() * sizeof(uint64_t);
			assert ( static_cast<int64_t>(metapos) == metaOut->tellp() );			
			
			nseq += 1;

			if ( verbose )
				std::cerr << std::endl;
		}
	}

	metaOut->seekp(0);
	libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,nseq);
	metaOut->flush();
	metaOut.reset();
	metametaOut->flush();
	metametaOut.reset();
	
	seqtmpOut->flush();
	seqtmpOut.reset();

	libmaus2::aio::InputStreamInstance::unique_ptr_type metametaISI(new libmaus2::aio::InputStreamInstance(metametaoutputfilename));
	libmaus2::aio::InputStreamInstance::unique_ptr_type metaISI(new libmaus2::aio::InputStreamInstance(metaoutputfilename));
	libmaus2::aio::InputStreamInstance::unique_ptr_type seqtmpISI(new libmaus2::aio::InputStreamInstance(seqtmpfilename));

	libmaus2::bitio::CompactArrayWriterFile::unique_ptr_type Pcompactout(
		new libmaus2::bitio::CompactArrayWriterFile(outputfilename,2 /* bits per symbol */)
	);
	libmaus2::bitio::CompactArrayWriterFile::unique_ptr_type Pcompactforwout(
		new libmaus2::bitio::CompactArrayWriterFile(outputfilenameforw,2 /* bits per symbol */)
	);
	libmaus2::bitio::CompactArrayWriterFile::unique_ptr_type Pcompactrecoout(
		new libmaus2::bitio::CompactArrayWriterFile(outputfilenamereco,2 /* bits per symbol */)
	);
	
	libmaus2::aio::OutputStreamInstance::unique_ptr_type finalmetaOut(new libmaus2::aio::OutputStreamInstance(finalmetaoutputfilename));
	libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,2*nseq);

	libmaus2::autoarray::AutoArray<char> B(64*1024,false);
	for ( uint64_t i = 0; i < nseq; ++i )
	{
		metametaISI->clear();
		metametaISI->seekg(i * 2 * sizeof(uint64_t));
		
		uint64_t const metapos = libmaus2::util::NumberSerialisation::deserialiseNumber(*metametaISI);
		uint64_t const seqoff = libmaus2::util::NumberSerialisation::deserialiseNumber(*metametaISI);
		
		metaISI->clear();
		metaISI->seekg(metapos);

		uint64_t const seqlen = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		uint64_t const replintv = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		
		std::vector < std::pair<uint64_t,uint64_t> > Vreplace(replintv);
		for ( uint64_t j = 0; j < replintv; ++j )
		{
			Vreplace[j].first = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
			Vreplace[j].second = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		}
		
		// write meta data
		libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,seqlen);
		libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,replintv);
		for ( uint64_t j = 0; j < replintv; ++j )
		{
			libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,Vreplace[j].first);
			libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,Vreplace[j].second);
		}
				
		
		seqtmpISI->clear();
		seqtmpISI->seekg(seqoff);
		
		uint64_t todo = seqlen;
		
		while ( todo )
		{
			uint64_t const tocopy = std::min(todo,B.size());
			
			seqtmpISI->read(B.begin(),tocopy);
			assert ( seqtmpISI->gcount() == static_cast<int64_t>(tocopy) );
			Pcompactout->write(B.begin(),tocopy);
			Pcompactforwout->write(B.begin(),tocopy);
			
			todo -= tocopy;
		}

		if ( verbose )
			std::cerr << "[V] rewriting forward of sequence " << i << " length " << seqlen << " replintv=" << replintv << std::endl;
	}
	seqtmpISI.reset();
	Pcompactforwout->flush();
	Pcompactforwout.reset();

	// reverse complement
	for ( uint64_t ii = 0; ii < nseq; ++ii )
	{
		uint64_t const i = nseq-ii-1;
	
		metametaISI->clear();
		metametaISI->seekg(i * 2 * sizeof(uint64_t));
		
		uint64_t const metapos = libmaus2::util::NumberSerialisation::deserialiseNumber(*metametaISI);
		uint64_t const seqoff = libmaus2::util::NumberSerialisation::deserialiseNumber(*metametaISI);
		
		metaISI->clear();
		metaISI->seekg(metapos);

		uint64_t const seqlen = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		uint64_t const replintv = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		
		std::vector < std::pair<uint64_t,uint64_t> > Vreplace(replintv);
		for ( uint64_t j = 0; j < replintv; ++j )
		{
			Vreplace[j].first = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
			Vreplace[j].second = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		}
		
		// write meta data
		libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,seqlen);
		libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,replintv);
		for ( uint64_t jj = 0; jj < replintv; ++jj )
		{
			uint64_t const j = replintv - jj - 1;
			
			libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,seqlen-Vreplace[j].second);
			libmaus2::util::NumberSerialisation::serialiseNumber(*finalmetaOut,seqlen-Vreplace[j].first);
		}
				
		libmaus2::aio::CircularReverseWrapper seqtmp(seqtmpfilename,seqoff+seqlen);

		uint64_t todo = seqlen;		
		while ( todo )
		{
			uint64_t const tocopy = std::min(todo,B.size());
						
			seqtmp.read(B.begin(),tocopy);
			
			for ( uint64_t j = 0; j < tocopy; ++j )
				B[j] ^= 3;
			
			Pcompactout->write(B.begin(),tocopy);
			Pcompactrecoout->write(B.begin(),tocopy);
			
			todo -= tocopy;
		}

		if ( verbose )
			std::cerr << "[V] rewriting RC of sequence " << i << " length " << seqlen << " replintv=" << replintv << std::endl;
	}

	Pcompactrecoout->flush();
	Pcompactrecoout.reset();
	
	metametaISI.reset();
	libmaus2::aio::FileRemoval::removeFile(metametaoutputfilename);
	metaISI.reset();
	libmaus2::aio::FileRemoval::removeFile(metaoutputfilename);
	seqtmpISI.reset();
	libmaus2::aio::FileRemoval::removeFile(seqtmpfilename);

	Pcompactout->flush();
	Pcompactout.reset();
	
	finalmetaOut->flush();
	finalmetaOut.reset();

	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::unique_ptr_type Pindex(
		libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::load(finalmetaoutputfilename)
	);
	
	assert ( Pindex->S.size() == 2*nseq );
		
	for ( uint64_t i = 0; i < nseq; ++i )
	{
		libmaus2::fastx::DNAIndexMetaDataSequence const & forw = Pindex->S[i];
		libmaus2::fastx::DNAIndexMetaDataSequence const & reve = Pindex->S[2*nseq-i-1];

		assert ( forw.l == reve.l );
		assert ( forw.nblocks.size() == reve.nblocks.size() );
		
		for ( uint64_t j = 0; j < forw.nblocks.size(); ++j )
		{
		
			std::pair<uint64_t,uint64_t> const & fblock = forw.nblocks[j];
			std::pair<uint64_t,uint64_t> const & rblock = reve.nblocks[reve.nblocks.size()-j-1];
			
			assert ( rblock.second == forw.l - fblock.first );
			assert ( rblock.first  == forw.l - fblock.second );
		}
	}
	
	libmaus2::bitio::CompactArray::unique_ptr_type CA(libmaus2::bitio::CompactArray::load(outputfilename));

	uint64_t suml = 0;
	for ( uint64_t i = 0; i < nseq; ++i )
	{
		libmaus2::fastx::DNAIndexMetaDataSequence const & forw = Pindex->S[i];
		// libmaus2::fastx::DNAIndexMetaDataSequence const & reve = Pindex->S[2*nseq-i-1];
		
		uint64_t pf = Pindex->L[i];
		uint64_t pr = Pindex->L[2*nseq-i-1] + forw.l;
		
		for ( uint64_t j = 0; j < forw.l; ++j )
			assert ( (*CA)[pf++] == ((*CA)[--pr]^3) );
		
		suml += forw.l;
	}
	
	for ( uint64_t i = 0; i < suml; ++i )
		assert (  (*CA)[i] == ((*CA)[2*suml-i-1]^3) );
	
	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		return fagzToCompact4BigBandBiDir(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
