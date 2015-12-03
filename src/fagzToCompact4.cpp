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

int fagzToCompact4(libmaus2::util::ArgInfo const & arginfo)
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
	std::string const metaoutputfilename = outputfilename + ".meta";
	int const verbose = arginfo.getValue<int>("verbose",1);
	libmaus2::autoarray::AutoArray<char> B(8*1024,false);
	libmaus2::bitio::CompactArrayWriterFile compactout(outputfilename,2 /* bits per symbol */);

	if ( ! rc )
		std::cerr << "[V] not storing reverse complements" << std::endl;

	// forward mapping table
	libmaus2::autoarray::AutoArray<uint8_t> ftable(256,false);
	// rc mapping for mapped symbols
	libmaus2::autoarray::AutoArray<uint8_t> ctable(256,false);

	std::fill(ftable.begin(),ftable.end(),4);
	std::fill(ctable.begin(),ctable.end(),4);
	ftable['a'] = ftable['A'] = 0;
	ftable['c'] = ftable['C'] = 1;
	ftable['g'] = ftable['G'] = 2;
	ftable['t'] = ftable['T'] = 3;
	uint64_t insize = 0;

	ctable[0] = 3; // A->T
	ctable[1] = 2; // C->G
	ctable[2] = 1; // G->C
	ctable[3] = 0; // T->A

	libmaus2::aio::OutputStreamInstance::unique_ptr_type metaOut(new libmaus2::aio::OutputStreamInstance(metaoutputfilename));
	libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,0);
	uint64_t nseq = 0;
	std::vector<uint64_t> lvec;

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

			libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,pattern.spattern.size());
			lvec.push_back(pattern.spattern.size());
			libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,0);

			// map symbols
			for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
				pattern.spattern[j] = ftable[static_cast<uint8_t>(pattern.spattern[j])];

			uint64_t l = 0;
			uint64_t nr = 0;
			while ( l < pattern.spattern.size() )
			{
				while ( l < pattern.spattern.size() && pattern.spattern[l] < 4 )
					++l;
				assert ( l == pattern.spattern.size() || pattern.spattern[l] == 4 );

				uint64_t h = l;
				while ( h < pattern.spattern.size() && pattern.spattern[h] == 4 )
					++h;

				if ( h-l )
				{
					for ( uint64_t j = l; j < h; ++j )
						pattern.spattern[j] = (libmaus2::random::Random::rand8() & 3);

					libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,l);
					libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,h);
					nr += 1;
				}

				l = h;
			}

			for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
			{
				assert ( pattern.spattern[j] < 4 );
			}

			metaOut->seekp( - static_cast<int64_t>(2*nr+1)*sizeof(uint64_t), std::ios::cur );
			libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,nr);
			metaOut->seekp(   static_cast<int64_t>(2*nr  )*sizeof(uint64_t), std::ios::cur );

			// write
			compactout.write(pattern.spattern.c_str(),pattern.spattern.size());

			if ( rc )
			{
				// reverse complement
				std::reverse(pattern.spattern.begin(),pattern.spattern.end());
				for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
					pattern.spattern[j] = ctable[static_cast<uint8_t>(pattern.spattern[j])];

				// write
				compactout.write(pattern.spattern.c_str(),pattern.spattern.size());
			}

			insize += pattern.spattern.size()+1;
			nseq += 1;

			if ( verbose )
				std::cerr << "done, input size " << formatBytes(pattern.spattern.size()+1) << " acc " << formatBytes(insize) << std::endl;
		}
	}

	metaOut->seekp(0);
	libmaus2::util::NumberSerialisation::serialiseNumber(*metaOut,nseq);
	metaOut->flush();
	metaOut.reset();

	libmaus2::aio::InputStreamInstance::unique_ptr_type metaISI(new libmaus2::aio::InputStreamInstance(metaoutputfilename));
	// number of sequences
	uint64_t const rnseq = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
	assert ( nseq == rnseq );
	for ( uint64_t i = 0; i < nseq; ++i )
	{
		// length of sequence
		uint64_t const l = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		assert ( l == lvec[i] );
		uint64_t const nr = libmaus2::util::NumberSerialisation::deserialiseNumber(*metaISI);
		// skip replaced intervals
		metaISI->ignore(2*nr*sizeof(uint64_t));
	}
	assert ( metaISI->peek() == std::istream::traits_type::eof() );

	std::cerr << "Done, total input size " << insize << std::endl;

	compactout.flush();

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		return fagzToCompact4(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
