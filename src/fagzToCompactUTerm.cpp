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

template<typename writer_type>
void putSeqId(writer_type & writer, uint64_t const seqid, unsigned int const seqbits)
{
	for ( unsigned int i = 0; i < seqbits; ++i )
	{
		uint64_t const v = (seqid >> (seqbits - i - 1)) & 1;
		writer.put( v );
	}
}

int fagzToCompact(libmaus2::util::ArgInfo const & arginfo)
{
	bool const rc = arginfo.getValue<unsigned int>("rc",1);
	bool const gz = arginfo.getValue<unsigned int>("gz",1);
	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
	int const verbose = arginfo.getValue<int>("verbose",1);
	libmaus2::autoarray::AutoArray<char> B(8*1024,false);
	// A,C,G,T,N,0,1
	libmaus2::bitio::CompactArrayWriterFile compactout(outputfilename,3);

	if ( ! rc )
		std::cerr << "[V] not storing reverse complements" << std::endl;

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

	uint64_t numseq = 0;
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
			numseq += 1;
	}
	if ( rc )
		numseq *= 2;

	unsigned int const seqbits = numseq ? libmaus2::math::numbits(numseq-1) : 0;

	std::cerr << "[V] numseq=" << numseq << std::endl;

	// forward mapping table
	libmaus2::autoarray::AutoArray<uint8_t> ftable(256,false);
	// reverse complement mapping table
	libmaus2::autoarray::AutoArray<uint8_t> rtable(256,false);
	// rc mapping for mapped symbols
	libmaus2::autoarray::AutoArray<uint8_t> ctable(256,false);

	std::fill(ftable.begin(),ftable.end(),6);
	std::fill(rtable.begin(),rtable.end(),6);
	std::fill(ctable.begin(),ctable.end(),6);
	ftable['a'] = ftable['A'] = rtable['t'] = rtable['T'] = 2;
	ftable['c'] = ftable['C'] = rtable['g'] = rtable['G'] = 3;
	ftable['g'] = ftable['G'] = rtable['c'] = rtable['C'] = 4;
	ftable['t'] = ftable['T'] = rtable['a'] = rtable['A'] = 5;

	ctable[2] = 5; // A->T
	ctable[3] = 4; // C->G
	ctable[4] = 3; // G->C
	ctable[5] = 2; // T->A

	uint64_t seqid = 0;
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

			// map symbols
			for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
				pattern.spattern[j] = ftable[static_cast<uint8_t>(pattern.spattern[j])];

			// write
			compactout.write(pattern.spattern.c_str(),pattern.spattern.size());
			putSeqId(compactout, seqid++, seqbits);

			if ( rc )
			{
				// reverse complement
				std::reverse(pattern.spattern.begin(),pattern.spattern.end());
				for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
					pattern.spattern[j] = ctable[static_cast<uint8_t>(pattern.spattern[j])];

				// write
				compactout.write(pattern.spattern.c_str(),pattern.spattern.size());
				putSeqId(compactout, seqid++, seqbits);
			}

			if ( verbose )
				std::cerr << "done, input size " << formatBytes(pattern.spattern.size()) << std::endl;
		}
	}

	compactout.flush();

	assert ( seqid == numseq );

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		return fagzToCompact(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
