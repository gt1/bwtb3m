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
#include <libmaus/bitio/CompactArrayWriterFile.hpp>
#include <libmaus/fastx/CharBuffer.hpp>
#include <libmaus/fastx/StreamFastAReader.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>
#include <libmaus/util/ArgInfo.hpp>

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

int fagzToCompact(libmaus::util::ArgInfo const & arginfo)
{
	bool const rc = arginfo.getValue<unsigned int>("rc",1);
	bool const gz = arginfo.getValue<unsigned int>("gz",1);
	uint64_t const limit = arginfo.getValueUnsignedNumeric<uint64_t>("limit",std::numeric_limits<uint64_t>::max());
	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
	int const verbose = arginfo.getValue<int>("verbose",1);
	libmaus::autoarray::AutoArray<char> B(8*1024,false);
	libmaus::bitio::CompactArrayWriterFile compactout(outputfilename,3);
	
	if ( ! rc )
		std::cerr << "[V] not storing reverse complements" << std::endl;

	// forward mapping table		
	libmaus::autoarray::AutoArray<uint8_t> ftable(256,false);
	// reverse complement mapping table
	libmaus::autoarray::AutoArray<uint8_t> rtable(256,false);
	// rc mapping for mapped symbols
	libmaus::autoarray::AutoArray<uint8_t> ctable(256,false);
	
	std::fill(ftable.begin(),ftable.end(),5);
	std::fill(rtable.begin(),rtable.end(),5);
	std::fill(ctable.begin(),ctable.end(),5);
	ftable['a'] = ftable['A'] = rtable['t'] = rtable['T'] = 1;
	ftable['c'] = ftable['C'] = rtable['g'] = rtable['G'] = 2;
	ftable['g'] = ftable['G'] = rtable['c'] = rtable['C'] = 3;
	ftable['t'] = ftable['T'] = rtable['a'] = rtable['A'] = 4;
	uint64_t insize = 0;
	
	ctable[1] = 4; // A->T
	ctable[2] = 3; // C->G
	ctable[3] = 2; // G->C
	ctable[4] = 1; // T->A
			
	for ( uint64_t i = 0; i < arginfo.restargs.size() && insize < limit; ++i )
	{
		std::string const fn = arginfo.stringRestArg(i);
		libmaus::aio::CheckedInputStream CIS(fn);
		libmaus::lz::BufferedGzipStream::unique_ptr_type BGS;
		std::istream * istr = 0;
		if ( gz )
		{
			libmaus::lz::BufferedGzipStream::unique_ptr_type tBGS(
				new libmaus::lz::BufferedGzipStream(CIS));
			BGS = UNIQUE_PTR_MOVE(tBGS);
			istr = BGS.get();
		}
		else
		{
			istr = &CIS;			
		}
		libmaus::fastx::StreamFastAReaderWrapper fain(*istr);
		libmaus::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		
		while ( fain.getNextPatternUnlocked(pattern) )
		{
			if ( verbose )
				std::cerr << (i+1) << " " << stripAfterDot(basename(fn)) << " " << pattern.sid << "...";
		
			// map symbols
			for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
				pattern.spattern[j] = ftable[static_cast<uint8_t>(pattern.spattern[j])];

			// write
			compactout.write(pattern.spattern.c_str(),pattern.spattern.size());
			compactout.put(0);
	
			if ( rc )
			{
				// reverse complement
				std::reverse(pattern.spattern.begin(),pattern.spattern.end());
				for ( uint64_t j = 0; j < pattern.spattern.size(); ++j )
					pattern.spattern[j] = ctable[static_cast<uint8_t>(pattern.spattern[j])];

				// write
				compactout.write(pattern.spattern.c_str(),pattern.spattern.size());
				compactout.put(0);
			}
						
			insize += pattern.spattern.size()+1;
			
			if ( verbose )
				std::cerr << "done, input size " << formatBytes(pattern.spattern.size()+1) << " acc " << formatBytes(insize) << std::endl;
		}
	}
	
	std::cerr << "Done, total input size " << insize << std::endl;
	
	compactout.flush();
	
	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);
		return fagzToCompact(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
