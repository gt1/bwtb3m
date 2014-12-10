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
#include <libmaus/aio/PosixFdInputStream.hpp>
#include <libmaus/bitio/CompactArray.hpp>
#include <libmaus/bitio/CompactArrayWriterFile.hpp>
#include <libmaus/lz/BufferedGzipStream.hpp>
#include <libmaus/util/ArgInfo.hpp>

int digitsToCompact(libmaus::util::ArgInfo const & arginfo)
{
	// is input file gzipped?
	bool const gz = arginfo.getValue<unsigned int>("gz",0);
	// add terminator symbol?
	bool const addterm = arginfo.getValue<unsigned int>("term",0);
	// name of output file
	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
	libmaus::autoarray::AutoArray<char> B(8*1024,false);
	libmaus::bitio::CompactArrayWriterFile compactout(outputfilename,4);
	
	// error table
	libmaus::autoarray::AutoArray<uint8_t> etable(256,false);
	std::fill(etable.begin(),etable.end(),1);
	for ( int i = '0'; i <= '9'; ++i )
		etable[i] = 0;
	uint8_t const termadd = addterm ? 1 : 0;
	
	libmaus::lz::BufferedGzipStream::unique_ptr_type BGS;
	std::istream * istr = 0;
	if ( gz )
	{
		libmaus::lz::BufferedGzipStream::unique_ptr_type tBGS(
			new libmaus::lz::BufferedGzipStream(std::cin));
		BGS = UNIQUE_PTR_MOVE(tBGS);
		istr = BGS.get();
	}
	else
	{
		istr = &std::cin;
	}
	
	while ( *istr )
	{
		istr->read(B.begin(),B.size());
		uint64_t const num = istr->gcount();
		
		uint8_t err = 0;
		for ( uint64_t i = 0; i < num; ++i )
		{
			err = err | etable[B[i]];
			B[i] = (B[i] - '0')+termadd;
		}
			
		if ( err )
		{
			libmaus::exception::LibMausException lme;
			lme.getStream() << "Input file contains non decimal digit symbols." << std::endl;
			lme.finish();
			throw lme;
		}
	
		// write
		compactout.write(B.begin(),num);
	}
	
	// add terminator if requested
	if ( addterm )
	{
		uint8_t c = 0;
		compactout.write(&c,1);
	}
	
	compactout.flush();
	
	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus::util::ArgInfo const arginfo(argc,argv);
		int const r = digitsToCompact(arginfo);

		#if 0
		// read array and output data
		std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
		libmaus::aio::PosixFdInputStream PFIS(outputfilename);
		libmaus::bitio::CompactArray C(PFIS);
		for ( uint64_t i = 0; i < C.n; ++i )
			std::cout.put(C[i]+'0');
		#endif
		
		return r;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
