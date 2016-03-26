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
#include <config.h>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/bitio/CompactArray.hpp>
#include <libmaus2/bitio/CompactArrayWriterFile.hpp>
#include <libmaus2/lz/BufferedGzipStream.hpp>
#include <libmaus2/util/ArgInfo.hpp>

int digitsToCompact(libmaus2::util::ArgInfo const & arginfo)
{
	// is input file gzipped?
	bool const gz = arginfo.getValue<unsigned int>("gz",0);
	// add terminator symbol?
	bool const addterm = arginfo.getValue<unsigned int>("term",0);
	// name of output file
	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
	libmaus2::autoarray::AutoArray<char> B(8*1024,false);
	libmaus2::bitio::CompactArrayWriterFile compactout(outputfilename,4);

	// error table
	libmaus2::autoarray::AutoArray<uint8_t> etable(256,false);
	std::fill(etable.begin(),etable.end(),1);
	for ( int i = '0'; i <= '9'; ++i )
		etable[i] = 0;
	uint8_t const termadd = addterm ? 1 : 0;

	libmaus2::lz::BufferedGzipStream::unique_ptr_type BGS;
	std::istream * istr = 0;
	if ( gz )
	{
		libmaus2::lz::BufferedGzipStream::unique_ptr_type tBGS(
			new libmaus2::lz::BufferedGzipStream(std::cin));
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
			libmaus2::exception::LibMausException lme;
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
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( arginfo.helpRequested() )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			std::cerr << std::endl;
			std::cerr << "usage: " << arginfo.progname << " [options] < input" << std::endl;
			std::cerr << std::endl;
			std::cerr << "options:" << std::endl;
			std::cerr << "gz=[0|1] (input file is plain (0)/gzip compressed (1), default is uncompressed)" << std::endl;
			std::cerr << "outputfilename=<output.compact> (name of output file, default is output.compact)" << std::endl;
			std::cerr << "term=[0|1] (0: map digits to symbols 0-9, 1: map digits to symbols 1-10 and append 0 at the end)" << std::endl;

			return EXIT_SUCCESS;
		}
		else
		{
			int const r = digitsToCompact(arginfo);

			#if 0
			// read array and output data
			std::string const outputfilename = arginfo.getUnparsedValue("outputfilename","output.compact");
			libmaus2::aio::InputStreamInstance PFIS(outputfilename);
			libmaus2::bitio::CompactArray C(PFIS);
			for ( uint64_t i = 0; i < C.n; ++i )
				std::cout.put(C[i]+'0');
			#endif

			return r;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
