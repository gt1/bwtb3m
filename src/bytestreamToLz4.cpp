/**
    bwtb3m
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

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
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <libmaus2/lz/Lz4CompressStream.hpp>
#include <libmaus2/util/ArgInfo.hpp>

unsigned int getDefaultBlockSize()
{
	return 64*1024;
}

void bytestreamToLz4(libmaus2::util::ArgInfo const & arginfo)
{
	uint64_t const blocksize = arginfo.getValueUnsignedNumeric<uint64_t>("blocksize",getDefaultBlockSize());
	libmaus2::autoarray::AutoArray<char> B(8*1024,false);
	libmaus2::lz::Lz4CompressStream lzout(std::cout,blocksize);
	while ( std::cin )
	{
		std::cin.read(B.begin(),B.size());
		uint64_t const r = std::cin.gcount();
		if ( r )
			lzout.write(B.begin(),r);
	}

	lzout.writeIndex();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( arginfo.helpRequested() )
		{
			::libmaus2::exception::LibMausException se;

			std::ostream & str = se.getStream();

			str << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			str << std::endl;
			str << "usage: " << arginfo.progname << " [options]" << std::endl;
			str << std::endl;
			str << "options:" << std::endl;
			str << "blocksize=[" << getDefaultBlockSize() << "] block size" << std::endl;
			str << std::endl;
			str << "bytestreamToLz4 reads a byte stream from standard input and produces an lz4 compressed stream for bwtb3m on standard output." << std::endl;

			se.finish();
			throw se;
		}


		bytestreamToLz4(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
