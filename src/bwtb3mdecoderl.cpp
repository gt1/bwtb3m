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
#include <libmaus/huffman/RLDecoder.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <iostream>

void bwtb3mdecode(::libmaus::util::ArgInfo const & arginfo)
{
	std::string const infn = arginfo.getRestArg<std::string>(0);
	libmaus::huffman::RLDecoder dec(std::vector<std::string>(1,infn));
	std::string const inputtype = arginfo.getUnparsedValue("inputtype","bytestream");
	
	if ( inputtype == "bytestream" )
	{
		std::pair<int64_t,uint64_t> P;
		libmaus::aio::SynchronousGenericOutput<uint8_t> SGO(std::cout,64*1024);
		while ( (P = dec.decodeRun()).first >= 0 )
		{
			if ( P.first >= 256 )
			{
				::libmaus::exception::LibMausException se;
				se.getStream() << "unsupported character value " << P.first << " for inputtype=bytestream" << std::endl;
				se.finish();
				throw se;
			}
			for ( uint64_t i = 0; i < P.second; ++i )
				SGO.put(P.first);
		}
		
		SGO.flush();
	}
	else
	{
		::libmaus::exception::LibMausException se;
		se.getStream() << "input type " << inputtype << " is currently not supported for bwtb3mdecode" << std::endl;
		se.finish();
		throw se;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		
		if ( arginfo.helpRequested() || arginfo.restargs.size() < 1 )
		{
			::libmaus::exception::LibMausException se;
			std::ostream & str = se.getStream();
			str << "usage: " << argv[0] << " <in.bwt>" << std::endl;
			str << std::endl;
			
			se.finish();
			throw se;
		}
		
		bwtb3mdecode(arginfo);

		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
