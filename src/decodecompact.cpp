/*
    bwtb3m
    Copyright (C) 2016 German Tischler

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
*/
#include <libmaus2/bitio/CompactDecoderBuffer.hpp>
#include <libmaus2/util/ArgParser.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		libmaus2::autoarray::AutoArray<char> B(8*1024,false);
		
		for ( uint64_t i = 0; i < arg.size(); ++i )
		{
			libmaus2::bitio::CompactDecoderWrapper CDW(arg[i]);
			
			while ( CDW )
			{
				CDW.read(B.begin(),B.size());
				std::cout.write(B.begin(),CDW.gcount());
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}