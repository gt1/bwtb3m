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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/aio/SynchronousGenericOutput.hpp>
#include <libmaus2/math/ilog.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		uint64_t const sub = arg.uniqueArgPresent("s") ? arg.getParsedArg<uint64_t>("s") : 1;

		assert ( (1ull << libmaus2::math::ilog(sub)) == sub );

		uint64_t const submask = sub-1;

		libmaus2::aio::SynchronousGenericInput<uint64_t> SGI(std::cin,8*1024);
		uint64_t inrate = 0;
		bool const inrateok = SGI.getNext(inrate);
		assert ( inrateok );
		uint64_t incnt = 0;
		bool const incntok = SGI.getNext(incnt);
		assert ( incntok );

		libmaus2::aio::SynchronousGenericOutput<uint64_t> SGO(std::cout,8*1024);
		uint64_t outrate = inrate * sub;
		SGO.put(outrate);
		uint64_t outcnt = (incnt + sub - 1)/sub;
		SGO.put(outcnt);

		for ( uint64_t i = 0; i < incnt; ++i )
		{
			uint64_t v = 0;
			bool const vok = SGI.getNext(v);
			assert ( vok );

			if ( !(i & submask) )
			{
				SGO.put(v);
			}
		}

		SGO.flush();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
