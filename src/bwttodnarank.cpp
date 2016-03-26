/*
    bwtb3m
    Copyright (C) 2015 German Tischler

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
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/rank/DNARank.hpp>
#include <libmaus2/parallel/NumCpus.hpp>

/**
 * convert a BWT over 4 symbols (0,1,2,3) to a serialised DNARank object
 **/
int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		unsigned int const numthreads = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
		std::string const fn = arg[0];
		std::string const outfn = libmaus2::util::OutputFileNameTools::clipOff(fn,".bwt") + ".dnarank";
		libmaus2::rank::DNARank::unique_ptr_type Prank(libmaus2::rank::DNARank::loadFromRunLength(fn,numthreads));
		libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(outfn));
		Prank->serialise(*OSI);
		OSI->flush();
		OSI.reset();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
