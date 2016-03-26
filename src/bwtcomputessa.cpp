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
#include <libmaus2/suffixsort/bwtb3m/BwtComputeSSA.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( arginfo.getNumRestArgs() < 1 )
		{
			std::cerr << "usage: " << argv[0] << " <bwt>" << std::endl;
			return EXIT_FAILURE;
		}

		#if defined(_OPENMP)
		unsigned int const maxthreads = omp_get_max_threads();
		#else
		unsigned int const maxthreads = 1;
		#endif

		std::string bwt = arginfo.getUnparsedRestArg(0);
		uint64_t const sasamplingrate = arginfo.getValueUnsignedNumeric<uint64_t>("sasamplingrate",32);
		uint64_t const isasamplingrate = arginfo.getValueUnsignedNumeric<uint64_t>("isasamplingrate",32);
		std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpprefix",arginfo.getDefaultTmpFileName());
		bool const copyinputtomemory = arginfo.getValue<uint64_t>("copyinputtomemory",false);
		// get number of threads
		unsigned int const numthreads = arginfo.getValue<unsigned int>("threads", maxthreads);
		uint64_t const maxsortmem = arginfo.getValueUnsignedNumeric<uint64_t>("maxsortmem",2*1024ull*1024ull*1024ull);
		uint64_t const maxtmpfiles = arginfo.getValueUnsignedNumeric<uint64_t>("maxtmpfiles",1024);
		std::ostream * logstr = &(std::cerr);
		std::string const ref_isa_fn = arginfo.getUnparsedValue("ref_isa",std::string());
		std::string const ref_sa_fn = arginfo.getUnparsedValue("ref_sa",std::string());

		libmaus2::suffixsort::bwtb3m::BwtComputeSSA::computeSSA(bwt,sasamplingrate,isasamplingrate,tmpfilenamebase,copyinputtomemory,numthreads,maxsortmem,maxtmpfiles,logstr,ref_isa_fn,ref_sa_fn);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
