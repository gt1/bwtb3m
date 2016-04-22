/**
    bwtb3m
    Copyright (C) 2009-2016 German Tischler
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

#include <libmaus2/suffixsort/bwtb3m/BwtMergeSort.hpp>

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::timing::RealTimeClock rtc; rtc.start();

		if ( arginfo.helpRequested() || ! arginfo.restargs.size() )
		{
			::libmaus2::exception::LibMausException se;

			std::ostream & str = se.getStream();

			str << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			str << std::endl;
			str << "usage: " << arginfo.progname << " [options] <inputfile>" << std::endl;
			str << std::endl;
			str << "options:" << std::endl;
			str << "inputtype=[<bytestream>] (bytestream,compactstream,pac,pacterm,lz4,utf-8)" << std::endl;
			str << "outputfilename=[<"<< arginfo.getDefaultTmpFileName()+".bwt" << ">] (name of output .bwt file)" << std::endl;
			str << "sasamplingrate=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultSaSamplingRate() << "] sampling rate for sampled suffix array"<< std::endl;
			str << "isasamplingrate=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultIsaSamplingRate() << "] sampling rate for sampled inverse suffix array"<< std::endl;
			str << "mem=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultMem() << "] memory target" << std::endl;
			#if defined(_OPENMP)
			str << "numthreads=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultNumThreads() << "] number of threads" << std::endl;
			#endif
			str << "bwtonly=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultBWTOnly() << "] compute BWT only (no sampled suffix array and reverse)" << std::endl;
			str << std::string("tmpprefix=[") + arginfo.getDefaultTmpFileName() + "] (prefix for tmp files)" << std::endl;
			str << "sparsetmpprefix=[tmpprefix] (prefix for sparse gap tmp files)" << std::endl;
			str << "copyinputtomemory=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultCopyInputToMemory() << "] (copy input file to memory)" << std::endl;
			str << "largelcpthres=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultLargeLCPThres() << "] (large LCP value threshold)" << std::endl;
			str << "verbose=[" << libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultVerbose() << "] (verbosity level)" << std::endl;

			se.finish();
			throw se;
		}

		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions options(arginfo);
		libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt(options,&std::cerr);

		std::cerr << "[M] " << libmaus2::util::MemUsage() << " runtime " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
