/**
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
**/
#include <libmaus2/lcp/PLCPBitDecoder.hpp>
#include <iostream>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/fm/SampledSA.hpp>
#include <libmaus2/fm/SampledISA.hpp>
#include <libmaus2/lcp/SuccinctLCP.hpp>
#include <libmaus2/lcp/WaveletLCP.hpp>
#include <libmaus2/lf/ImpCompactHuffmanWaveletLF.hpp>
#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/rmq/RMMTree.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/FileTempFileContainer.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/wavelet/ImpHuffmanWaveletTree.hpp>
#include <libmaus2/parallel/NumCpus.hpp>


int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo arginfo(argc,argv);
		::libmaus2::timing::RealTimeClock rtc;
		std::string const tmpfilenamebase = arginfo.getDefaultTmpFileName();
		libmaus2::util::TempFileNameGenerator tmpgen(tmpfilenamebase+"_lcpdir",3);
		libmaus2::util::FileTempFileContainer tmpcont(tmpgen);

		uint64_t const numthreads = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
		std::string const namebase = arginfo.getRestArg<std::string>(0);
		std::string const hwtname = namebase + ".hwt";
		std::string const isaname = namebase + ".isa";
		std::string const saname = namebase + ".sa";
		std::string const lcpbitname = namebase + ".lcpbit";

		std::cerr << "[V] Loading index...";
		::libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PIHWLF =
			UNIQUE_PTR_MOVE(::libmaus2::lf::ImpCompactHuffmanWaveletLF::load(hwtname,numthreads));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF const & IHWLF = *PIHWLF;
		std::cerr << "done." << std::endl;

		// compute LCP array
		std::cerr << "[V] Computing LCP array...";
		rtc.start();
		::libmaus2::lcp::WaveletLCPResult::unique_ptr_type LCP =
			::libmaus2::lcp::WaveletLCP::computeLCP(&IHWLF,numthreads,false /* zero symbols are not distinct */,&(std::cerr));
		std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;

		uint64_t const n = IHWLF.getN();
		libmaus2::lcp::PLCPBitDecoder SLCP(lcpbitname,n);

		std::cerr << "[V] Loading sa...";
		::libmaus2::aio::InputStreamInstance sain(saname);
		::libmaus2::fm::SimpleSampledSA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF > SA(&IHWLF,sain);
		std::cerr << "done." << std::endl;

		std::cerr << "[V] checking...";
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < n; ++i )
		{
			// std::cerr << (*LCP)[i] << " " << SLCP.get(SA[i]) << std::endl;
			assert ( (*LCP)[i] == SLCP.get(SA[i]) );
		}
		std::cerr << "done.\n";
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
