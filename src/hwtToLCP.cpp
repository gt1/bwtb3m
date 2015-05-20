/**
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
**/
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

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo arginfo(argc,argv);
		::libmaus2::timing::RealTimeClock rtc;
		std::string const tmpfilenamebase = arginfo.getDefaultTmpFileName();
		libmaus2::util::TempFileNameGenerator tmpgen(tmpfilenamebase+"_lcpdir",3);
		libmaus2::util::FileTempFileContainer tmpcont(tmpgen);
		
		std::string const namebase = arginfo.getRestArg<std::string>(0);
		bool const showavg = arginfo.getValue<bool>("showavg",false);
		bool const checklcp = arginfo.getValue<bool>("checklcp",false);
		bool const creatermmtree = arginfo.getValue<bool>("rmmtree",true);
		std::string const hwtname = namebase + ".hwt";
		std::string const isaname = namebase + ".isa";
		std::string const saname = namebase + ".sa";
		std::string const lcpname = namebase + ".lcp";
		std::string const ulcpname = namebase + ".ulcp";
		std::string const rmmname = namebase + ".rmm";

		std::cerr << "[V] Loading index...";
		::libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PIHWLF =
			UNIQUE_PTR_MOVE(::libmaus2::lf::ImpCompactHuffmanWaveletLF::load(hwtname));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF const & IHWLF = *PIHWLF;
		std::cerr << "done." << std::endl;

		// compute LCP array
		std::cerr << "[V] Computing LCP array...";
		rtc.start();
		::libmaus2::lcp::WaveletLCPResult::unique_ptr_type LCP = 
			::libmaus2::lcp::WaveletLCP::computeLCP(&IHWLF,false /* zero symbols are not distinct */);
		std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
		
		std::cerr << "[V] Writing LCP array...";
		::libmaus2::aio::CheckedOutputStream ulcpCOS(ulcpname);
		LCP->serialise(ulcpCOS);
		ulcpCOS.flush();
		ulcpCOS.close();
		std::cerr << "done." << std::endl;

		typedef ::libmaus2::lcp::SuccinctLCP<
			::libmaus2::lf::ImpCompactHuffmanWaveletLF,
			::libmaus2::fm::SimpleSampledSA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF >,
			::libmaus2::fm::SampledISA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF >
		> succinct_lcp_type;

		std::cerr << "[V] Loading inverse sampled suffix array...";
		::libmaus2::aio::CheckedInputStream isain(isaname);
		::libmaus2::fm::SampledISA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF > ISA(&IHWLF,isain);
		std::cerr << "done." << std::endl;
				
		std::cerr << "[V] Serialising succinct PLCP array...";
		rtc.start();
		::libmaus2::aio::CheckedOutputStream lcpCOS(lcpname);
		succinct_lcp_type::writeSuccinctLCP(IHWLF,ISA,*LCP,lcpCOS,tmpcont,false /* verbose */);
		lcpCOS.flush();
		lcpCOS.close();
		std::cerr << "done, time " << rtc.getElapsedSeconds()  << std::endl;

		#if defined(_OPENMP)		
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif

		// create rmm tree
		if ( creatermmtree )
		{
			::libmaus2::lcp::WaveletLCPResult const & LCPref = *LCP;
			typedef libmaus2::rmq::RMMTree< ::libmaus2::lcp::WaveletLCPResult ,3> rmm_tree_type;
			rmm_tree_type rmmtree(LCPref,PIHWLF->n);
			
			libmaus2::aio::OutputStream::unique_ptr_type Prmmout(libmaus2::aio::OutputStreamFactoryContainer::constructUnique(rmmname));
			rmmtree.serialise(*Prmmout);
		}

		if ( showavg )
		{
			std::cerr << "[V] Loading sa...";
			::libmaus2::aio::CheckedInputStream sain(saname);
			::libmaus2::fm::SimpleSampledSA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF > SA(&IHWLF,sain);
			std::cerr << "done." << std::endl;

			std::cerr << "[V] Loading succinct LCP...";
			::libmaus2::aio::CheckedInputStream lcpCIS(lcpname);
			succinct_lcp_type const SLCP(lcpCIS,SA);
			std::cerr << "done." << std::endl;
		
			uint64_t const symsperblock = (IHWLF.n + numthreads - 1)/numthreads;
			::libmaus2::parallel::OMPLock lock;
			uint64_t gacc = 0;
			::libmaus2::util::Histogram ghist;
			uint64_t const lcpthres = 1024;

			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif		
			for ( int64_t b = 0; b < static_cast<int64_t>(numthreads); ++b )
			{
				uint64_t const low = std::min(b*symsperblock,IHWLF.n);
				uint64_t const high = std::min(low+symsperblock,IHWLF.n);
				
				uint64_t acc = 0;
				::libmaus2::util::Histogram lhist;
				for ( uint64_t i = low; i < high; ++i )
				{
					uint64_t const lcp = (SLCP)[i];
					uint64_t const caplcp = std::min(lcpthres,lcp);
					acc += lcp; 
					lhist(caplcp);
				}
				
				lock.lock();
				gacc += acc;
				ghist.merge(lhist);
				lock.unlock();
			}
			
			std::cerr << "[V] avg LCP is " << static_cast<double>(gacc) / IHWLF.n << std::endl;
			ghist.print(std::cerr);
		}
			
		if ( checklcp )
		{
			std::cerr << "[V] Loading sa...";
			::libmaus2::aio::CheckedInputStream sain(saname);
			::libmaus2::fm::SimpleSampledSA< ::libmaus2::lf::ImpCompactHuffmanWaveletLF > SA(&IHWLF,sain);
			std::cerr << "done." << std::endl;

			std::cerr << "[V] Loading succinct LCP...";
			::libmaus2::aio::CheckedInputStream lcpCIS(lcpname);
			succinct_lcp_type const SLCP(lcpCIS,SA);
			std::cerr << "done." << std::endl;

			std::cerr << "[V] Checking succinct LCP against non succinct version...";			
			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t r = 0; r < static_cast<int64_t>(IHWLF.n); ++r )
				assert ( (*LCP)[r] == SLCP[r] );
			std::cerr << "done." << std::endl;
			
			uint64_t const checkranks = IHWLF.n ? (IHWLF.n-1) : 0;
			uint64_t const checkpacksize = (checkranks + numthreads-1)/numthreads;
			uint64_t const numcheckpacks = (checkranks + checkpacksize-1)/checkpacksize;
			::libmaus2::parallel::SynchronousCounter<uint64_t> cnt(0);
			
			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t packid = 0; packid < static_cast<int64_t>(numcheckpacks); ++packid )
			{
				uint64_t const checklow = std::min(packid * checkpacksize,checkranks);
				uint64_t const checkhigh = std::min(checklow + checkpacksize,checkranks);
				
				for ( uint64_t r = checklow+1; r < (checkhigh+1); ++r )
				{
					uint64_t l = SLCP[r];
					uint64_t r1 = ISA[(SA[r-1] + l + 1)%IHWLF.n];
					uint64_t r0 = ISA[(SA[r-0] + l + 1)%IHWLF.n];
					
					assert ( IHWLF[r0] != IHWLF[r1] );
					
					while ( l-- )
					{
						r0 = IHWLF(r0);
						r1 = IHWLF(r1);
						assert ( IHWLF[r0] == IHWLF[r1] );
					}
					
					uint64_t const lcnt = ++cnt;
					if ( lcnt % (1024*1024) == 0 || lcnt == checkranks )
						std::cerr << static_cast<double>(lcnt)/checkranks << std::endl;			
				}
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
