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
#include <iostream>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/fm/SampledISA.hpp>
#include <libmaus2/lf/ImpCompactHuffmanWaveletLF.hpp>
#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/huffman/RLDecoder.hpp>
#include <libmaus2/sorting/InPlaceParallelSort.hpp>

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo arginfo(argc,argv);
		::libmaus2::timing::RealTimeClock rtc;
		std::string const tmpfilenamebase = arginfo.getDefaultTmpFileName();

		uint64_t const numthreads = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
		std::string const namebase = arginfo.getUnparsedRestArg(0);
		std::string const hwtname = namebase + ".hwt";
		std::string const preisaname = namebase + ".preisa";
		std::string const metaname = namebase + ".preisa.meta";
		std::string const isaname = namebase + ".isa";
		// std::string const saname = namebase + ".sa";
		
		#if 0
		uint64_t rate = 0;
		{
			libmaus2::aio::InputStreamInstance ISI(metaname);
			rate = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		}
		#endif
		
		uint64_t n = libmaus2::huffman::RLDecoder::getLength(namebase + ".bwt",numthreads);
		
		uint64_t const preisafs = libmaus2::util::GetFileSize::getFileSize(preisaname);
		assert ( preisafs % (2*sizeof(uint64_t)) == 0 );
		uint64_t const numpreisasamples = preisafs / (2*sizeof(uint64_t));
		
		uint64_t const samplesperthread = (numpreisasamples + numthreads - 1)/numthreads;
		uint64_t const numpackages = (numpreisasamples + samplesperthread - 1)/samplesperthread;
		
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > Asamples(numpreisasamples,false);
		
		{
			libmaus2::aio::InputStreamInstance ISI(preisaname);
			libmaus2::aio::SynchronousGenericInput<uint64_t> SGI(ISI,4096);
			
			uint64_t r, p;
			uint64_t i = 0;
			while ( SGI.getNext(r) )
			{
				bool const ok = SGI.getNext(p);
				assert ( ok );
				
				Asamples[i++] = std::pair<uint64_t,uint64_t>(p,r);
			}
		}
		
		libmaus2::sorting::InPlaceParallelSort::inplacesort(Asamples.begin(),Asamples.end(),numthreads);
		
		#if 0
		for ( uint64_t i = 0; i < numpreisasamples; ++i )
			std::cerr << Asamples[i].first << "\t" << Asamples[i].second << std::endl;
		#endif
		
		std::cerr << "[V] n=" << n << std::endl;
		
		uint64_t isasamplingrate = arginfo.getValueUnsignedNumeric<uint64_t>("isasamplingrate",64);
		
		unsigned int bcnt = 0;
		for ( uint64_t tt = isasamplingrate; tt; tt>>= 1)
			if ( tt&1 )
				bcnt++;
				
		assert ( bcnt == 1 );
		
		uint64_t const isasamplingmask = isasamplingrate-1;
		
		std::cerr << "[V] using output sampling rate " << isasamplingrate << std::endl;

		std::cerr << "[V] Loading index...";
		::libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PIHWLF =
			UNIQUE_PTR_MOVE(::libmaus2::lf::ImpCompactHuffmanWaveletLF::load(hwtname,numthreads));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF const & IHWLF = *PIHWLF;
		std::cerr << "done." << std::endl;
		
		uint64_t const numisa = (n + isasamplingrate - 1)/isasamplingrate;
		libmaus2::autoarray::AutoArray<uint64_t> ISA(numisa,false);
		std::fill(ISA.begin(),ISA.end(),std::numeric_limits<uint64_t>::max());

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t t = 0; t < numpackages; ++t )
		{
			uint64_t const thispackage = t;
			uint64_t const prevpackage = (t + numpackages - 1)%numpackages;
			
			uint64_t const thispackageid = thispackage * samplesperthread;
			uint64_t const prevpackageid = prevpackage * samplesperthread;
			
			uint64_t const thisrank = Asamples[thispackageid].second;
			uint64_t const thispos = Asamples[thispackageid].first;
			uint64_t const prevrank = Asamples[prevpackageid].second;
			uint64_t const prevpos = Asamples[prevpackageid].first;
			
			//std::cerr << "thispos=" << thispos << " prevpos=" << prevpos << std::endl;
			
			uint64_t todo;
			
			if ( thispos > prevpos )
				todo = thispos-prevpos;
			else if ( thispos < prevpos )
			{
				assert ( thispos == 0 );
				todo = n - prevpos;
			}
			else
			{
				assert ( thispos == prevpos );
				assert ( thispos == 0 );
				todo = n;
			}
			
			//std::cerr << "todo=" << todo << std::endl;
			
			uint64_t curpos = thispos;
			uint64_t currank = thisrank;
			
			while ( todo-- )
			{
				if ( (curpos & isasamplingmask) == 0 )
				{
					ISA [ curpos / isasamplingrate ] = currank;
					// std::cerr << "curpos=" << curpos << " currank=" << currank << std::endl;
				}
				
				curpos = (curpos + n - 1) % n;
				currank = (*PIHWLF)(currank);
			}
		}
		
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < ISA.size(); ++i )
			assert ( ISA[i] != std::numeric_limits<uint64_t>::max() );

		{
			libmaus2::aio::OutputStreamInstance OSI(isaname);
			libmaus2::fm::SampledISABase::writeSampledInverseSuffixArray(OSI,ISA.begin(),n,isasamplingrate,1);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
