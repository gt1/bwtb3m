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
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSortResult.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSortOptions.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/suffixsort/BwtMergeBlockSortRequest.hpp>

template<typename io_type>
int checkBwt(libmaus2::util::ArgParser const & arg)
{
	std::string const bwtname = arg[0];
	std::string const textname = arg[1];
	std::string const prefix = libmaus2::util::OutputFileNameTools::clipOff(bwtname,".bwt");
	// std::cerr << "prefix " << prefix << std::endl;
	std::string const preisa = prefix + ".preisa";
	std::string const checkinfofn = prefix + ".preisa.checkinfo";
	//std::string const textname = libmaus2::util::OutputFileNameTools::clipOff(prefix,".bwt");
	std::string const tmpprefix = libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	libmaus2::suffixsort::bwtb3m::BwtMergeSortResult res;
	res.textfn = textname;
	res.bwtfn = bwtname;
	res.preisafn = preisa;
	res.histfn = prefix + ".hist";
	res.metafn = prefix + ".meta";
	
	if ( 
		libmaus2::util::GetFileSize::fileExists(prefix + ".hwt") 
		&&
		!libmaus2::util::GetFileSize::isOlder(prefix+".hwt",bwtname)
	)
	{
		res.hwtfn = prefix + ".hwt";
	}

	uint64_t const numthreads = libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	//uint64_t const numthreads = 3;
	uint64_t const n = libmaus2::huffman::RLDecoder::getLength(bwtname,numthreads);
	uint64_t const o = (n + numthreads - 1)/numthreads;
	
	struct Info
	{
		uint64_t pq;
		uint64_t p;
		uint64_t r;
		
		static uint64_t getUndef()
		{
			return std::numeric_limits<uint64_t>::max();
		}
		
		Info(uint64_t const rpq = getUndef(), uint64_t const rp = getUndef(), uint64_t const rr = getUndef())
		: pq(rpq), p(rp), r(rr) {}
		
		bool valid() const
		{
			return p != getUndef();
		}
		
		std::ostream & serialise(std::ostream & out) const
		{
			libmaus2::util::NumberSerialisation::serialiseNumber(out,pq);
			libmaus2::util::NumberSerialisation::serialiseNumber(out,p);
			libmaus2::util::NumberSerialisation::serialiseNumber(out,r);
			return out;
		}
		
		std::istream & deserialise(std::istream & in)
		{
			pq = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
			p = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
			r = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
			return in;
		}
		
		Info(std::istream & in)
		{
			deserialise(in);
		}
		
		static void serialiseVector(std::string const & fn, std::vector<Info> const & V)
		{
			libmaus2::aio::OutputStreamInstance OSI(fn);
			for ( uint64_t i = 0; i < V.size(); ++i )
				V[i].serialise(OSI);
		}
		
		static std::vector<Info> deserialiseVector(std::string const & fn)
		{
			libmaus2::aio::InputStreamInstance ISI(fn);
			std::vector<Info> V;
			while ( ISI.peek() != std::istream::traits_type::eof() )
				V.push_back(Info(ISI));
			return V;
		}
	};
	
		
	struct InfoPQComp
	{
		bool operator()(Info const & A, Info const & B) const
		{
			return A.pq < B.pq;
		}
	};
	
	{
		std::vector<Info> I(numthreads);
		for ( uint64_t i = 0; i < I.size(); ++i )
			I[i].pq = i*o;

		libmaus2::aio::InputStreamInstance ISI(preisa);
		libmaus2::aio::SynchronousGenericInput< std::pair<uint64_t,uint64_t> > SGI(ISI,4*1024);
		std::pair<uint64_t,uint64_t> P;
		while ( SGI.getNext(P) )
		{
			typename std::vector<Info>::iterator it = std::lower_bound(I.begin(),I.end(),Info(P.second),InfoPQComp());
			if ( (it == I.end()) || (it->pq != P.second) )
			{
				it -= 1;
			}

			assert ( P.second >= it->pq );
			
			if ( P.second < it->p )
			{
				it->p = P.second;
				it->r = P.first;
			}
		}

		uint64_t o = 0;
		for ( uint64_t i = 0; i < I.size(); ++i )
			if ( I[i].valid() )
				I[o++] = I[i];
		I.resize(o);
		
		Info::serialiseVector(checkinfofn,I);
	}
	
	std::vector<Info> I = Info::deserialiseVector(checkinfofn);
	
	for ( uint64_t i = 0; i < I.size(); ++i )
		std::cerr << "(" << I[i].pq << "," << I[i].p << "," << I[i].r << ")" << std::endl;

	libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PLF(res.loadLF(tmpprefix,numthreads));
	libmaus2::lf::ImpCompactHuffmanWaveletLF & LF = *PLF;

	typedef typename  io_type::circular_reverse_wrapper crw_type;
	uint64_t volatile gc = 0;
	uint64_t volatile gok = 1;
	libmaus2::parallel::PosixSpinLock gclock;
	uint64_t const verb = 1024*1024ull;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1)
	#endif	
	for ( uint64_t t = 0; t < I.size(); ++t )
	{
		uint64_t t1 = (t + 1) % I.size();
		
		std::pair<uint64_t,uint64_t> const PL(I[t].p,I[t].r);
		std::pair<uint64_t,uint64_t> const PH(I[t1].p,I[t1].r);
				
		crw_type C(textname,PH.first);
		
		uint64_t const bu = PH.first ? PH.first : n;
		uint64_t const bl = PL.first;
		uint64_t const b = bu-bl;
		bool tok = true;
		
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr 
				<< "[V] t=" << t << " PL=(" << PL.first << "," << PL.second 
					<< ") PH=(" << PH.first << "," << PH.second << ")" << std::endl;
		}

		uint64_t r = PH.second;
		uint64_t c = 0;
		
		for ( uint64_t i = 0; i < b; ++i )
		{
			std::pair<int64_t,uint64_t> L = LF.extendedLF(r);
			bool const ok = ( L.first == C.get() );
			if ( ! ok )
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[E] failure in thread " << t << " rank " << r << std::endl;
				tok = false;
			}
			r = L.second;
			
			if ( ++c == verb )
			{
				libmaus2::parallel::ScopePosixSpinLock slock(gclock);
				gc += c;
				c = 0;
				
				{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] " << gc << "/" << n << " " << static_cast<double>(gc)/static_cast<double>(n) << std::endl;
				}
			}
		}

		if ( c )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(gclock);
			gc += c;
			c = 0;

			{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "[V] " << gc << "/" << n << " " << static_cast<double>(gc)/static_cast<double>(n) << std::endl;
			}
		}
		
		gclock.lock();
		if ( ! tok )
			gok = 0;
		gclock.unlock();
	}
	
	std::cerr << "[V] gok=" << gok << std::endl;

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);
		std::string const inputtype = arg.uniqueArgPresent("i") ? arg["i"] : "bytestream";
		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type itype = 
			libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::parseInputType(inputtype);

		switch ( itype )
		{
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_bytestream:
				return checkBwt<libmaus2::suffixsort::ByteInputTypes>(arg);
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_compactstream:
				return checkBwt<libmaus2::suffixsort::CompactInputTypes>(arg);
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_pac:
				return checkBwt<libmaus2::suffixsort::PacInputTypes>(arg);
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_pacterm:
				return checkBwt<libmaus2::suffixsort::PacTermInputTypes>(arg);
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_lz4:
				return checkBwt<libmaus2::suffixsort::Lz4InputTypes>(arg);
			case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_utf_8:
				return checkBwt<libmaus2::suffixsort::Utf8InputTypes>(arg);
			default:
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt: unknown/unsupported input type" << std::endl;
				lme.finish();
				throw lme;
			}
			break;
		}

	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
