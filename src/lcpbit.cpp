/*
    libmaus2
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
#include <libmaus2/aio/CircularWrapper.hpp>
#include <libmaus2/aio/PairFileDecoder.hpp>
#include <libmaus2/bitio/BitVectorInput.hpp>
#include <libmaus2/bitio/BitVectorOutput.hpp>
#include <libmaus2/gamma/GammaFlaggedIntervalDecoder.hpp>
#include <libmaus2/gamma/GammaFlaggedIntervalEncoder.hpp>
#include <libmaus2/gamma/GammaFlaggedPartitionDecoder.hpp>
#include <libmaus2/gamma/GammaFlaggedPartitionEncoder.hpp>
#include <libmaus2/gamma/GammaPDDecoder.hpp>
#include <libmaus2/gamma/GammaPDEncoder.hpp>
#include <libmaus2/huffman/LFPhiPairDecoder.hpp>
#include <libmaus2/huffman/LFPhiPairEncoder.hpp>
#include <libmaus2/huffman/LFPhiPairLCPDecoder.hpp>
#include <libmaus2/huffman/LFPhiPairLCPEncoder.hpp>
#include <libmaus2/huffman/LFRankLCPDecoder.hpp>
#include <libmaus2/huffman/LFRankLCPEncoder.hpp>
#include <libmaus2/huffman/LFRankPosDecoder.hpp>
#include <libmaus2/huffman/LFRankPosEncoder.hpp>
#include <libmaus2/huffman/LFSetBitDecoder.hpp>
#include <libmaus2/huffman/LFSetBitEncoder.hpp>
#include <libmaus2/huffman/LFSymRankPosDecoder.hpp>
#include <libmaus2/huffman/LFSymRankPosEncoder.hpp>
#include <libmaus2/huffman/RLDecoder.hpp>
#include <libmaus2/huffman/RLEncoder.hpp>
#include <libmaus2/huffman/SymBitDecoder.hpp>
#include <libmaus2/huffman/SymBitEncoder.hpp>
#include <libmaus2/lcp/LCP.hpp>
#include <libmaus2/lz/XzInputStream.hpp>
#include <libmaus2/math/ilog.hpp>
#include <libmaus2/random/Random.hpp>
#include <libmaus2/sorting/PairFileSorting.hpp>
#include <libmaus2/sorting/ParallelExternalRadixSort.hpp>
#include <libmaus2/sorting/ParallelRunLengthRadixSort.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSortOptions.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtSelectSSA.hpp>
#include <libmaus2/suffixsort/BwtMergeBlockSortRequest.hpp>
#include <libmaus2/suffixsort/ByteInputTypes.hpp>
#include <libmaus2/suffixsort/divsufsort.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/NumberMapSerialisation.hpp>
#include <libmaus2/util/PrefixSums.hpp>

/**
 * linear array accessor through file
 **/
template<typename _stream_type>
struct LinearAccessor
{
	typedef _stream_type stream_type;
	typedef typename stream_type::traits_type::char_type char_type;

	stream_type & in;
	uint64_t p;

	LinearAccessor(stream_type & rin) : in(rin), p(0) {}

	char_type operator[](uint64_t i)
	{
		assert ( i >= p );
		if ( i > p )
			in.ignore(i-p);

		p = i;
		return in.peek();
	}
};

// get size of decompressed XZ file
uint64_t xzSize(std::istream & ISI)
{
	libmaus2::lz::XzInputStream xzin(ISI);
	uint64_t c = 0;
	while ( xzin )
	{
		xzin.ignore(64*10124);
		c += xzin.gcount();
	}

	return c;
}

// read XZ compressed file and append # as terminator
std::string readXzFile(std::string const & fn)
{
	libmaus2::aio::InputStreamInstance ISI(fn);
	ISI.clear();
	ISI.seekg(0,std::ios::beg);
	uint64_t const s = xzSize(ISI);
	libmaus2::autoarray::AutoArray<char> A(s+1,false);
	ISI.clear();
	ISI.seekg(0,std::ios::beg);
	libmaus2::lz::XzInputStream xzin(ISI);
	xzin.read(A.begin(),s);
	assert ( xzin.gcount() == static_cast<std::streamsize>(s) );
	int64_t minc = std::numeric_limits<int64_t>::max();
	for ( uint64_t i = 0; i < s; ++i )
		minc = std::min(minc,static_cast<int64_t>(A[i]));
	A[s] = '#';
	return std::string(A.begin(),A.end());
}

// read a plain text file and append 0
std::string readPlainFile(std::string const & fn)
{
	libmaus2::aio::InputStreamInstance ISI(fn);
	uint64_t const s = libmaus2::util::GetFileSize::getFileSize(ISI);
	ISI.clear();
	ISI.seekg(0);
	libmaus2::autoarray::AutoArray<char> A(s+1,false);
	ISI.read(A.begin(),s);
	assert ( static_cast<int64_t>(ISI.gcount()) == static_cast<int64_t>(s) );
	for ( uint64_t i = 0; i < s; ++i )
		assert ( A[s] != 0 );
	A[s] = 0;
	return std::string(A.begin(),A.end());
}

// merge queue files so they contain complete (non fractional) intervals only
void mergeQueueFiles(std::vector<std::string> & qfiletmp, uint64_t const numthreads, libmaus2::util::TempFileNameGenerator & tmpgen)
{
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < qfiletmp.size(); ++i )
	{
		bool const empty = libmaus2::gamma::GammaFlaggedIntervalDecoder::isEmpty(qfiletmp[i]);

		if ( empty )
		{
			libmaus2::aio::FileRemoval::removeFile(qfiletmp[i]);
			qfiletmp[i] = std::string();
		}
	}

	uint64_t o = 0;
	for ( uint64_t i = 0; i < qfiletmp.size(); ++i )
		if ( qfiletmp[i].size() )
			qfiletmp[o++] = qfiletmp[i];
	qfiletmp.resize(o);

	bool merging = true;

	while ( merging )
	{
		merging = false;

		libmaus2::autoarray::AutoArray < libmaus2::gamma::FlaggedInterval::interval_type > Atype(qfiletmp.size());

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < qfiletmp.size(); ++i )
		{
			libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type dec(new libmaus2::gamma::GammaFlaggedIntervalDecoder(std::vector<std::string>(1,qfiletmp[i]),0 /* offset */));
			libmaus2::gamma::FlaggedInterval QP;
			bool const ok = dec->getNext(QP);
			assert ( ok );
			Atype[i] = QP.type;
		}

		uint64_t j = 0;
		while ( j+1 < Atype.size() )
		{
			if ( libmaus2::gamma::FlaggedInterval::mergeable(Atype[j],Atype[j+1]) )
			{
				libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type dec0(new libmaus2::gamma::GammaFlaggedIntervalDecoder(std::vector<std::string>(1,qfiletmp[j]),0 /* offset */));
				libmaus2::gamma::FlaggedInterval QP0;
				bool const ok0 = dec0->getNext(QP0);
				assert ( ok0 );
				assert ( QP0.type == Atype[j] );
				dec0.reset();

				libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type dec1(new libmaus2::gamma::GammaFlaggedIntervalDecoder(std::vector<std::string>(1,qfiletmp[j+1]),0 /* offset */));
				libmaus2::gamma::FlaggedInterval QP1;
				bool const ok1 = dec1->getNext(QP1);
				assert ( ok1 );
				assert ( QP1.type == Atype[j+1] );
				dec1.reset();

				assert ( QP0.to == QP1.from );

				libmaus2::gamma::FlaggedInterval QPm;
				QPm.from = QP0.from;
				QPm.to = QP1.to;
				QPm.type = libmaus2::gamma::FlaggedInterval::merge(QP0.type,QP1.type);
				assert ( QP0.active == QP1.active );
				QPm.active = QP0.active;

				libmaus2::aio::FileRemoval::removeFile(qfiletmp[j]);
				libmaus2::aio::FileRemoval::removeFile(qfiletmp[j+1]);

				if ( QPm.type == libmaus2::gamma::FlaggedInterval::interval_type_complete && QPm.from == QPm.to )
				{

					qfiletmp[j] = std::string();
					qfiletmp[j+1] = std::string();
				}
				else
				{
					std::string const outfn = tmpgen.getFileName(true) + ".q";
					libmaus2::gamma::GammaFlaggedIntervalEncoder::unique_ptr_type GPE(new libmaus2::gamma::GammaFlaggedIntervalEncoder(outfn));
					GPE->put(QPm);
					GPE->flush();
					GPE.reset();

					qfiletmp[j] = outfn;
					qfiletmp[j+1] = std::string();

					//std::cerr << "mergeable " << QP0 << "," << QP1 << " -> " << QPm << std::endl;
				}

				j += 2;
				merging = true;
			}
			else
			{
				//std::cerr << "not mergeable " << Atype[j] << "," << Atype[j+1] << std::endl;
				j += 1;
			}
		}

		uint64_t o = 0;
		for ( uint64_t i = 0; i < qfiletmp.size(); ++i )
			if ( qfiletmp[i].size() )
				qfiletmp[o++] = qfiletmp[i];
		qfiletmp.resize(o);
	}

	#if 0
	{
		libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type dec(new libmaus2::gamma::GammaFlaggedIntervalDecoder(qfiletmp,0 /* offset */));
		libmaus2::gamma::FlaggedInterval QP;
		while ( dec->getNext(QP) )
		{
			//std::cerr << "kept (2) " << QP << std::endl;
			assert ( QP.type == libmaus2::gamma::FlaggedInterval::interval_type_complete );
		}
	}
	#endif
}

// compute number of bits set in the two bit vectors Sred and S
uint64_t computeBitsSet(
	std::vector<std::string> const & Sredfn,
	std::vector<std::string> const & Sfn,
	uint64_t const threadpacksize,
	uint64_t const threadpacks,
	uint64_t const n
)
{
	uint64_t volatile s_set = 0;
	libmaus2::parallel::PosixSpinLock s_set_lock;
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		libmaus2::bitio::BitVectorInput Sredin(Sredfn,low);
		libmaus2::bitio::BitVectorInput Sin(Sfn,low);

		uint64_t lc = 0;
		for ( uint64_t i = low; i < high; ++i )
		{
			bool const s0 = Sredin.readBit();
			bool const s1 = Sin.readBit();
			lc += (s0 || s1) ? 1 : 0;
		}

		libmaus2::parallel::ScopePosixSpinLock slock(s_set_lock);
		s_set += lc;
	}

	return s_set;
}

struct AlphabetSymbol
{
	int64_t sym;
	uint64_t id;
	uint64_t subid;
	uint64_t freq;
	bool frequent;

	AlphabetSymbol() {}
	AlphabetSymbol(
		int64_t rsym,
		uint64_t rid,
		uint64_t rfreq,
		bool rfrequent
	) : sym(rsym), id(rid), freq(rfreq), frequent(rfrequent) {}
};

// compute number of files required if freq sum per file is bm
uint64_t computeNF(std::vector<AlphabetSymbol> const & V, uint64_t const bm)
{
	uint64_t nf = 0;
	uint64_t s = 0;

	for ( uint64_t i = 0; i < V.size(); ++i )
	{
		if ( s && s >= bm )
		{
			nf += 1;
			s = 0;
		}

		s += V[i].freq;
	}

	if ( s )
		nf += 1;

	return nf;
}

// compute alphabet split into at most tnumfiles output intervals
static uint64_t computeSplit(uint64_t const tnumfiles, std::vector<AlphabetSymbol> const & V)
{
	uint64_t n = 0;
	for ( uint64_t i = 0; i < V.size(); ++i )
		n += V[i].freq;

	// look for smallest split value s.t. number of files is <= tnumfiles
	uint64_t bl = 0, bh = n;
	while ( bh - bl > 2 )
	{
		uint64_t const bm = (bh + bl)>>1;
		uint64_t const nf = computeNF(V,bm);

		// std::cerr << "bl=" << bl << ",bh=" << bh << ",bm=" << bm << " nf=" << nf << std::endl;

		// number of files too large? bm is not a valid solution
		if ( nf > tnumfiles )
			bl = bm+1;
		// number of files small enough, bm is a valid solution
		else
			bh = bm+1;
	}

	uint64_t split = bl;
	for ( uint64_t i = bl; i < bh; ++i )
		if ( computeNF(V,i) <= tnumfiles )
		{
			split = i;
			break;
		}
	// std::cerr << "bl=" << bl << ",bh=" << bh << " split=" << split << " numfiles=" << computeNF(H,NZH,split) << std::endl;

	return split;
}

/**
 * compute succinct LCP bit vector
 *
 * @param rBWTfn file name of run length encoded BWT
 * @param isafn file name of sampled inverse suffix array
 * @param textfn file name of the text
 * @param tmpgen temporary file name generator
 * @param numthreads number of threads used
 * @param maxrounds maximum number of rounds
 * @param maxmem maximum memory used for sorting and semi external memory LCP computation
 * @param SA (used for debugging purposes if not NULL)
 **/
template<typename input_types_type>
std::vector<std::string> computeSuccinctLCP(
	std::string const & rBWTfn,
	bool const bwtdeleteable,
	std::string const & isafn,
	std::string const & textfn,
	libmaus2::util::TempFileNameGenerator & tmpgen,
	uint64_t const numthreads,
	uint64_t const maxrounds,
	uint64_t const maxmem = 16*1024*1024,
	int32_t * SA = 0
)
{
	libmaus2::timing::RealTimeClock fullrtc; fullrtc.start();

	uint64_t const unsortthreads = numthreads;
	std::vector<std::string> BWTin(1,rBWTfn);
	uint64_t const n = libmaus2::huffman::RLDecoder::getLength(BWTin,numthreads);
	uint64_t const threadpacksize = (n + numthreads - 1)/numthreads;
	uint64_t const threadpacks = (n + threadpacksize - 1)/threadpacksize;
	bool deleteBWTin = false;

	uint64_t const verbthres = 1024;
	bool const verbose = n >= verbthres;

	if ( verbose )
		std::cerr << "[V] processing file is length " << n << std::endl;


	if ( verbose )
		std::cerr << "[V] writing Q to file...";
	std::vector<std::string> qfile(1,tmpgen.getFileName(true) + ".qfile");

	{
		libmaus2::gamma::GammaFlaggedPartitionEncoder GPE(qfile[0]);
		GPE.put(libmaus2::gamma::FlaggedInterval(0,n,libmaus2::gamma::FlaggedInterval::interval_type_complete,true /* active */));
	}
	if ( verbose )
		std::cerr << "done." << std::endl;

	struct RLRadixSortProjector
	{
		static uint64_t project(libmaus2::huffman::RLRun const & P)
	        {
	        	return P.sym;
		}
	};

	// run sort on the BWT to get the key sequences for unsorting
	if ( verbose )
		std::cerr << "[V] computing unsort key sequences...";

	uint64_t maxsym = 0;
	libmaus2::sorting::ParallelRunLengthRadixUnsort::unique_ptr_type unsortinfo(new libmaus2::sorting::ParallelRunLengthRadixUnsort);
	std::vector<std::string> BWTsort = libmaus2::sorting::ParallelRunLengthRadixSort::parallelRadixSort<
		libmaus2::huffman::RLDecoder,
		libmaus2::huffman::RLEncoderStd,
		RLRadixSortProjector,
		libmaus2::huffman::RLEncoderStd,
		libmaus2::huffman::RLDecoder,
		true>
	(
		BWTin,
		numthreads /* num threads */,
		4096 /* max files */,
		false /* delete input */,
		tmpgen,
		4096 /* bs */,
		0,
		libmaus2::huffman::RLDecoder::getLength(BWTin,numthreads),
		0,
		false /* max sym valid */,
		4096 /* keybs */,
		&maxsym,
		unsortinfo.get(),
		unsortthreads /* numthreads */ /* unsort threads */
	);

	if ( verbose )
		std::cerr << "[V] maxsym=" << maxsym << std::endl;

	std::vector < uint64_t > Vsyms;
	libmaus2::parallel::PosixSpinLock Vsymslock;
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		libmaus2::huffman::RLDecoder dec(BWTsort,low,1 /* numthreads */);

		uint64_t todo = high-low;

		libmaus2::huffman::RLDecoder::run_type R;

		std::vector<uint64_t> Lsyms;
		Lsyms.push_back(std::numeric_limits<uint64_t>::max());

		while ( todo )
		{
			dec.decodeRun(R);
			uint64_t const av = std::min(todo,R.rlen);

			if ( static_cast<uint64_t>(R.sym) != Lsyms.back() )
				Lsyms.push_back(R.sym);

			todo -= av;
		}

		std::sort(Lsyms.begin(),Lsyms.end());
		Lsyms.pop_back();

		libmaus2::parallel::ScopePosixSpinLock slock(Vsymslock);
		for ( uint64_t i = 0; i < Lsyms.size(); ++i )
			Vsyms.push_back(Lsyms[i]);
	}

	std::sort(Vsyms.begin(),Vsyms.end());
	Vsyms.resize ( std::unique(Vsyms.begin(),Vsyms.end()) - Vsyms.begin() );

	if ( Vsyms.size() && Vsyms.back()+1 != Vsyms.size() )
	{
		if ( verbose )
		{
			#if 0
			for ( uint64_t i = 0; i < Vsyms.size(); ++i )
				std::cerr << Vsyms[i] << ";";
			std::cerr << std::endl;
			#endif

			std::cerr << "[V] alphabet is not dense, reducing to [0," << Vsyms.size()-1 << "]" << std::endl;
		}

		if ( bwtdeleteable )
			libmaus2::aio::FileRemoval::removeFile(rBWTfn);

		std::vector < std::string > Vbwtrewrite(threadpacks);
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(threadpacks)
		#endif
		for ( uint64_t i = 0; i < threadpacks; ++i )
		{
			uint64_t const low = i*threadpacksize;
			uint64_t const high = std::min(low+threadpacksize,n);

			std::string const outfn = tmpgen.getFileName(true) + ".bwt.rewr";
			Vbwtrewrite[i] = outfn;

			uint64_t symrank = 0;
			int64_t firstsym = -1;

			{
				libmaus2::huffman::RLDecoder dec(BWTsort,low,1 /* numthreads */);
				firstsym = dec.decode();
				assert ( firstsym >= 0 );

				std::vector < uint64_t >::const_iterator ita = ::std::lower_bound(Vsyms.begin(),Vsyms.end(),firstsym);
				assert ( ita != Vsyms.end() );
				assert ( *ita == firstsym );

				symrank = ita - Vsyms.begin();
			}

			uint64_t todo = high-low;

			libmaus2::huffman::RLDecoder::run_type R;

			int64_t prevsym = firstsym;
			libmaus2::huffman::RLDecoder dec(BWTsort,low,1 /* numthreads */);
			libmaus2::huffman::RLEncoderStd enc(outfn,16*1024);
			while ( todo )
			{
				dec.decodeRun(R);
				uint64_t const av = std::min(todo,R.rlen);

				if ( R.sym != prevsym )
				{
					symrank++;
					prevsym = R.sym;
				}

				enc.encodeRun(libmaus2::huffman::RLEncoderStd::run_type(symrank,av));

				todo -= av;
			}

			enc.flush();
		}

		assert ( libmaus2::huffman::RLDecoder::getLength(Vbwtrewrite,numthreads) == n );

		std::vector<std::string> Vun =
			unsortinfo->unsort<
				libmaus2::huffman::RLDecoder,
				libmaus2::huffman::RLEncoderStd,
				libmaus2::huffman::RLDecoder
			>
			(
				Vbwtrewrite,
				true /* delete */,
				tmpgen,
				16*1024 /* enc bs */
			);

		BWTin = Vun;
		deleteBWTin = true;

		unsortinfo.reset();

		maxsym = Vsyms.size()-1;
		libmaus2::sorting::ParallelRunLengthRadixUnsort::unique_ptr_type Tunsortinfo(new libmaus2::sorting::ParallelRunLengthRadixUnsort);
		unsortinfo = UNIQUE_PTR_MOVE(Tunsortinfo);
		std::vector<std::string> BWTsort = libmaus2::sorting::ParallelRunLengthRadixSort::parallelRadixSort<
			libmaus2::huffman::RLDecoder,
			libmaus2::huffman::RLEncoderStd,
			RLRadixSortProjector,
			libmaus2::huffman::RLEncoderStd,
			libmaus2::huffman::RLDecoder,
			true>
		(
			BWTin,
			numthreads /* num threads */,
			4096 /* max files */,
			false /* delete input */,
			tmpgen,
			4096 /* bs */,
			0,
			libmaus2::huffman::RLDecoder::getLength(BWTin,numthreads),
			maxsym,
			true /* max sym valid */,
			4096 /* keybs */,
			0, /* max sym pointer */
			unsortinfo.get(),
			unsortthreads /* numthreads */ /* unsort threads */
		);

		for ( uint64_t i = 0; i < BWTsort.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(BWTsort[i]);

		if ( verbose )
			std::cerr << "[V] minimised alphabet, new maxsym=" << maxsym << std::endl;
	}

	if ( verbose )
		std::cerr << "done." << std::endl;

	for ( uint64_t i = 0; i < BWTsort.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(BWTsort[i]);

	if ( verbose )
		std::cerr << "[V] initialising S bit vector...";

	std::vector<std::string> Sfn(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);
		Sfn[i] = tmpgen.getFileName(true) + ".Sfn";

		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(Sfn[i]));
		for ( uint64_t i = low; i < high; ++i )
			PBout->writeBit(false);
		PBout->flush();
		PBout.reset();
	}

	if ( verbose )
		std::cerr << "done." << std::endl;

	if ( verbose )
		std::cerr << "[V] initialising S reducible bit vector...";

	std::vector<std::string> Sredfn(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);
		Sredfn[i] = tmpgen.getFileName(true) + ".Sfn";

		int64_t prevsym = -1;
		{
			libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,(low+n-1)%n/* offset */,1));
			prevsym = rldec->decode();
		}

		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(Sredfn[i]));
		libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,low/* offset */,1));
		for ( uint64_t i = low; i < high; ++i )
		{
			int64_t cursym = rldec->decode();
			PBout->writeBit(cursym == prevsym);
			prevsym = cursym;
		}
		PBout->flush();
		PBout.reset();
	}

	if ( verbose )
		std::cerr << "done." << std::endl;

	if ( verbose )
		std::cerr << "[V] initialising active bit vector...";

	std::vector<std::string> activefn(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);
		activefn[i] = tmpgen.getFileName(true) + ".active";

		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(activefn[i]));
		for ( uint64_t i = low; i < high; ++i )
			PBout->writeBit(false);
		PBout->flush();
		PBout.reset();
	}

	if ( verbose )
		std::cerr << "done." << std::endl;

	if ( verbose )
		std::cerr << "[V] initialising PD vector...";

	std::vector<std::string> pdfn(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);
		pdfn[i] = tmpgen.getFileName(true) + ".pd";

		libmaus2::gamma::GammaPDEncoder::unique_ptr_type PDout(new libmaus2::gamma::GammaPDEncoder(pdfn[i]));
		for ( uint64_t i = low; i < high; ++i )
			PDout->encode(0);
		PDout->flush();
		PDout.reset();
	}

	if ( verbose )
		std::cerr << "done." << std::endl;

	uint64_t volatile numset = computeBitsSet(Sredfn, Sfn, threadpacksize, threadpacks, n);
	uint64_t volatile unset = n - numset;
	libmaus2::parallel::PosixSpinLock unsetlock;

	uint64_t const logn = std::max(static_cast<uint64_t>(1),static_cast<uint64_t>(::libmaus2::math::ilog(n)));
	uint64_t const d = 2;
	uint64_t const dlogn = d*logn;
	uint64_t const n_div_dlogn = (n + dlogn - 1) / dlogn;

	uint64_t const symbitencoderblocksize = 16*1024;

	if ( verbose )
		std::cerr << "[V] running succinct algorithm until no more than " << n_div_dlogn << " (" << static_cast<double>(n_div_dlogn)/static_cast<double>(n) << ") ranks are unset" << std::endl;

	uint64_t round = 0;
	for ( ; unset && unset > n_div_dlogn && round < maxrounds; round++ )
	{
		if ( verbose )
			std::cerr << "[V] entering round " << round << " fraction of unset ranks " << static_cast<double>(unset)/static_cast<double>(n) << " full " << fullrtc.getElapsedSeconds() << std::endl;

		libmaus2::timing::RealTimeClock rtc;
		rtc.start();

		libmaus2::timing::RealTimeClock roundrtc;
		roundrtc.start();

		libmaus2::timing::RealTimeClock splitrtc; splitrtc.start();
		std::vector < std::string > qflagged(3*threadpacks);
		std::vector < std::string > symbit(threadpacks);

		unsigned int const sc_file_bits = 3;
		uint64_t const sc_output_files = 1ull<<sc_file_bits;
		uint64_t const sc_file_mask = sc_output_files-1;

		unsigned int const sc_sym_bits = libmaus2::math::numbits(maxsym);
		unsigned int const sc_sort_runs = (sc_sym_bits + sc_file_bits - 1) / sc_file_bits;

		libmaus2::autoarray::AutoArray < uint64_t > sc_sym_hist((threadpacks+1) * sc_output_files + 1, false);

		// first round
		{
			libmaus2::timing::RealTimeClock histclock;
			histclock.start();
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < threadpacks; ++t )
			{
				uint64_t * const histp_out = sc_sym_hist.begin() + t * sc_output_files;
				std::fill(histp_out,histp_out+sc_output_files,0ull);

				uint64_t const low = t * threadpacksize;
				uint64_t const high = std::min(low+threadpacksize,n);

				libmaus2::huffman::RLDecoder rldec(BWTin,low/* offset */,1);

				assert ( high > low );
				uint64_t todo = high-low;

				libmaus2::huffman::RLDecoder::run_type SBR;

				while ( todo )
				{
					bool const ok = rldec.decodeRun(SBR);
					assert ( ok );

					uint64_t const av = std::min(todo,SBR.rlen);

					uint64_t const key = SBR.sym;
					uint64_t const nextkey = key & sc_file_mask;
					histp_out[nextkey] += av;
					todo -= av;
				}
			}

			assert ( std::accumulate(sc_sym_hist.begin(),sc_sym_hist.begin() + sc_output_files * threadpacks,0ull) == n );

			std::vector < std::string > symbittmp(sc_output_files * threadpacks);
			std::vector < std::string > qflaggedtmp(3 * sc_output_files * threadpacks);

			std::fill(sc_sym_hist.begin()+sc_output_files * threadpacks,sc_sym_hist.end(),0ull);
			for ( uint64_t i = 0; i < sc_output_files; ++i )
				libmaus2::util::PrefixSums::prefixSums(sc_sym_hist.begin()+i,sc_sym_hist.begin() + sc_output_files * (threadpacks+1) + i,sc_output_files);
			libmaus2::util::PrefixSums::prefixSums(sc_sym_hist.begin()+sc_output_files * threadpacks,sc_sym_hist.end(),1);
			assert ( sc_sym_hist.end()[-1] == n );

			for ( uint64_t i = 0; i < sc_output_files; ++i )
			{
				uint64_t const c = (sc_sym_hist.begin()+sc_output_files * threadpacks)[i];
				for ( uint64_t j = 0; j < threadpacks; ++j )
					sc_sym_hist[ j * sc_output_files + i ] += c;
			}

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < threadpacks; ++t )
			{
				uint64_t * const histp_in = sc_sym_hist.begin() + t * sc_output_files;

				uint64_t const low = t * threadpacksize;
				uint64_t const high = std::min(low+threadpacksize,n);

				std::vector<std::string> lsymbittmp(sc_output_files);
				libmaus2::autoarray::AutoArray < libmaus2::huffman::SymBitEncoderStd::unique_ptr_type > ASBE(sc_output_files);

				std::vector<std::string> lqflaggedtmp(3*sc_output_files);
				libmaus2::autoarray::AutoArray < libmaus2::gamma::GammaFlaggedIntervalEncoder::unique_ptr_type > AGPE(3*sc_output_files);

				for ( uint64_t i = 0; i < sc_output_files; ++i )
				{
					// sym 0 thread 0, sym 0 thread 1, ..., sym 0 thread threadpacks-1, sym 1 thread 0 ...
					lsymbittmp[i] = symbittmp[ i * threadpacks + t ] = tmpgen.getFileName(true) + ".symbit";

					libmaus2::huffman::SymBitEncoderStd::unique_ptr_type SBE(new libmaus2::huffman::SymBitEncoderStd(lsymbittmp[i],symbitencoderblocksize));
					ASBE[i] = UNIQUE_PTR_MOVE(SBE);

					for ( uint64_t j = 0; j < 3; ++j )
					{
						lqflaggedtmp[3*i + j] = qflaggedtmp [ 3 * (i * threadpacks + t) + j ] = tmpgen.getFileName(true) + ".q";
						libmaus2::gamma::GammaFlaggedIntervalEncoder::unique_ptr_type GPE(new libmaus2::gamma::GammaFlaggedIntervalEncoder(lqflaggedtmp[3*i + j]));
						AGPE[3*i + j] = UNIQUE_PTR_MOVE(GPE);
					}
				}

				libmaus2::gamma::GammaFlaggedPartitionDecoder GPD(qfile,low/* offset */,numthreads);
				libmaus2::huffman::RLDecoder rldec(BWTin,low/* offset */,1);
				libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(Sfn,low));

				libmaus2::autoarray::AutoArray<uint64_t> lsymcnt(sc_output_files);

				libmaus2::gamma::FlaggedInterval QP;
				{
					bool const firstok = GPD.getNext(QP);

					assert ( firstok );
					libmaus2::gamma::FlaggedInterval::interval_type const firsttype =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					uint64_t const firstlow = std::max(low,QP.from);
					uint64_t const firsthigh = std::min(high,QP.to);
					bool const firstactive = QP.active;
					assert ( firsthigh > firstlow );

					uint64_t firsttodo = firsthigh-firstlow;

					while ( firsttodo )
					{
						libmaus2::huffman::RLDecoder::run_type R;
						rldec.decodeRun(R);

						uint64_t const use = std::min(R.rlen,firsttodo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = key & sc_file_mask;

						lsymcnt[thiskey] += use;

						for ( uint64_t i = 0; i < use; ++i )
							ASBE[thiskey]->encode(libmaus2::huffman::SymBitRun(R.sym,PBin->readBit()));

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec.putBack(R);
						}

						firsttodo -= use;
					}

					if ( firsttype == libmaus2::gamma::FlaggedInterval::interval_type_complete )
					{
						for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
							if ( lsymcnt[sym] )
							{
								libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],firsttype,firstactive);
								AGPE[3*sym+1]->put(intv);
								histp_in[sym] += lsymcnt[sym];
								lsymcnt[sym] = 0;
							}
					}
					else
					{
						for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
						{
							libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],firsttype,firstactive);
							AGPE[3*sym+0]->put(intv);
							histp_in[sym] += lsymcnt[sym];
							lsymcnt[sym] = 0;
						}
					}
				}

				bool qpok = false;

				while ( (qpok = GPD.getNext(QP)) && QP.to <= high )
				{
					libmaus2::gamma::FlaggedInterval::interval_type type =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					assert ( type == libmaus2::gamma::FlaggedInterval::interval_type_complete );

					uint64_t todo = QP.to-QP.from;
					bool const active = QP.active;

					while ( todo )
					{
						libmaus2::huffman::RLDecoder::run_type R;
						rldec.decodeRun(R);

						uint64_t const use = std::min(R.rlen,todo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = key & sc_file_mask;

						for ( uint64_t i = 0; i < use; ++i )
							ASBE[thiskey]->encode(libmaus2::huffman::SymBitRun(R.sym,PBin->readBit()));

						lsymcnt[thiskey] += use;

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec.putBack(R);
						}

						todo -= use;
					}

					for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
						if ( lsymcnt[sym] )
						{
							libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],type,active);
							AGPE[3*sym+1]->put(intv);
							histp_in[sym] += lsymcnt[sym];
							lsymcnt[sym] = 0;
						}
				}

				if ( qpok && QP.from < high )
				{
					assert ( QP.to > high );

					libmaus2::gamma::FlaggedInterval::interval_type lasttype =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					assert ( lasttype != libmaus2::gamma::FlaggedInterval::interval_type_complete );

					uint64_t const lastlow = std::max(low,QP.from);
					uint64_t const lasthigh = std::min(high,QP.to);
					assert ( lasthigh > lastlow );
					bool const lastactive = QP.active;

					uint64_t lasttodo = lasthigh-lastlow;

					while ( lasttodo )
					{
						libmaus2::huffman::RLDecoder::run_type R;
						rldec.decodeRun(R);

						uint64_t const use = std::min(R.rlen,lasttodo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = key & sc_file_mask;

						for ( uint64_t i = 0; i < use; ++i )
							ASBE[thiskey]->encode(libmaus2::huffman::SymBitRun(R.sym,PBin->readBit()));

						lsymcnt[thiskey] += use;

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec.putBack(R);
						}

						lasttodo -= use;
					}

					for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
					{
						libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],lasttype,lastactive);
						AGPE[3*sym+2]->put(intv);
						histp_in[sym] += lsymcnt[sym];
						lsymcnt[sym] = 0;
					}
				}

				for ( uint64_t i = 0; i < ASBE.size(); ++i )
				{
					ASBE[i]->flush();
					ASBE[i].reset();
				}
				for ( uint64_t i = 0; i < AGPE.size(); ++i )
				{
					AGPE[i]->flush();
					AGPE[i].reset();
				}
			}

			mergeQueueFiles(qflaggedtmp,numthreads,tmpgen);

			for ( uint64_t i = 0; i < symbit.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(symbit[i]);
			symbit = symbittmp;

			for ( uint64_t i = 0; i < qflagged.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(qflagged[i]);
			qflagged = qflaggedtmp;
		}

		if ( verbose )
			std::cerr << "\t[V] rewrote to SymBit in time " << splitrtc.getElapsedSeconds() << std::endl;

		libmaus2::timing::RealTimeClock qmergertc; qmergertc.start();
		mergeQueueFiles(qflagged,numthreads,tmpgen);

		if ( verbose )
			std::cerr << "\t[V] merged queue files in time " << qmergertc.getElapsedSeconds() << std::endl;

		// sort rounds over sc
		libmaus2::timing::RealTimeClock sortscrtc; sortscrtc.start();
		for ( uint64_t sc_round = 1; sc_round < sc_sort_runs; ++sc_round )
		{
			libmaus2::timing::RealTimeClock histclock;
			histclock.start();
			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < threadpacks; ++t )
			{
				uint64_t * const histp_out = sc_sym_hist.begin() + t * sc_output_files;
				std::fill(histp_out,histp_out+sc_output_files,0ull);

				uint64_t const low = t * threadpacksize;
				uint64_t const high = std::min(low+threadpacksize,n);

				libmaus2::huffman::SymBitDecoder::unique_ptr_type rldec(new libmaus2::huffman::SymBitDecoder(symbit,low/* offset */,1));

				assert ( high > low );
				uint64_t todo = high-low;

				libmaus2::huffman::SymBitRun SBR;

				while ( todo )
				{
					bool const ok = rldec->decodeRun(SBR);
					assert ( ok );

					uint64_t const av = std::min(todo,SBR.rlen);

					uint64_t const key = SBR.sym;
					uint64_t const nextkey = (key >> ((sc_round)*sc_file_bits)) & sc_file_mask;
					histp_out[nextkey] += av;
					todo -= av;
				}
			}

			assert ( std::accumulate(sc_sym_hist.begin(),sc_sym_hist.begin() + sc_output_files * threadpacks,0ull) == n );

			std::vector < std::string > symbittmp(sc_output_files * threadpacks);
			std::vector < std::string > qflaggedtmp(3 * sc_output_files * threadpacks);

			std::fill(sc_sym_hist.begin()+sc_output_files * threadpacks,sc_sym_hist.end(),0ull);
			for ( uint64_t i = 0; i < sc_output_files; ++i )
				libmaus2::util::PrefixSums::prefixSums(sc_sym_hist.begin()+i,sc_sym_hist.begin() + sc_output_files * (threadpacks+1) + i,sc_output_files);
			libmaus2::util::PrefixSums::prefixSums(sc_sym_hist.begin()+sc_output_files * threadpacks,sc_sym_hist.end(),1);
			assert ( sc_sym_hist.end()[-1] == n );

			for ( uint64_t i = 0; i < sc_output_files; ++i )
			{
				uint64_t const c = (sc_sym_hist.begin()+sc_output_files * threadpacks)[i];
				for ( uint64_t j = 0; j < threadpacks; ++j )
					sc_sym_hist[ j * sc_output_files + i ] += c;
			}

			#if defined(_OPENMP)
			#pragma omp parallel for num_threads(numthreads)
			#endif
			for ( uint64_t t = 0; t < threadpacks; ++t )
			{
				uint64_t * const histp_in = sc_sym_hist.begin() + t * sc_output_files;

				uint64_t const low = t * threadpacksize;
				uint64_t const high = std::min(low+threadpacksize,n);

				std::vector<std::string> lsymbittmp(sc_output_files);
				libmaus2::autoarray::AutoArray < libmaus2::huffman::SymBitEncoderStd::unique_ptr_type > ASBE(sc_output_files);

				std::vector<std::string> lqflaggedtmp(3*sc_output_files);
				libmaus2::autoarray::AutoArray < libmaus2::gamma::GammaFlaggedIntervalEncoder::unique_ptr_type > AGPE(3*sc_output_files);

				for ( uint64_t i = 0; i < sc_output_files; ++i )
				{
					// sym 0 thread 0, sym 0 thread 1, ..., sym 0 thread threadpacks-1, sym 1 thread 0 ...
					lsymbittmp[i] = symbittmp[ i * threadpacks + t ] = tmpgen.getFileName(true) + ".symbit";

					libmaus2::huffman::SymBitEncoderStd::unique_ptr_type SBE(new libmaus2::huffman::SymBitEncoderStd(lsymbittmp[i],symbitencoderblocksize));
					ASBE[i] = UNIQUE_PTR_MOVE(SBE);

					for ( uint64_t j = 0; j < 3; ++j )
					{
						lqflaggedtmp[3*i + j] = qflaggedtmp [ 3 * (i * threadpacks + t) + j ] = tmpgen.getFileName(true) + ".q";
						libmaus2::gamma::GammaFlaggedIntervalEncoder::unique_ptr_type GPE(new libmaus2::gamma::GammaFlaggedIntervalEncoder(lqflaggedtmp[3*i + j]));
						AGPE[3*i + j] = UNIQUE_PTR_MOVE(GPE);
					}
				}

				libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type GPD(new libmaus2::gamma::GammaFlaggedIntervalDecoder(qflagged,low/* offset */));
				libmaus2::huffman::SymBitDecoder::unique_ptr_type rldec(new libmaus2::huffman::SymBitDecoder(symbit,low/* offset */,1));

				libmaus2::autoarray::AutoArray<uint64_t> lsymcnt(sc_output_files);

				libmaus2::gamma::FlaggedInterval QP;
				{
					bool const firstok = GPD->getNext(QP);

					assert ( firstok );
					libmaus2::gamma::FlaggedInterval::interval_type const firsttype =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					uint64_t const firstlow = std::max(low,QP.from);
					uint64_t const firsthigh = std::min(high,QP.to);
					bool const firstactive = QP.active;
					assert ( firsthigh > firstlow );

					uint64_t firsttodo = firsthigh-firstlow;

					while ( firsttodo )
					{
						libmaus2::huffman::SymBitRun R;
						rldec->decodeRun(R);

						uint64_t const use = std::min(R.rlen,firsttodo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = (key >> ((sc_round+0)*sc_file_bits)) & sc_file_mask;

						lsymcnt[thiskey] += use;

						ASBE[thiskey]->encodeRun(libmaus2::huffman::SymBitRun(R.sym,R.sbit,use));

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec->putBack(R);
						}

						firsttodo -= use;
					}

					if ( firsttype == libmaus2::gamma::FlaggedInterval::interval_type_complete )
					{
						for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
							if ( lsymcnt[sym] )
							{
								libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],firsttype,firstactive);
								AGPE[3*sym+1]->put(intv);
								histp_in[sym] += lsymcnt[sym];
								lsymcnt[sym] = 0;
							}
					}
					else
					{
						for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
						{
							libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],firsttype,firstactive);
							AGPE[3*sym+0]->put(intv);
							histp_in[sym] += lsymcnt[sym];
							lsymcnt[sym] = 0;
						}
					}
				}

				bool qpok = false;

				while ( (qpok = GPD->getNext(QP)) && QP.to <= high )
				{
					libmaus2::gamma::FlaggedInterval::interval_type type =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					assert ( type == libmaus2::gamma::FlaggedInterval::interval_type_complete );

					uint64_t todo = QP.to-QP.from;
					bool const active = QP.active;

					while ( todo )
					{
						libmaus2::huffman::SymBitRun R;
						rldec->decodeRun(R);

						uint64_t const use = std::min(R.rlen,todo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = (key >> ((sc_round+0)*sc_file_bits)) & sc_file_mask;

						ASBE[thiskey]->encodeRun(libmaus2::huffman::SymBitRun(R.sym,R.sbit,use));

						lsymcnt[thiskey] += use;

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec->putBack(R);
						}

						todo -= use;
					}

					for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
						if ( lsymcnt[sym] )
						{
							libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],type,active);
							AGPE[3*sym+1]->put(intv);
							histp_in[sym] += lsymcnt[sym];
							lsymcnt[sym] = 0;
						}
				}

				if ( qpok && QP.from < high )
				{
					assert ( QP.to > high );

					libmaus2::gamma::FlaggedInterval::interval_type lasttype =
						libmaus2::gamma::FlaggedInterval::getType(low,high,QP.from,QP.to);

					assert ( lasttype != libmaus2::gamma::FlaggedInterval::interval_type_complete );

					uint64_t const lastlow = std::max(low,QP.from);
					uint64_t const lasthigh = std::min(high,QP.to);
					assert ( lasthigh > lastlow );
					bool const lastactive = QP.active;

					uint64_t lasttodo = lasthigh-lastlow;

					while ( lasttodo )
					{
						libmaus2::huffman::SymBitRun R;
						rldec->decodeRun(R);

						uint64_t const use = std::min(R.rlen,lasttodo);

						uint64_t const key = R.sym;
						uint64_t const thiskey = (key >> ((sc_round+0)*sc_file_bits)) & sc_file_mask;

						ASBE[thiskey]->encodeRun(libmaus2::huffman::SymBitRun(R.sym,R.sbit,use));

						lsymcnt[thiskey] += use;

						if ( use != R.rlen )
						{
							R.rlen -= use;
							rldec->putBack(R);
						}

						lasttodo -= use;
					}

					for ( uint64_t sym = 0; sym < sc_output_files; ++sym )
					{
						libmaus2::gamma::FlaggedInterval intv(histp_in[sym],histp_in[sym]+lsymcnt[sym],lasttype,lastactive);
						AGPE[3*sym+2]->put(intv);
						histp_in[sym] += lsymcnt[sym];
						lsymcnt[sym] = 0;
					}
				}

				for ( uint64_t i = 0; i < ASBE.size(); ++i )
				{
					ASBE[i]->flush();
					ASBE[i].reset();
				}
				for ( uint64_t i = 0; i < AGPE.size(); ++i )
				{
					AGPE[i]->flush();
					AGPE[i].reset();
				}
			}

			mergeQueueFiles(qflaggedtmp,numthreads,tmpgen);

			for ( uint64_t i = 0; i < symbit.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(symbit[i]);
			symbit = symbittmp;

			for ( uint64_t i = 0; i < qflagged.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(qflagged[i]);
			qflagged = qflaggedtmp;
		}

		if ( verbose && sc_sort_runs > 1 )
			std::cerr << "\t[V] sorted SymBit in time " << sortscrtc.getElapsedSeconds() << std::endl;

		// compute T vector
		// extract and combine S bits from sorted and unsorted whenever there is a count set
		libmaus2::timing::RealTimeClock trtc; trtc.start();
		std::vector<std::string> Tfn(threadpacks);
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t t = 0; t < threadpacks; ++t )
		{
			uint64_t const low = t * threadpacksize;
			uint64_t const high = std::min(low+threadpacksize,n);

			libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type GPD(new libmaus2::gamma::GammaFlaggedIntervalDecoder(qflagged,low/* offset */));
			libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(Sfn,low));
			libmaus2::huffman::SymBitDecoder::unique_ptr_type rldec(new libmaus2::huffman::SymBitDecoder(symbit,low/* offset */,1));
			libmaus2::gamma::FlaggedInterval QP;

			std::string const outfn = tmpgen.getFileName(true) + ".T";
			Tfn[t] = outfn;
			libmaus2::huffman::RLEncoderStd::unique_ptr_type BVO(new libmaus2::huffman::RLEncoderStd(outfn,4096 /* bufsize */));

			{
				bool const firstok = GPD->getNext(QP);
				assert ( firstok );

				uint64_t const qlow = std::max(low,QP.from);
				uint64_t const qhigh = std::min(high,QP.to);
				assert ( qhigh > qlow );
				uint64_t qtodo = qhigh-qlow;

				if ( QP.from == low )
				{
					// read S bit
					bool const Stgtbit = PBin->readBit();
					libmaus2::huffman::SymBit SB;
					bool const ok = rldec->decode(SB);
					assert ( ok );

					bool const Tbit = (!(SB.sbit)) && (!Stgtbit);

					// write T bit
					BVO->encode(Tbit);

					qtodo -= 1;
				}

				while ( qtodo )
				{
					libmaus2::huffman::SymBitRun SBR;
					bool const ok = rldec->decodeRun(SBR);
					assert ( ok );

					uint64_t const av = std::min(SBR.rlen,qtodo);
					SBR.rlen -= av;

					for ( uint64_t i = 0; i < av; ++i )
						PBin->readBit();

					BVO->encodeRun(libmaus2::huffman::RLEncoderStd::run_type(false,av));

					if ( SBR.rlen )
						rldec->putBack(SBR);

					qtodo -= av;
				}
			}

			bool qpok;
			while ( (qpok=GPD->getNext(QP)) && QP.to <= high )
			{
				assert ( QP.from >= low );
				assert ( QP.to <= high );

				uint64_t qtodo = QP.to-QP.from;

				{
					// read S bit
					bool const Stgtbit = PBin->readBit();
					libmaus2::huffman::SymBit SB;
					bool const ok = rldec->decode(SB);
					assert ( ok );

					bool const Tbit = (!(SB.sbit)) && (!Stgtbit);

					// write T bit
					BVO->encode(Tbit);

					qtodo -= 1;
				}

				while ( qtodo )
				{
					libmaus2::huffman::SymBitRun SBR;
					bool const ok = rldec->decodeRun(SBR);
					assert ( ok );

					uint64_t const av = std::min(SBR.rlen,qtodo);
					SBR.rlen -= av;

					for ( uint64_t i = 0; i < av; ++i )
						PBin->readBit();
					BVO->encodeRun(libmaus2::huffman::RLEncoderStd::run_type(false,av));

					if ( SBR.rlen )
						rldec->putBack(SBR);

					qtodo -= av;
				}
			}

			if ( qpok && QP.from < high )
			{
				assert ( QP.to > high );

				uint64_t qtodo = high-QP.from;

				{
					// read S bit
					bool const Stgtbit = PBin->readBit();
					libmaus2::huffman::SymBit SB;
					bool const ok = rldec->decode(SB);
					assert ( ok );

					bool const Tbit = (!(SB.sbit)) && (!Stgtbit);

					// write T bit
					BVO->encode(Tbit);

					qtodo -= 1;
				}

				while ( qtodo )
				{
					libmaus2::huffman::SymBitRun SBR;
					bool const ok = rldec->decodeRun(SBR);
					assert ( ok );

					uint64_t const av = std::min(SBR.rlen,qtodo);
					SBR.rlen -= av;

					for ( uint64_t i = 0; i < av; ++i )
						PBin->readBit();
					BVO->encodeRun(libmaus2::huffman::RLEncoderStd::run_type(false,av));

					if ( SBR.rlen )
						rldec->putBack(SBR);

					qtodo -= av;
				}
			}

			BVO->flush();
			BVO.reset();
		}

		assert ( libmaus2::huffman::RLDecoder::getLength(Tfn,1) == n );

		if ( verbose )
			std::cerr << "\t[V] computed T vector in time " << trtc.getElapsedSeconds() << std::endl;

		for ( uint64_t i = 0; i < symbit.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(symbit[i]);

		libmaus2::timing::RealTimeClock unsortrtc; unsortrtc.start();
		// unradixsort resulting bit vector
		std::vector<std::string> const Toutfn = unsortinfo->unsort<libmaus2::huffman::RLDecoder,libmaus2::huffman::RLEncoderStd,libmaus2::huffman::RLDecoder>(
			Tfn,true /* delete input */,tmpgen,4096 /* enc bs */);

		for ( uint64_t i = 0; i < Tfn.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Tfn[i]);

		if ( verbose )
			std::cerr << "\t[V] T bit vector unsorted in time " << unsortrtc.getElapsedSeconds() << std::endl;

		// update
		libmaus2::timing::RealTimeClock outintvrtc; outintvrtc.start();
		std::vector<std::string> pdfntmp(threadpacks);
		std::vector<std::string> soutfntmp(threadpacks);
		std::vector<std::string> qoutfntmp(threadpacks);
		std::vector<std::string> activefntmp(threadpacks);
		uint64_t volatile g_outintv = 0;
		libmaus2::parallel::PosixSpinLock g_outintv_lock;
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t t = 0; t < threadpacks; ++t )
		{
			uint64_t const low = t * threadpacksize;
			uint64_t const high = std::min(low+threadpacksize,n);

			std::string const activeoutfn = tmpgen.getFileName(true) + ".active";
			activefntmp[t] = activeoutfn;

			std::string const soutfn = tmpgen.getFileName(true) + ".S";
			soutfntmp[t] = soutfn;
			std::string const qoutfn = tmpgen.getFileName(true) + ".q";
			qoutfntmp[t] = qoutfn;

			libmaus2::gamma::GammaPDDecoder::unique_ptr_type PDin(new libmaus2::gamma::GammaPDDecoder(pdfn,low));
			std::string const pdoutfn = tmpgen.getFileName(true) + ".pdfn";
			libmaus2::gamma::GammaPDEncoder::unique_ptr_type PDout(new libmaus2::gamma::GammaPDEncoder(pdoutfn));
			pdfntmp[t] = pdoutfn;

			libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(soutfn));
			libmaus2::bitio::BitVectorOutput::unique_ptr_type PactiveOut(new libmaus2::bitio::BitVectorOutput(activeoutfn));
			libmaus2::gamma::GammaFlaggedPartitionEncoder::unique_ptr_type GPE(new libmaus2::gamma::GammaFlaggedPartitionEncoder(qoutfn));

			libmaus2::huffman::RLDecoder::unique_ptr_type TBVI(new libmaus2::huffman::RLDecoder(Toutfn,low,1 /* numthreads */));

			// active vector input
			libmaus2::bitio::BitVectorInput::unique_ptr_type PactiveIn(new libmaus2::bitio::BitVectorInput(activefn,low));
			// S vector
			libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(Sfn,low));
			// pre Q
			libmaus2::gamma::GammaFlaggedIntervalDecoder::unique_ptr_type GPD(new libmaus2::gamma::GammaFlaggedIntervalDecoder(qflagged,low/* offset */));

			libmaus2::gamma::FlaggedInterval prevP;
			bool prevPvalid = false;
			uint64_t outintv = 0;

			libmaus2::gamma::FlaggedInterval QP;
			{

				bool const firstok = GPD->getNext(QP);
				assert ( firstok );
				uint64_t const firstlow = std::max(low,QP.from);
				uint64_t const firsthigh = std::min(high,QP.to);
				assert ( firsthigh > firstlow );
				uint64_t qtodo = firsthigh - firstlow;

				if ( QP.from == low )
				{
					// S
					bool const Sinbit = PBin->readBit();
					// active
					bool const Tinbit = TBVI->decode();
					bool const activeInBit = PactiveIn->readBit() || Tinbit;

					uint64_t const p = PDin->decode();
					PDout->encode(p + (activeInBit?1:0));

					// next interval
					libmaus2::gamma::FlaggedInterval nextP(QP.from,QP.to,libmaus2::gamma::FlaggedInterval::interval_type_complete,(!Sinbit) /* active XXX check */);

					prevP = nextP;
					prevPvalid = true;

					if ( (!Sinbit) )
					{
						PactiveOut->writeBit(false); // deactivate
						PBout->writeBit(true); // S[i] = true
					}
					else
					{
						PactiveOut->writeBit(activeInBit);
						PBout->writeBit(Sinbit);
					}

					qtodo -= 1;
				}

				while ( qtodo-- )
				{
					bool const Tinbit = TBVI->decode();
					bool const activeInBit = PactiveIn->readBit() || Tinbit;

					uint64_t const p = PDin->decode();
					PDout->encode(p + (activeInBit?1:0));

					PactiveOut->writeBit(activeInBit);
					PBout->writeBit(PBin->readBit());
				}
			}

			bool qok;
			while ( (qok=GPD->getNext(QP)) && QP.to <= high )
			{
				assert ( QP.from >= low );

				uint64_t qtodo = QP.to-QP.from;

				// S
				bool const Sinbit = PBin->readBit();
				// active
				bool const Tinbit = TBVI->decode();
				bool const activeInBit = PactiveIn->readBit() || Tinbit;

				uint64_t const p = PDin->decode();
				PDout->encode(p + (activeInBit?1:0));

				libmaus2::gamma::FlaggedInterval nextP(QP.from,QP.to,libmaus2::gamma::FlaggedInterval::interval_type_complete,(!Sinbit) /* active XXX check */);

				if ( prevPvalid )
				{
					if ( prevP.active || (nextP.active != prevP.active) )
					{
						GPE->put(prevP);
						outintv += 1;
						prevP = nextP;
					}
					else
					{
						assert ( ! prevP.active );
						assert ( ! nextP.active );

						prevP.to = nextP.to;
					}
				}
				else
				{
					prevP = nextP;
					prevPvalid = true;
				}

				if ( (!Sinbit) )
				{
					PactiveOut->writeBit(false); // deactivate
					PBout->writeBit(true); // S[i] = true
				}
				else
				{
					PactiveOut->writeBit(activeInBit);
					PBout->writeBit(Sinbit);
				}

				qtodo -= 1;

				while ( qtodo-- )
				{
					bool const Tinbit = TBVI->decode();
					bool const activeInBit = PactiveIn->readBit() || Tinbit;

					uint64_t const p = PDin->decode();
					PDout->encode(p + (activeInBit?1:0));

					PactiveOut->writeBit(activeInBit);
					PBout->writeBit(PBin->readBit());
				}
			}

			if ( qok && QP.from < high )
			{
				assert ( QP.from >= low );

				uint64_t qtodo = high - QP.from;

				// S
				bool const Sinbit = PBin->readBit();
				// active
				bool const Tinbit = TBVI->decode();
				bool const activeInBit = PactiveIn->readBit() || Tinbit;

				uint64_t const p = PDin->decode();
				PDout->encode(p + (activeInBit?1:0));

				libmaus2::gamma::FlaggedInterval nextP(QP.from,QP.to,libmaus2::gamma::FlaggedInterval::interval_type_complete,(!Sinbit) /* active XXX check */);

				if ( prevPvalid )
				{
					if ( prevP.active || (nextP.active != prevP.active) )
					{
						GPE->put(prevP);
						outintv += 1;
						prevP = nextP;
					}
					else
					{
						assert ( ! prevP.active );
						assert ( ! nextP.active );

						prevP.to = nextP.to;
					}
				}
				else
				{
					prevP = nextP;
					prevPvalid = true;
				}

				if ( (!Sinbit) )
				{
					PactiveOut->writeBit(false); // deactivate
					PBout->writeBit(true); // S[i] = true
				}
				else
				{
					PactiveOut->writeBit(activeInBit);
					PBout->writeBit(Sinbit);
				}

				qtodo -= 1;

				while ( qtodo-- )
				{
					bool const Tinbit = TBVI->decode();
					bool const activeInBit = PactiveIn->readBit() || Tinbit;

					uint64_t const p = PDin->decode();
					PDout->encode(p + (activeInBit?1:0));

					PactiveOut->writeBit(activeInBit);
					PBout->writeBit(PBin->readBit());
				}
			}

			if ( prevPvalid && prevP.from >= low && prevP.from < high )
			{
				//std::cerr << "put last " << prevP << std::endl;
				GPE->put(prevP);
				outintv += 1;
			}

			g_outintv_lock.lock();
			g_outintv += outintv;
			g_outintv_lock.unlock();

			PBout->flush();
			PactiveOut->flush();
			GPE->flush();

			PDout->flush();
			PDout.reset();
		}

		// remove old files
		for ( uint64_t i = 0; i < Toutfn.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Toutfn[i]);

		for ( uint64_t i = 0; i < pdfn.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(pdfn[i]);
		pdfn = pdfntmp;

		for ( uint64_t i = 0; i < qflagged.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(qflagged[i]);
		if ( verbose )
			std::cerr << "\t[V] cleaned active/updated S/updated PD/produced queue for next round in time " << outintvrtc.getElapsedSeconds() << std::endl;

		if ( verbose )
			std::cerr << "\t[V] output intervals " << g_outintv << std::endl;

		for ( uint64_t i = 0; i < Sfn.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Sfn[i]);
		Sfn = soutfntmp;

		for ( uint64_t i = 0; i < activefn.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(activefn[i]);
		activefn = activefntmp;

		for ( uint64_t i = 0; i < qfile.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(qfile[i]);
		qfile = qoutfntmp;

		numset = computeBitsSet(Sredfn, Sfn, threadpacksize, threadpacks, n);
		unset = n - numset;

		if ( verbose )
			std::cerr << "\t[V] round completed in time " << roundrtc.getElapsedSeconds() << std::endl;
	}

	// remove unsorting key files
	unsortinfo.reset();

	// remove queue file(s)
	for ( uint64_t i = 0; i < qfile.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(qfile[i]);

	for ( uint64_t i = 0; i < Sredfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Sredfn[i]);

	// filter values from PD which are for ranks still in the active set
	if ( verbose )
		std::cerr << "[V] filtering PD values still in active set" << std::endl;
	std::vector < std::string > pdfntmp(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < threadpacks; ++t )
	{
		uint64_t const low = t * threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		std::string const outfn = tmpgen.getFileName(true) + ".pd";
		pdfntmp[t] = outfn;

		libmaus2::gamma::GammaPDDecoder::unique_ptr_type PDin(new libmaus2::gamma::GammaPDDecoder(pdfn,low));
		libmaus2::bitio::BitVectorInput::unique_ptr_type PactiveIn(new libmaus2::bitio::BitVectorInput(activefn,low));
		libmaus2::gamma::GammaPDEncoder::unique_ptr_type PDout(new libmaus2::gamma::GammaPDEncoder(outfn));

		for ( uint64_t i = low; i < high; ++i )
		{
			uint64_t v;
			bool const pdok = PDin->decode(v);
			assert ( pdok );
			bool const active = PactiveIn->readBit();

			if ( active )
			{
				PDout->encode(0);
			}
			else
			{
				PDout->encode(v);
			}
		}

		PDout->flush();
		PDout.reset();
	}
	for ( uint64_t i = 0; i < activefn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(activefn[i]);
	for ( uint64_t i = 0; i < pdfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(pdfn[i]);
	pdfn = pdfntmp;

	// mark all reducible ranks as set
	if ( verbose )
		std::cerr << "[V] marking reducible ranks as set" << std::endl;
	std::vector < std::string > sirredtmp(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < threadpacks; ++t )
	{
		uint64_t const low = t * threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		std::string const outfn = tmpgen.getFileName(true) + ".sirred";
		sirredtmp[t] = outfn;

		int64_t prevsym;

		{
			libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,(low+n-1)%n/* offset */,1/*numthreads*/));
			prevsym = rldec->decode();
			assert ( prevsym >= 0 );
		}

		// S bit input
		libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,low/* offset */,1/*numthreads*/));
		libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(Sfn,low));
		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(outfn));

		for ( uint64_t i = low; i < high; ++i )
		{
			int64_t const sym = rldec->decode();
			bool const inbit = PBin->readBit();
			bool const sbit = inbit || (sym == prevsym);
			PBout->writeBit(sbit);

			#if 0
			if ( sbit )
			{
				uint64_t const pos = SA[i];
				uint64_t const prevpos = (pos+n-1)%n;
				uint64_t const prevrank = ISA[prevpos];
				//std::cerr << "rank " << i << " marked as set" << " sym=" << sym << " prevsym=" << prevsym << " inbit=" << inbit << " LCP[]=" << LCP[i] << " LCP[prev]=" << LCP[prevrank] << std::endl;
			}
			else
			{
				//std::cerr << "rank " << i << " missing" << std::endl;
			}
			#endif

			prevsym = sym;
		}

		PBout->flush();
		PBout.reset();
	}
	for ( uint64_t i = 0; i < Sfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Sfn[i]);
	Sfn = sirredtmp;

	// extract missing ranks
	if ( verbose )
		std::cerr << "[V] extracting bit vector of missing ranks" << std::endl;
	std::vector < std::string > smisstmp(threadpacks);
	std::vector < std::string > slftmp(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < threadpacks; ++t )
	{
		uint64_t const low = t * threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		std::string const outfn = tmpgen.getFileName(true) + ".slftmp";
		slftmp[t] = outfn;

		std::string const misoutfn = tmpgen.getFileName(true) + ".smisstmp";
		smisstmp[t] = misoutfn;

		// S bit input
		libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,low/* offset */,1/*numthreads*/));
		libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(Sfn,low));
		libmaus2::huffman::LFSetBitEncoder::unique_ptr_type lfenc(new libmaus2::huffman::LFSetBitEncoder(outfn,4096));
		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(misoutfn));

		for ( uint64_t i = low; i < high; ++i )
		{
			bool const sbit = PBin->readBit();
			int64_t const sym = rldec->decode();
			bool const missing = !sbit;
			lfenc->encode(libmaus2::huffman::LFSetBit(sym,missing));
			PBout->writeBit(missing);
			#if 0
			if ( missing )
				std::cerr << "RANK MIS " << i << std::endl;
			#endif
		}

		lfenc->flush();
		lfenc.reset();
		PBout->flush();
		PBout.reset();
	}
	for ( uint64_t i = 0; i < Sfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Sfn[i]);

	// perform lf on missing ranks to get ranks for previous positions
	if ( verbose )
		std::cerr << "[V] performing LF on missing ranks" << std::endl;
	{
		struct LFSetBitProjector
		{
			static uint64_t project(libmaus2::huffman::LFSetBit const & S)
			{
				return S.sym;
			}
		};

		slftmp = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFSetBitDecoder,
			libmaus2::huffman::LFSetBitEncoder,
			LFSetBitProjector
		>(
			slftmp,
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			n, /* ihigh */
			maxsym, /* max sym */
			true /* max sym valid */
		);
	}

	// merge
	if ( verbose )
		std::cerr << "[V] merging missing ranks with their LF and adding previous ranks" << std::endl;

	std::vector < std::string > smissmergetmp(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < threadpacks; ++t )
	{
		uint64_t const low = t * threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		std::string const misoutfn = tmpgen.getFileName(true) + ".smissmergetmp";
		smissmergetmp[t] = misoutfn;

		libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(smisstmp,low));
		libmaus2::huffman::LFSetBitDecoder lfdec(slftmp,low);
		libmaus2::bitio::BitVectorOutput::unique_ptr_type PBout(new libmaus2::bitio::BitVectorOutput(misoutfn));

		bool curbit = false;
		if ( low < high )
		{
			bool const sbit = PBin->readBit();
			libmaus2::huffman::LFSetBit lfbit;
			bool const lfok = lfdec.decode(lfbit);
			assert ( lfok );

			curbit = sbit || lfbit.sbit;
		}

		for ( uint64_t i = low+1; i < high; ++i )
		{
			bool const sbit = PBin->readBit();
			libmaus2::huffman::LFSetBit lfbit;
			bool const lfok = lfdec.decode(lfbit);
			assert ( lfok );

			bool const nextbit = sbit || lfbit.sbit;

			PBout->writeBit(curbit || nextbit);

			curbit = nextbit;
		}

		if ( high != n )
		{
			bool const sbit = PBin->readBit();
			libmaus2::huffman::LFSetBit lfbit;
			bool const lfok = lfdec.decode(lfbit);
			assert ( lfok );

			bool const nextbit = sbit || lfbit.sbit;

			PBout->writeBit(curbit || nextbit);
		}
		else
		{
			libmaus2::bitio::BitVectorInput::unique_ptr_type PBin(new libmaus2::bitio::BitVectorInput(smisstmp,0));
			libmaus2::huffman::LFSetBitDecoder lfdec(slftmp,0);

			bool const sbit = PBin->readBit();
			libmaus2::huffman::LFSetBit lfbit;
			bool const lfok = lfdec.decode(lfbit);
			assert ( lfok );

			bool const nextbit = sbit || lfbit.sbit;

			PBout->writeBit(curbit || nextbit);
		}

		PBout->flush();
		PBout.reset();
	}

	for ( uint64_t i = 0; i < smisstmp.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(smisstmp[i]);
	// remove LF files
	for ( uint64_t i = 0; i < slftmp.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(slftmp[i]);

	if ( verbose )
		std::cerr << "[V] computing character histograms for thread blocks" << std::endl;

	std::map<int64_t,uint64_t> histmap;
	libmaus2::parallel::PosixMutex histmaplock;
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < threadpacks; ++t )
	{
		uint64_t const low = t * threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		libmaus2::huffman::RLDecoder::unique_ptr_type rldec(new libmaus2::huffman::RLDecoder(BWTin,low/* offset */,1/*numthreads*/));
		libmaus2::util::Histogram lhist;

		uint64_t todo = high-low;
		libmaus2::huffman::RLDecoder::run_type R;
		while ( todo )
		{
			rldec->decodeRun(R);
			uint64_t const av = std::min(R.rlen,todo);
			lhist.add(R.sym,av);
			todo -= av;
		}

		libmaus2::parallel::ScopePosixMutex slock(histmaplock);
		std::map<int64_t,uint64_t> lmap = lhist.getByType<int64_t>();
		for ( std::map<int64_t,uint64_t>::const_iterator ita = lmap.begin(); ita != lmap.end(); ++ita )
			histmap[ita->first] += ita->second;
	}
	std::string const histfn = tmpgen.getFileName(true) + ".hist";
	{
		libmaus2::aio::OutputStreamInstance OSI(histfn);
		libmaus2::util::NumberMapSerialisation::serialiseMap(OSI,histmap);
		OSI.flush();
	}

	if ( verbose )
		std::cerr << "[V] extracting/selecting SA values required for missing ranks" << std::endl;

	libmaus2::timing::RealTimeClock selectrtc; selectrtc.start();
	std::string const sasub = tmpgen.getFileName(true) + ".sasub";
	std::ostream * nullstr = 0;
	libmaus2::suffixsort::bwtb3m::BwtSelectSSA::computeSSA(BWTin,
		histfn,
		isafn,
		sasub,
		smissmergetmp,
		tmpgen.getFileName(true) + ".selectssa",
		false,
		numthreads,
		maxmem,
		4096 /* max temp files */,
		verbose ? (&(std::cerr)) : nullstr
	);

	for ( uint64_t i = 0; i < smissmergetmp.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(smissmergetmp[i]);

	if ( verbose )
		std::cerr << "[V] selected in time " << selectrtc.getElapsedSeconds() << std::endl;

	uint64_t const numprephi = libmaus2::aio::PairFileDecoder::getLength(std::vector<std::string>(1,sasub));
	uint64_t const numprephiperthread = (numprephi + numthreads - 1)/numthreads;
	uint64_t const numprephipackages = numprephi ? ((numprephi + numprephiperthread - 1)/numprephiperthread) : 0;

	std::vector<std::string> phipairfn(numprephipackages);

	if ( verbose )
		std::cerr << "[V] generating Phi pairs (p1,r0,p0)" << std::endl;

	// generate phi pairs
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < numprephipackages; ++t )
	{
		uint64_t const low = t * numprephiperthread;
		uint64_t const high = std::min(low+numprephiperthread,numprephi);

		std::string const outfn = tmpgen.getFileName() + ".phipairfn";
		phipairfn[t] = outfn;

		libmaus2::huffman::LFPhiPairEncoder::unique_ptr_type lfenc(new libmaus2::huffman::LFPhiPairEncoder(outfn,4096));

		std::pair<uint64_t,uint64_t> Pprev;
		{
			libmaus2::aio::PairFileDecoder PFD(std::vector<std::string>(1,sasub));
			PFD.seekg((low + numprephi - 1)%numprephi);
			bool const ok = PFD.getNext(Pprev);
			assert ( ok );
		}

		libmaus2::aio::PairFileDecoder PFD(std::vector<std::string>(1,sasub));
		PFD.seekg(low);
		for ( uint64_t i = low; i < high; ++i )
		{
			std::pair<uint64_t,uint64_t> P;
			bool const ok = PFD.getNext(P);
			assert ( ok );

			if ( P.first == (Pprev.first + 1)%n )
			{
				// std::cerr << Pprev.first << "," << Pprev.second << "\t" << P.first << "," << P.second << std::endl;
				lfenc->encode(libmaus2::huffman::LFPhiPair(/*Pprev.first,*/Pprev.second,P.first,P.second));
			}

			Pprev = P;
		}

		lfenc->flush();
		lfenc.reset();
	}

	libmaus2::aio::FileRemoval::removeFile(histfn);
	libmaus2::aio::FileRemoval::removeFile(sasub);

	if ( verbose )
		std::cerr << "[V] sorting Phi pairs by p1" << std::endl;

	// sort phi pairs by component p1
	{
		struct LFPhiPairProjector
		{
			static uint64_t project(libmaus2::huffman::LFPhiPair const & S)
			{
				return S.p1;
			}
		};

		phipairfn = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFPhiPairDecoder,
			libmaus2::huffman::LFPhiPairEncoder,
			LFPhiPairProjector
		>(
			phipairfn,
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			libmaus2::huffman::LFPhiPairDecoder::getLength(phipairfn,numthreads), /* ihigh */
			(n-1), /* max sym */
			true /* max sym valid */
		);
	}

	struct LFPhiPairDecoderAdapter
	{
		typedef LFPhiPairDecoderAdapter this_type;
		typedef typename libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

		mutable libmaus2::huffman::LFPhiPairDecoder decoder;
		uint64_t const n;

		typedef libmaus2::util::ConstIterator<this_type,libmaus2::huffman::LFPhiPair> const_iterator;

		LFPhiPairDecoderAdapter(std::vector<std::string> const & Vfn, uint64_t const numthreads)
		: decoder(Vfn,0), n(libmaus2::huffman::LFPhiPairDecoder::getLength(Vfn,numthreads))
		{

		}

		libmaus2::huffman::LFPhiPair get(uint64_t const i) const
		{
			decoder.init(i);
			libmaus2::huffman::LFPhiPair P;
			bool const ok = decoder.decode(P);
			assert ( ok );
			return P;
		}

		libmaus2::huffman::LFPhiPair operator[](uint64_t const i) const
		{
			return get(i);
		}

		const_iterator begin() const
		{
			return const_iterator(this,0);
		}

		const_iterator end() const
		{
			return const_iterator(this,n);
		}

		struct P1Comp
		{
			bool operator()(libmaus2::huffman::LFPhiPair const & PA, libmaus2::huffman::LFPhiPair const & PB) const
			{
				return PA.p1 < PB.p1;
			}
		};

		uint64_t getOffset(uint64_t const p1) const
		{
			const_iterator it = ::std::lower_bound(begin(),end(),libmaus2::huffman::LFPhiPair(0,0,p1),P1Comp());
			return it - begin();
		}
	};

	if ( verbose )
	{
		typename LFPhiPairDecoderAdapter::unique_ptr_type Padp(new LFPhiPairDecoderAdapter(phipairfn,numthreads));
		std::cerr << "[V] computing LCP values via Karkkainen and Kempa algorithm for " << Padp->n << " ranks" << std::endl;
	}

	typedef typename input_types_type::linear_wrapper linear_wrapper;
	typedef typename input_types_type::circular_wrapper circular_wrapper;
	typedef typename input_types_type::base_input_stream::traits_type::char_type char_type;

	uint64_t const phiblocksize = (maxmem + sizeof(char_type) - 1)/sizeof(char_type);
	uint64_t const numphiblocks = (n+phiblocksize-1)/phiblocksize;
	std::string overflow_fn;
	uint64_t b_set = 0;

	std::string const philcpfn = tmpgen.getFileName(true) + ".philcp";
	libmaus2::huffman::LFPhiPairLCPEncoder::unique_ptr_type Pphilcpenc(new libmaus2::huffman::LFPhiPairLCPEncoder(philcpfn,4096));

	for ( uint64_t bi = 0; (bi < numphiblocks) || overflow_fn.size(); ++bi )
	{
		if ( verbose )
			std::cerr << "[V] computing LCP via Karkkainen and Kempa algorithm, block " << bi+1 << "/" << numphiblocks << std::endl;

		// block low
		uint64_t const b_low = bi * phiblocksize;
		// block high
		uint64_t const b_high = b_low + phiblocksize;

		typename LFPhiPairDecoderAdapter::unique_ptr_type Padp(new LFPhiPairDecoderAdapter(phipairfn,numthreads));
		uint64_t const off_low = Padp->getOffset(b_low);
		uint64_t const off_high = Padp->getOffset(b_high);
		uint64_t const off_range = off_high - off_low;

		struct LFPhiPairProjector
		{
			static uint64_t project(libmaus2::huffman::LFPhiPair const & S)
			{
				return S.p0;
			}
		};

		uint64_t const off_per_thread = (off_range + numthreads - 1)/numthreads;
		uint64_t const off_packs = off_per_thread ? ((off_range + off_per_thread -1)/off_per_thread) : 0;

		if ( verbose )
			std::cerr << "\t[V] extracting Phi range" << std::endl;

		std::vector<std::string> Vrange(off_packs);
		for ( uint64_t p = 0; p < off_packs; ++p )
		{
			std::string const outfn = tmpgen.getFileName(true) + ".Vrange";
			Vrange[p] = outfn;

			uint64_t const r_low = off_low + p * off_per_thread;
			uint64_t const r_high = std::min(r_low + off_per_thread, off_high);

			libmaus2::huffman::LFPhiPairEncoder::unique_ptr_type enc(new libmaus2::huffman::LFPhiPairEncoder(outfn,4096));
			libmaus2::huffman::LFPhiPairDecoder lfdec(phipairfn,r_low); // ZZZ

			for ( uint64_t i = r_low; i < r_high; ++i )
			{
				libmaus2::huffman::LFPhiPair P;
				bool const ok = lfdec.decode(P);
				assert ( ok );
				enc->encode(P);
			}

			enc->flush();
			enc.reset();
		}
		uint64_t ov_in = 0;
		if ( overflow_fn.size() )
		{
			Vrange.push_back(overflow_fn);
			ov_in = libmaus2::huffman::LFPhiPairDecoder::getLength(overflow_fn,numthreads);
		}

		if ( verbose )
			std::cerr << "\t[V] sorting Phi range" << std::endl;

		Vrange = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFPhiPairDecoder,
			libmaus2::huffman::LFPhiPairEncoder,
			LFPhiPairProjector
		>(
			Vrange,
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			libmaus2::huffman::LFPhiPairDecoder::getLength(Vrange,numthreads), /* ihigh */
			n-1, /* max sym */
			true /* max sym valid */
		);

		assert ( libmaus2::huffman::LFPhiPairDecoder::getLength(Vrange,numthreads) == (off_high-off_low) + ov_in );

		if ( verbose )
			std::cerr << "\t[V] loading text block" << std::endl;

		libmaus2::autoarray::AutoArray<char_type> A(b_high-b_low,false);

		typename circular_wrapper::unique_ptr_type circularBlockISI(new circular_wrapper(textfn,b_low));
		circularBlockISI->read(A.begin(),b_high-b_low);
		assert ( static_cast<int64_t>(circularBlockISI->gcount()) == static_cast<int64_t>(b_high-b_low) );
		circularBlockISI.reset();

		if ( verbose )
			std::cerr << "\t[V] processing block" << std::endl;

		// set up decoder
		libmaus2::huffman::LFPhiPairDecoder lfdec(Vrange,0);

		uint64_t l = 0;
		uint64_t prevp = 0;

		typename circular_wrapper::unique_ptr_type circularTextISI(new circular_wrapper(textfn,0));
		LinearAccessor<std::basic_istream<char_type> > LA(*circularTextISI);

		uint64_t overflow_cnt = 0;
		overflow_fn = tmpgen.getFileName(true) + ".phioverflow";
		libmaus2::huffman::LFPhiPairEncoder::unique_ptr_type Poverflow(new libmaus2::huffman::LFPhiPairEncoder(overflow_fn,4096));

		libmaus2::huffman::LFPhiPair P;
		while ( lfdec.decode(P) )
		{
			uint64_t const p0 = P.p0;
			uint64_t const p1 = P.p1;

			assert ( p0 >= prevp );
			l -= std::min(l,p0-prevp);

			if ( p1 < b_low )
				l = std::max(l,b_low - p1);

			while ( p1 + l < b_high && A[p1+l-b_low] == LA[p0+l] )
				++l;

			if ( p1+l >= b_high )
			{
				overflow_cnt += 1;
				Poverflow->encode(P);
			}
			else
			{
				Pphilcpenc->encode(libmaus2::huffman::LFPhiPairLCP(P.r1,P.p1,l));
				b_set += 1;
			}

			prevp = p0;
		}

		Poverflow->flush();
		Poverflow.reset();

		if ( ! overflow_cnt )
		{
			libmaus2::aio::FileRemoval::removeFile(overflow_fn);
			overflow_fn = std::string();
		}

		for ( uint64_t i = 0; i < Vrange.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Vrange[i]);
	}
	for ( uint64_t i = 0; i < phipairfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(phipairfn[i]);

	Pphilcpenc->flush();
	Pphilcpenc.reset();

	if ( verbose )
		std::cerr << "[V] sorting resulting data by p1" << std::endl;
	// sort by p1
	std::vector<std::string> Vphilcpfn;
	{
		struct LFPhiPairLCPProjector
		{
			static uint64_t project(libmaus2::huffman::LFPhiPairLCP const & S)
			{
				return S.p1;
			}
		};

		Vphilcpfn = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFPhiPairLCPDecoder,
			libmaus2::huffman::LFPhiPairLCPEncoder,
			LFPhiPairLCPProjector
		>(
			std::vector<std::string>(1,philcpfn),
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			libmaus2::huffman::LFPhiPairLCPDecoder::getLength(std::vector<std::string>(1,philcpfn),numthreads), /* ihigh */
			n-1, /* max sym */
			true /* max sym valid */
		);
	}

	uint64_t const numphilcp = libmaus2::huffman::LFPhiPairLCPDecoder::getLength(Vphilcpfn,numthreads);
	uint64_t const numphilcpperthread = (numphilcp + numthreads - 1)/numthreads;
	uint64_t const numphilcppackages = numphilcp ? ((numphilcp+numphilcpperthread-1)/numphilcpperthread) : 0;

	std::vector<std::string> Vranklcpfn(numphilcppackages);

	if ( verbose )
		std::cerr << "[V] computing differences p1,p0" << std::endl;
	// compute differences
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < numphilcppackages; ++t )
	{
		uint64_t const p_low = t * numphilcpperthread;
		uint64_t const p_high = std::min(p_low + numphilcpperthread,numphilcp);

		std::string const outfn = tmpgen.getFileName(true) + ".ranklcp";
		Vranklcpfn[t] = outfn;

		libmaus2::huffman::LFRankLCPEncoder::unique_ptr_type Pranklcpenc(new libmaus2::huffman::LFRankLCPEncoder(outfn,4096));

		libmaus2::huffman::LFPhiPairLCP Pprev;
		{
			libmaus2::huffman::LFPhiPairLCPDecoder dec(Vphilcpfn, (p_low + numphilcp - 1) % numphilcp );
			bool const ok = dec.decode(Pprev);
			assert ( ok );
		}

		libmaus2::huffman::LFPhiPairLCPDecoder dec(Vphilcpfn, p_low);
		for ( uint64_t i = p_low; i < p_high; ++i )
		{
			libmaus2::huffman::LFPhiPairLCP P;
			bool const ok = dec.decode(P);
			assert ( ok );

			if ( (Pprev.p1 + 1)%n == P.p1 )
			{
				assert ( P.lcp+1 >= Pprev.lcp );
				Pranklcpenc->encode(libmaus2::huffman::LFRankLCP(P.r1,P.lcp+1-Pprev.lcp));
			}

			Pprev = P;
		}

		Pranklcpenc->flush();
		Pranklcpenc.reset();
	}

	if ( verbose )
		std::cerr << "[V] sorting difference tuples by rank" << std::endl;
	// sort RankLCP files by rank
	{
		struct LFRankLCPProjector
		{
			static uint64_t project(libmaus2::huffman::LFRankLCP const & S)
			{
				return S.r;
			}
		};

		Vranklcpfn = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFRankLCPDecoder,
			libmaus2::huffman::LFRankLCPEncoder,
			LFRankLCPProjector
		>(
			Vranklcpfn,
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			libmaus2::huffman::LFRankLCPDecoder::getLength(Vranklcpfn,numthreads), /* ihigh */
			n-1, /* max sym */
			true /* max sym valid */
		);
	}

	struct LFRankLCPDecoderAdapter
	{
		typedef LFRankLCPDecoderAdapter this_type;
		typedef typename libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

		mutable libmaus2::huffman::LFRankLCPDecoder decoder;
		uint64_t const n;

		typedef libmaus2::util::ConstIterator<this_type,libmaus2::huffman::LFRankLCP> const_iterator;

		LFRankLCPDecoderAdapter(std::vector<std::string> const & Vfn, uint64_t const numthreads)
		: decoder(Vfn,0), n(libmaus2::huffman::LFRankLCPDecoder::getLength(Vfn,numthreads))
		{

		}

		libmaus2::huffman::LFRankLCP get(uint64_t const i) const
		{
			decoder.init(i);
			libmaus2::huffman::LFRankLCP P;
			bool const ok = decoder.decode(P);
			assert ( ok );
			return P;
		}

		libmaus2::huffman::LFRankLCP operator[](uint64_t const i) const
		{
			return get(i);
		}

		const_iterator begin() const
		{
			return const_iterator(this,0);
		}

		const_iterator end() const
		{
			return const_iterator(this,n);
		}

		struct RComp
		{
			bool operator()(libmaus2::huffman::LFRankLCP const & PA, libmaus2::huffman::LFRankLCP const & PB) const
			{
				return PA.r < PB.r;
			}
		};

		uint64_t getOffset(uint64_t const r) const
		{
			const_iterator it = ::std::lower_bound(begin(),end(),libmaus2::huffman::LFRankLCP(r,0),RComp());
			return it - begin();
		}
	};

	if ( verbose )
		std::cerr << "[V] filling differences into PD vector" << std::endl;
	// fill newly computed values into PD
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		typename LFRankLCPDecoderAdapter::unique_ptr_type adp(new LFRankLCPDecoderAdapter(Vranklcpfn,1));
		uint64_t const r_low = adp->getOffset(low);
		//uint64_t const r_high = adp->getOffset(high);

		assert ( r_low == adp->n || adp->get(r_low).r >= low );

		adp.reset();

		std::string const outfn = tmpgen.getFileName(true) + ".pdfntmp.1";
		pdfntmp[i] = outfn;

		libmaus2::gamma::GammaPDDecoder::unique_ptr_type PDin(new libmaus2::gamma::GammaPDDecoder(pdfn,low));
		libmaus2::gamma::GammaPDEncoder::unique_ptr_type PDout(new libmaus2::gamma::GammaPDEncoder(outfn,4096));

		uint64_t rnext = low;
		libmaus2::huffman::LFRankLCPDecoder lfdec(Vranklcpfn,r_low);

		libmaus2::huffman::LFRankLCP P;
		while ( lfdec.decode(P) && P.r < high )
		{
			while ( rnext < P.r )
			{
				uint64_t v;
				bool const ok = PDin->decode(v);
				assert ( ok );

				// use v
				PDout->encode(v);

				rnext += 1;
			}

			// decode unused value from PDin
			uint64_t v;
			bool const ok = PDin->decode(v);
			assert ( ok );

			// use P
			PDout->encode(P.lcp);

			rnext += 1;
		}

		while ( rnext < high )
		{
			uint64_t v;
			bool const ok = PDin->decode(v);
			assert ( ok );

			// use v
			PDout->encode(v);

			rnext += 1;
		}

		PDout->flush();
		PDout.reset();

		assert ( libmaus2::gamma::GammaPDDecoder::getNumValues(outfn) == (high-low) );
	}

	for ( uint64_t i = 0; i < pdfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(pdfn[i]);
	pdfn = pdfntmp;

	libmaus2::aio::FileRemoval::removeFile(philcpfn);
	for ( uint64_t t = 0; t < Vphilcpfn.size(); ++t )
		libmaus2::aio::FileRemoval::removeFile(Vphilcpfn[t]);
	for ( uint64_t t = 0; t < Vranklcpfn.size(); ++t )
		libmaus2::aio::FileRemoval::removeFile(Vranklcpfn[t]);

	struct AlphabetSymbolFreqComparator
	{
		bool operator()(AlphabetSymbol const & A, AlphabetSymbol const & B) const
		{
			return A.freq > B.freq;
		}
	};
	struct AlphabetSymbolSymComparator
	{
		bool operator()(AlphabetSymbol const & A, AlphabetSymbol const & B) const
		{
			return A.sym < B.sym;
		}
	};
	std::vector < AlphabetSymbol > Valphabet;
	uint64_t cid = 0;
	uint64_t freqsum = 0;
	for ( std::map<int64_t,uint64_t>::const_iterator ita = histmap.begin(); ita != histmap.end(); ++ita )
	{
		Valphabet.push_back(AlphabetSymbol(ita->first,cid++,ita->second,false));
		freqsum += ita->second;
	}
	//std::sort(Valphabet.begin(),Valphabet.end(),AlphabetSymbolFreqComparator());
	// std::sort(Valphabet.begin(),Valphabet.end(),AlphabetSymbolSymComparator());
	uint64_t const tbuckets = 32;

	uint64_t const bucketsplit = computeSplit(tbuckets, Valphabet);
	if ( verbose )
		std::cerr << "[V] bucketsplit=" << bucketsplit << std::endl;

	uint64_t alow = 0;
	uint64_t bucketid = 0;
	std::vector<bool> bucketfrequent;
	std::vector<uint64_t> Vbucketsize;
	while ( alow < Valphabet.size() )
	{
		uint64_t ahigh = alow;
		uint64_t s = Valphabet[ahigh++].freq;

		while ( ahigh < Valphabet.size() && s < bucketsplit )
			s += Valphabet[ahigh++].freq;

		Vbucketsize.push_back(ahigh-alow);
		bool const frequent = (ahigh-alow)==1;
		for ( uint64_t i = alow; i < ahigh; ++i )
		{
			Valphabet[i].id = bucketid;
			Valphabet[i].subid = i-alow;
			Valphabet[i].frequent = frequent;
		}
		bucketid++;

		bucketfrequent.push_back(frequent);

		alow = ahigh;
	}

	if ( verbose )
		std::cerr << "[V] generating SA pairs from sampled inverse suffix array" << std::endl;
	// generate pair SA
	std::string const pairisa = tmpgen.getFileName(true) + ".pairisa";
	uint64_t isasamplingrate = 0;
	{
		libmaus2::aio::InputStreamInstance ISAIN(isafn);
		libmaus2::aio::SynchronousGenericInput<uint64_t> ISASGI(isafn,4096);

		bool const rateok = ISASGI.getNext(isasamplingrate);
		assert ( rateok );

		uint64_t numisasamples;
		bool const numisasamplesok = ISASGI.getNext(numisasamples);
		assert ( numisasamplesok );

		libmaus2::aio::OutputStreamInstance RISAOUT(pairisa);
		libmaus2::aio::SynchronousGenericOutput<uint64_t> RISASGO(RISAOUT,4096);

		for ( uint64_t ri = 0; ri < numisasamples; ++ri )
		{
			uint64_t v;
			bool const ok = ISASGI.getNext(v);
			assert ( ok );
			RISASGO.put(v /* rank */);
			RISASGO.put(ri * isasamplingrate);
		}

		RISASGO.flush();
		RISAOUT.flush();
	}

	if ( verbose )
		std::cerr << "[V] sorting SA pairs by rank" << std::endl;
	// sort pair SA by rank
	std::string const resortedisa = tmpgen.getFileName(true) + ".resortedisa";
	{
		std::string const tmp = tmpgen.getFileName(true) + ".resortedisatmp";

		libmaus2::sorting::PairFileSorting::sortPairFile(
			std::vector<std::string>(1,pairisa),
			tmp,
			false /* sort by first component */,
			true /* keep first */,
			true /* keep second */,
			resortedisa,
			maxmem,
			true /* parallel */,
			true /* delete input */
		);

		libmaus2::aio::FileRemoval::removeFile(pairisa);
		libmaus2::aio::FileRemoval::removeFile(tmp);
	}

	uint64_t const resortedisasize = libmaus2::util::GetFileSize::getFileSize(resortedisa);
	assert ( resortedisasize % (2*sizeof(uint64_t)) == 0 );
	uint64_t const numisasamples = resortedisasize / (2*sizeof(uint64_t));
	uint64_t const lfpacksize = (numisasamples + numthreads - 1)/numthreads;
	uint64_t const lfpacks = numisasamples ? ((numisasamples + lfpacksize - 1)/lfpacksize) : 0;

	std::vector<std::string> Vlfrankpos(lfpacks);

	if ( verbose )
		std::cerr << "[V] converting SA pairs to RankPos" << std::endl;
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < lfpacks; ++t )
	{
		std::string const lfrankpos = tmpgen.getFileName(true) + ".lfrankpos";
		Vlfrankpos[t] = lfrankpos;
		libmaus2::huffman::LFRankPosEncoder::unique_ptr_type Plfrankposenc(new libmaus2::huffman::LFRankPosEncoder(lfrankpos,16*1024));
		uint64_t const low = t * lfpacksize;
		uint64_t const high = std::min(low+lfpacksize,numisasamples);
		libmaus2::aio::InputStreamInstance ISI(resortedisa);
		ISI.clear();
		ISI.seekg(low * 2 * sizeof(uint64_t));
		libmaus2::aio::SynchronousGenericInput<uint64_t> SGI(ISI,4096);

		for ( uint64_t i = low; i < high; ++i )
		{
			uint64_t r, p;
			bool const rok = SGI.getNext(r);
			assert ( rok );
			bool const pok = SGI.getNext(p);
			assert ( pok );
			uint64_t * v = 0;

			if ( SA )
			{
				assert ( SA[r] == p );
			}

			Plfrankposenc->encode(libmaus2::huffman::LFRankPos(r,p,0,v,true /* active */));
		}

		Plfrankposenc->flush();
		Plfrankposenc.reset();
	}

	libmaus2::aio::FileRemoval::removeFile(resortedisa);

	if ( verbose )
		std::cerr << "[V] computing thread block wise symbol histogram" << std::endl;
	libmaus2::autoarray::AutoArray < uint64_t > Gsymhist((threadpacks+1) * (maxsym+1)+1, false);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(threadpacks)
	#endif
	for ( uint64_t i = 0; i < threadpacks; ++i )
	{
		uint64_t * hist = Gsymhist.begin() + i * (maxsym+1);
		std::fill(hist,hist + (maxsym+1), 0ull);

		uint64_t const low = i*threadpacksize;
		uint64_t const high = std::min(low+threadpacksize,n);

		uint64_t todo = high-low;

		libmaus2::huffman::RLDecoder rldec(BWTin,low/* offset */,1);
		libmaus2::huffman::RLDecoder::run_type R;

		while ( todo )
		{
			bool const ok = rldec.decodeRun(R);
			assert ( ok );
			uint64_t const av = std::min(R.rlen,todo);
			hist [ R.sym ] += av;
			todo -= av;
		}
	}

	std::fill(
		Gsymhist.begin() + (threadpacks+0) * (maxsym+1) + 0,
		Gsymhist.begin() + (threadpacks+1) * (maxsym+1) + 1,
		0ull
	);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < maxsym+1; ++i )
	{
		libmaus2::util::PrefixSums::prefixSums(
			Gsymhist.begin() + i,
			Gsymhist.begin() + (threadpacks+1) * (maxsym+1) + i,
			(maxsym+1)
		);
	}

	libmaus2::util::PrefixSums::prefixSums(Gsymhist.begin() + (threadpacks+0) * (maxsym+1) + 0,Gsymhist.begin() + (threadpacks+1) * (maxsym+1) + 1);
	assert ( Gsymhist.end()[-1] == n );

	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < maxsym+1; ++i )
	{
		uint64_t const add = Gsymhist[threadpacks * (maxsym+1) + i];
		for ( uint64_t j = 0; j < threadpacks; ++j )
			Gsymhist[ j * (maxsym+1) + i ] += add;
	}

	if ( verbose )
		std::cerr << "[V] filling PD values into position order" << std::endl;

	libmaus2::timing::RealTimeClock newrtc; newrtc.start();
	for ( uint64_t round = 0; round < isasamplingrate+1; ++round )
	{
		if ( verbose )
			std::cerr << "\t[V] round " << round+1 << "//" << isasamplingrate+1 << std::endl;

		libmaus2::autoarray::AutoArray<uint64_t> symhist(Gsymhist.size(),false);
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < symhist.size(); ++i )
			symhist[i] = Gsymhist[i];

		std::vector < std::string > Voutfn(bucketid * threadpacks);
		libmaus2::parallel::PosixSpinLock Vbucketnctlock;
		libmaus2::autoarray::AutoArray < uint64_t > Vbucketcnt(bucketid);

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(threadpacks)
		#endif
		for ( uint64_t t = 0; t < threadpacks; ++t )
		{
			uint64_t * H = symhist.begin() + t * (maxsym+1);

			libmaus2::autoarray::AutoArray < libmaus2::huffman::LFRankPosEncoder::unique_ptr_type > Aenc(bucketid);
			libmaus2::autoarray::AutoArray < libmaus2::huffman::LFSymRankPosEncoder::unique_ptr_type > Aresenc(bucketid);
			for ( uint64_t i = 0; i < bucketid; ++i )
			{
				std::string const fn = tmpgen.getFileName(true) + "lfrankpos";
				// Lfn [ i ] = fn;
				// bucket 0 thread 0, bucket 0 thread 1 ...
				Voutfn [ i * threadpacks + t ] = fn;

				if ( bucketfrequent[i] )
				{
					libmaus2::huffman::LFRankPosEncoder::unique_ptr_type Tenc(new libmaus2::huffman::LFRankPosEncoder(fn,16*1024));
					Aenc[i] = UNIQUE_PTR_MOVE(Tenc);
				}
				else
				{
					libmaus2::huffman::LFSymRankPosEncoder::unique_ptr_type Tresenc(new libmaus2::huffman::LFSymRankPosEncoder(fn,16*1024));
					Aresenc[i] = UNIQUE_PTR_MOVE(Tresenc);
				}
			}

			uint64_t const low = t*threadpacksize;
			uint64_t const high = std::min(low+threadpacksize,n);

			libmaus2::huffman::RLDecoder rldec(BWTin,low/* offset */,1);
			libmaus2::huffman::RLDecoder::run_type R;

			libmaus2::huffman::LFRankPosDecoder rpdec(Vlfrankpos,low,libmaus2::huffman::LFRankPosDecoder::init_type_rank);

			libmaus2::huffman::LFRankPos P;
			uint64_t rr = low;

			std::vector < uint64_t > bucketcnt(bucketid);

			if ( round == 0 )
			{
				while ( rpdec.decode(P) && P.r < high )
				{
					while ( rr < P.r )
					{
						bool const ok = rldec.decodeRun(R);
						assert ( ok );
						uint64_t const av = std::min(P.r-rr,R.rlen);

						H [ R.sym ] += av;
						rr += av;
						R.rlen -= av;

						if ( R.rlen )
							rldec.putBack(R);
					}

					assert ( rr == P.r );

					int64_t const sym = rldec.decode();
					if ( P.active )
						P.p = (P.p + n - 1) % n;
					P.r = H [ sym ]++;

					AlphabetSymbol const & A = Valphabet[sym];
					if ( A.frequent )
						Aenc[A.id]->encode(P);
					else
					{
						Aresenc[A.id]->encode(libmaus2::huffman::LFSymRankPos(A.subid,P.r,P.p,P.n,P.v,P.active));
						bucketcnt[A.id]++;
					}

					if ( SA )
					{
						assert ( SA[P.r] == P.p );
					}

					rr += 1;
				}
			}
			else
			{
				libmaus2::gamma::GammaPDDecoder PDin(pdfn,low);

				while ( rpdec.decode(P) && P.r < high )
				{
					while ( rr < P.r )
					{
						bool const ok = rldec.decodeRun(R);
						assert ( ok );
						uint64_t const av = std::min(P.r-rr,R.rlen);

						H [ R.sym ] += av;
						rr += av;
						R.rlen -= av;

						for ( uint64_t z = 0; z < av; ++z )
							PDin.decode();

						if ( R.rlen )
							rldec.putBack(R);
					}

					assert ( rr == P.r );

					bool const activein = P.active;
					bool const activeout = activein && ((P.p % isasamplingrate) != 0);

					uint64_t const pdv = PDin.decode();
					int64_t const sym = rldec.decode();
					if ( activeout )
						P.p = (P.p + n - 1) % n;
					P.r = H [ sym ]++;

					AlphabetSymbol const & A = Valphabet[sym];
					if ( A.frequent )
					{
						if ( activein )
							Aenc[A.id]->encode(libmaus2::huffman::LFRankPos(P.r,P.p,P.n,P.v,activeout),pdv);
						else
							Aenc[A.id]->encode(libmaus2::huffman::LFRankPos(P.r,P.p,P.n,P.v,activeout));
					}
					else
					{
						if ( activein )
							Aresenc[A.id]->encode(libmaus2::huffman::LFSymRankPos(A.subid,P.r,P.p,P.n,P.v,activeout),pdv);
						else
							Aresenc[A.id]->encode(libmaus2::huffman::LFSymRankPos(A.subid,P.r,P.p,P.n,P.v,activeout));
						bucketcnt[A.id]++;
					}

					if ( SA )
					{
						bool const ok = (!activeout) || (SA[P.r] == P.p);

						if ( ! ok )
							std::cerr << "P.r=" << P.r << " P.p=" << P.p << " activeout=" << activeout << std::endl;
						assert ( ok );
					}

					rr += 1;
				}
			}

			for ( uint64_t i = 0; i < bucketid; ++i )
				if ( bucketfrequent[i] )
				{
					Aenc[i]->flush();
					Aenc[i].reset();
				}
				else
				{
					Aresenc[i]->flush();
					Aresenc[i].reset();

					libmaus2::parallel::ScopePosixSpinLock slock(Vbucketnctlock);
					Vbucketcnt[i] += bucketcnt[i];
				}
		}

		for ( uint64_t i = 0; i < bucketid; ++i )
			if ( ! bucketfrequent[i] )
			{
				struct Projector
				{
					static uint64_t project(libmaus2::huffman::LFSymRankPos const & P)
					{
						return P.sym;
					}
				};

				std::vector<std::string> Rfnsorted = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
					libmaus2::huffman::LFSymRankPosDecoder,
					libmaus2::huffman::LFSymRankPosEncoder,
					Projector
				>
				(
					std::vector < std::string >(Voutfn.begin() + i * threadpacks,Voutfn.begin() + (i+1) * threadpacks),
					numthreads /* numthreads */,
					4096 /* maxfiles */,
					true /* delete input */,
					tmpgen,
					16*1024 /* bs */,
					0 /* low */,
					Vbucketcnt[i] /* high */,
					Vbucketsize[i]-1 /* maxsym */,
					true /* maxsymvalid */
				);

				uint64_t const packsize = (Vbucketcnt[i] + numthreads - 1)/numthreads;
				uint64_t const numpacks = Vbucketcnt[i] ? ((Vbucketcnt[i] + packsize - 1)/packsize) : 0;

				assert ( numpacks <= threadpacks );

				for ( uint64_t t = 0; t < threadpacks; ++t )
				{
					uint64_t const low = std::min(t * packsize,Vbucketcnt[i]);
					uint64_t const high = std::min(low+packsize,Vbucketcnt[i]);

					std::string const fn = tmpgen.getFileName(true) + ".resort";
					Voutfn [ i * threadpacks + t ] = fn;

					libmaus2::huffman::LFRankPosEncoder::unique_ptr_type Aenc(new libmaus2::huffman::LFRankPosEncoder(fn,16*1024));

					libmaus2::huffman::LFSymRankPosDecoder dec(Rfnsorted,low);
					libmaus2::huffman::LFSymRankPos P;
					for ( uint64_t j = low; j < high; ++j )
					{
						bool const ok = dec.decode(P);
						assert ( ok );
						Aenc->encode(libmaus2::huffman::LFRankPos(P.r,P.p,P.n,P.v,P.active));
					}

					Aenc->flush();
					Aenc.reset();
				}

				for ( uint64_t i = 0; i < Rfnsorted.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile(Rfnsorted[i]);
			}

		for ( uint64_t i = 0; i < Vlfrankpos.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Vlfrankpos[i]);
		Vlfrankpos = Voutfn;
	}

	if ( verbose )
		std::cerr << "[V] sorting resulting data tuples by position" << std::endl;

	{
		struct LFRankPosProjector
		{
			static uint64_t project(libmaus2::huffman::LFRankPos const & P)
			{
				return P.p;
			}
		};

		libmaus2::timing::RealTimeClock fillsortrtc; fillsortrtc.start();
		Vlfrankpos = libmaus2::sorting::ParallelExternalRadixSort::parallelRadixSort<
			libmaus2::huffman::LFRankPosDecoder,
			libmaus2::huffman::LFRankPosEncoder,
			LFRankPosProjector
		>(
			Vlfrankpos,
			numthreads,
			4096 /* max files */,
			true /* delete input */,
			tmpgen,
			4096 /* output block size */,
			0, /* ilow */
			libmaus2::huffman::LFRankPosDecoder::getLength(Vlfrankpos,numthreads), /* ihigh */
			n-1, /* max sym */
			true /* max sym valid */
		);

		if ( verbose )
			std::cerr << "\t[V] sorted in time " << fillsortrtc.getElapsedSeconds() << std::endl;
	}

	if ( verbose )
		std::cerr << "[V] filled and sorted in time " << newrtc.getElapsedSeconds() << std::endl;

	if ( verbose )
		std::cerr << "[V] producing succinct LCP bit vector" << std::endl;

	uint64_t const numentries = libmaus2::huffman::LFRankPosDecoder::getLength(Vlfrankpos,numthreads);
	uint64_t const entriesperthread = (numentries + numthreads - 1)/numthreads;
	uint64_t const entrypacks = numentries ? ((numentries + entriesperthread - 1)/entriesperthread) : 0;

	std::vector < std::string > Vbvfn(threadpacks);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads)
	#endif
	for ( uint64_t t = 0; t < entrypacks; ++t )
	{
		uint64_t const low = t * entriesperthread;
		uint64_t const high = std::min(low+entriesperthread,numentries);

		std::string const outfn = tmpgen.getFileName(true) + ".lcpbit";
		Vbvfn[t] = outfn;

		libmaus2::bitio::BitVectorOutput::unique_ptr_type BVout(new libmaus2::bitio::BitVectorOutput(outfn));

		libmaus2::huffman::LFRankPosDecoder dec(Vlfrankpos,low);
		libmaus2::huffman::LFRankPos P;

		for ( uint64_t i = low; i < high; ++i )
		{
			bool const ok = dec.decode(P);
			assert ( ok );

			for ( uint64_t j = 0; j < P.n; ++j )
			{
				uint64_t const v = P.v[j];
				for ( uint64_t k = 0; k < v; ++k )
					BVout->writeBit(false);
				BVout->writeBit(true);
			}
		}

		BVout->flush();
		BVout.reset();
	}

	for ( uint64_t i = 0; i < Vlfrankpos.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Vlfrankpos[i]);

	for ( uint64_t i = 0; i < pdfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(pdfn[i]);
	// libmaus2::aio::FileRemoval::removeFile(pairisa);
	if ( deleteBWTin )
		for ( uint64_t i = 0; i < BWTin.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(BWTin[i]);

	if ( verbose )
	{
		std::vector<std::string> memfn;
		libmaus2::aio::MemoryFileContainer::list(memfn);
		for ( uint64_t i = 0; i < memfn.size(); ++i )
			std::cerr << "left over " << memfn[i] << std::endl;
	}

	return Vbvfn;
}

// test for string s. String is assumed to have a unique minimal terminator symbol
void testn(std::string const & s, uint64_t const numthreads, uint64_t const maxrounds)
{
	uint64_t const verbthres = 1024;
	uint64_t const n = s.size();
	bool const verbose = n >= verbthres;

	std::string const tmpprefix = "mem://lcpbit_";
	libmaus2::util::TempFileNameGenerator tmpgen(tmpprefix,3);

	// compute suffix array
	if ( verbose )
		std::cerr << "[V] computing suffix array via divsufsort...";
	typedef libmaus2::suffixsort::DivSufSort<32,std::string::iterator,std::string::const_iterator,int32_t *, int32_t const *,256,false> sort_type;
	libmaus2::autoarray::AutoArray<int32_t> SA(n,false);
	sort_type::divsufsort(s.begin(),SA.begin(),n);
	if ( verbose )
		std::cerr << "done." << std::endl;

	// compute BWT
	if ( verbose )
		std::cerr << "[V] computing BWT from text and SA...";
	libmaus2::autoarray::AutoArray<char> BWT(n,false);
	for ( size_t i = 0; i < SA.size(); ++i )
		BWT[i] = s [ ((SA[i]+n)-1) % n ];
	if ( verbose )
		std::cerr << "done." << std::endl;

	// write BWT to memory file
	if ( verbose )
		std::cerr << "[V] writing BWT to file...";
	std::string const BWTfn = tmpgen.getFileName(true) + ".bwt";
	libmaus2::huffman::RLEncoderStd::unique_ptr_type Penc(new libmaus2::huffman::RLEncoderStd(BWTfn,n,8192));
	for ( uint64_t i = 0; i < n; ++i )
		Penc->encode(BWT[i]);
	Penc->flush();
	Penc.reset();
	if ( verbose )
		std::cerr << "done." << std::endl;

	// compute ISA
	if ( verbose )
		std::cerr << "[V] computing ISA from text and SA...";
	libmaus2::autoarray::AutoArray<int32_t> ISA(n,false);
	for ( uint64_t i = 0; i < n; ++i )
		ISA[SA[i]] = i;
	if ( verbose )
		std::cerr << "done." << std::endl;

	// write sampled inverse suffix array
	if ( verbose )
		std::cerr << "[V] writing sampled ISA to memory file...";
	std::string const isafn = tmpgen.getFileName(true) + ".isafn";
	{
		uint64_t const isasamplingrate = 8;
		uint64_t const numisasamples = (n + isasamplingrate - 1)/isasamplingrate;

		libmaus2::aio::OutputStreamInstance OSI(isafn);
		libmaus2::aio::SynchronousGenericOutput<uint64_t> SGOSA(OSI,4096);

		SGOSA.put(isasamplingrate);
		SGOSA.put(numisasamples);

		for ( uint64_t i = 0; i < n; i += isasamplingrate )
			SGOSA.put(ISA[i]);
		SGOSA.flush();
		OSI.flush();
	}
	if ( verbose )
		std::cerr << "done." << std::endl;

	if ( verbose )
		std::cerr << "[V] writing text to memory file...";
	std::string const textfn = tmpgen.getFileName(true) + ".text";
	libmaus2::aio::OutputStreamInstance::unique_ptr_type textOSI(new libmaus2::aio::OutputStreamInstance(textfn));
	textOSI->write(s.c_str(),n);
	textOSI->flush();
	textOSI.reset();
	if ( verbose )
		std::cerr << "done." << std::endl;

	typedef libmaus2::suffixsort::ByteInputTypes input_types_type;
	bool const bwtdeleteable = true;
	std::vector<std::string> const Vbvfn = computeSuccinctLCP<input_types_type>(BWTfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,16*1024*1024,SA.get());

	libmaus2::aio::FileRemoval::removeFile(BWTfn);
	libmaus2::aio::FileRemoval::removeFile(isafn);
	libmaus2::aio::FileRemoval::removeFile(textfn);

	if ( verbose )
		std::cerr << "[V] computing LCP array via Kasai et al algorithm...";
	libmaus2::autoarray::AutoArray<int32_t> LCP(n,false);
	libmaus2::lcp::computeLCPKasai(s.begin(),n,SA.begin(),LCP.begin());
	if ( verbose )
		std::cerr << "done." << std::endl;

	if ( verbose )
		std::cerr << "[V] checking succinct LCP...";
	{
		libmaus2::bitio::BitVectorInput BVin(Vbvfn);
		for ( uint64_t i = 0; i < n; ++i )
		{
			uint64_t v = 0;
			while ( !BVin.readBit() )
				++v;

			uint64_t p = i;
			uint64_t p1 = (i+n-1)%n;

			uint64_t r = ISA[p];
			uint64_t r1 = ISA[p1];

			uint64_t l0 = LCP[r];
			uint64_t l1 = LCP[r1];

			uint64_t e = (l0+1)-l1;

			// std::cerr << "p=" << p << " v=" << v << " l0=" << l0 << " l1=" << l1 << " e=" << e << std::endl;

			assert ( v == e );
		}
	}
	if ( verbose )
		std::cerr << "done." << std::endl;

	for ( uint64_t i = 0; i < Vbvfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Vbvfn[i]);
}

// test for all strings of length k with alphabet of albits bits width
void testnk(unsigned int const l, unsigned int const albits)
{
	std::cerr << "testing all strings of length l=" << l << " with albits=" << albits << "...";
	for ( uint64_t z = 0; z < (1ull<<(albits * l)); ++z )
	{
		std::string s(l+1,'#');
		for ( unsigned int i = 0; i < l; ++i )
			s[i] = ((z>>(i*albits))&((1ull<<albits)-1))+'a';
		testn(s,1 /* numthreads */,std::numeric_limits<uint64_t>::max() /* max rounds */);

		// std::cerr << s << std::endl;

		if ( (z & (1024-1)) == 0 )
			std::cerr << "(" << static_cast<double>(z) / (1ull<<(albits * l)) << ")";
	}

	std::cerr << "(1)" << std::endl;
}

// test for a random string of length l with alphabet width al bits
void testrandomn(uint64_t const l, unsigned int const al, uint64_t const numthreads, uint64_t const maxrounds)
{
	if ( l > 1024 )
		std::cerr << "[V] generating " << l << "," << al << "...";
	std::string s(l+1,'#');
	for ( uint64_t i = 0; i < l; ++i )
		s[i] = 'A' + (libmaus2::random::Random::rand64() % al);
	if ( l > 1024 )
		std::cerr << "done." << std::endl;

	if ( l > 1024 )
		std::cerr << "[V] running test for string of length " << l << std::endl;
	testn(s,numthreads,maxrounds);
	if ( l > 1024 )
		std::cerr << "[V] running test for string of length " << l << " done.\n";
}

// test for XZ compressed input file
void testnXz(std::string const & fn, uint64_t const numthreads, uint64_t maxrounds)
{
	if ( libmaus2::util::GetFileSize::fileExists(fn) )
	{
		std::cerr << "[V] loading " << fn << "...";
		std::string const sec = readXzFile(fn);
		std::cerr << "done." << std::endl;
		testn(sec,numthreads,maxrounds);
	}
}

// test for plain file (assumed not to contain the \0 symbol)
void testnPlain(std::string const & fn, uint64_t const numthreads, uint64_t maxrounds)
{
	if ( libmaus2::util::GetFileSize::fileExists(fn) )
	{
		std::cerr << "[V] loading " << fn << "...";
		std::string const sec = readPlainFile(fn);
		std::cerr << "done." << std::endl;
		testn(sec,numthreads,maxrounds);
	}
}

#include <libmaus2/aio/SizeMonitorThread.hpp>

void computeSuccinctLCP(libmaus2::util::ArgParser const & arg)
{
	std::string const inputtype = arg.uniqueArgPresent("i") ? arg["i"] : "bytestream";
	libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type itype =
		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::parseInputType(inputtype);
	std::string const tmpprefix = arg.uniqueArgPresent("T") ? arg["T"] : "tmp_prefix";
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t const maxrounds = arg.uniqueArgPresent("R") ? arg.getUnsignedNumericArg<uint64_t>("R") : 64;
	uint64_t maxmem = arg.uniqueArgPresent("M") ? arg.getUnsignedNumericArg<uint64_t>("M") : (1024ull*1024ull);

	std::string const bwtfn = arg[0];
	std::string const isafn = arg[1];
	std::string const textfn = arg[2];

	std::string const outputfn = arg.uniqueArgPresent("o") ? arg["o"] : (libmaus2::util::OutputFileNameTools::clipOff(bwtfn,".bwt")+".lcpbit");

	libmaus2::util::TempFileNameGenerator tmpgen(tmpprefix,5);

	libmaus2::aio::SizeMonitorThread SMT(tmpprefix,5/*sec*/,&std::cerr);
	SMT.start();

	std::string const tmpbwtfn = tmpgen.getFileName(true);

	{
		libmaus2::aio::OutputStreamInstance OSI(tmpbwtfn);
		libmaus2::aio::InputStreamInstance ISI(bwtfn);
		libmaus2::util::GetFileSize::copy(ISI,OSI,libmaus2::util::GetFileSize::getFileSize(ISI));
	}

	std::vector<std::string> resfn;
	bool const bwtdeleteable = true;

	switch ( itype )
	{
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_bytestream:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::ByteInputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_compactstream:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::CompactInputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_pac:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::PacInputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_pacterm:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::PacTermInputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_lz4:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::Lz4InputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		case libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::bwt_merge_input_type_utf_8:
			resfn = computeSuccinctLCP<libmaus2::suffixsort::Utf8InputTypes>(tmpbwtfn,bwtdeleteable,isafn,textfn,tmpgen,numthreads,maxrounds,maxmem);
			break;
		default:
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt: unknown/unsupported input type" << std::endl;
			lme.finish();
			throw lme;
		}
		break;
	}

	libmaus2::aio::FileRemoval::removeFile(tmpbwtfn);

	std::cerr << "[V] concatenating bit vectors to " << outputfn << std::endl;

	libmaus2::bitio::BitVectorOutput::unique_ptr_type BVO(new libmaus2::bitio::BitVectorOutput(outputfn));
	for ( uint64_t i = 0; i < resfn.size(); ++i )
	{
		libmaus2::bitio::BitVectorInput::unique_ptr_type BVI(new libmaus2::bitio::BitVectorInput(std::vector<std::string>(1,resfn[i]),0));
		uint64_t subn = BVI->bitsleft;

		while ( subn-- )
			BVO->writeBit(BVI->readBit());

		BVI.reset();

		libmaus2::aio::FileRemoval::removeFile(resfn[i]);
	}
	BVO->flush();
	BVO.reset();

	SMT.printSize(std::cerr);
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		bool const testmode = arg.uniqueArgPresent("test");

		if ( testmode )
		{
			testnk(6,1);

			uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
			uint64_t const maxrounds = arg.uniqueArgPresent("R") ? arg.getUnsignedNumericArg<uint64_t>("R") : 64;
			//uint64_t maxmem = arg.uniqueArgPresent("M") ? arg.getUnsignedNumericArg<uint64_t>("M") : (1024ull*1024ull);

			testn("abbab#",numthreads,maxrounds);

			testrandomn(1024,8,numthreads,maxrounds);
			testrandomn(2048,8,numthreads,maxrounds);
			testrandomn(16*1024,8,numthreads,maxrounds);
			testrandomn(1024*1024,8,numthreads,maxrounds);
			testrandomn(1024*1024*16,8,numthreads,maxrounds);
			testrandomn(1024*1024*128,8,numthreads,maxrounds);

			testnPlain("configure",numthreads,maxrounds);

			testnXz("testdata/hg19_000000.xz",numthreads,maxrounds);
			testnXz("testdata/dmel_test.xz",numthreads,maxrounds);
			testnXz("testdata/ecoli_test.xz",numthreads,maxrounds);

			testnk(6,3);
			testnk(8,2);
		}
		else
		{
			computeSuccinctLCP(arg);
		}

		return 0;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
