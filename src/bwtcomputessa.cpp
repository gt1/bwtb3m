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
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/aio/InputStreamFactoryContainer.hpp>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/util/iterator.hpp>
#include <libmaus2/huffman/RLDecoder.hpp>
#include <libmaus2/util/NumberMapSerialisation.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/bitio/CompactDecoderBuffer.hpp>
#include <libmaus2/aio/ConcatInputStream.hpp>
#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>

static std::vector<uint64_t> loadHMap(std::string const & hist)
{
	libmaus2::aio::InputStream::unique_ptr_type Phist(libmaus2::aio::InputStreamFactoryContainer::constructUnique(hist));
	std::map<int64_t,uint64_t> const hmap = ::libmaus2::util::NumberMapSerialisation::deserialiseMap<libmaus2::aio::InputStream,int64_t,uint64_t>(*Phist);
	std::vector<uint64_t> D;

	if ( hmap.size() && hmap.begin()->first >= 0 )
	{
		int64_t const topsym = hmap.rbegin()->first;
		D.resize(topsym+2);
		for ( std::map<int64_t,uint64_t>::const_iterator ita = hmap.begin(); ita != hmap.end(); ++ita )
		{
			D [ ita->first ] = ita->second;
			// std::cerr << "[HIN]\t" << ita->first << "\t" << ita->second << std::endl;
		}
		D [ D.size()-1 ] = 0;

		// compute prefix sums
		uint64_t s = 0;
		for ( uint64_t i = 0; i < D.size(); ++i )
		{
			uint64_t const t = D[i];
			D[i] = s;
			s += t;
		}

		#if 0
		for ( uint64_t i = 0; i < D.size(); ++i )
			std::cerr << "[HIST]\t" << i << "\t" << D[i] << std::endl;
		#endif
	}

	return D;
}


struct WordPairAccessor
{
	libmaus2::aio::InputStream::unique_ptr_type Pin;
	std::istream & in;
	uint64_t const n;

	static uint64_t getNumberOfEntries(std::istream & in)
	{
		in.seekg(0,std::ios::end);
		uint64_t const p = in.tellg();
		if ( (p % (2*sizeof(uint64_t))) != 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "WordPairAccessor::init(): File size " << p << " is not a multiple of " << (2*sizeof(uint64_t)) << std::endl;
			lme.finish();
			throw lme;
		}

		return p / (2*sizeof(uint64_t));
	}

	WordPairAccessor(std::string const & sortedisa)
	: Pin(libmaus2::aio::InputStreamFactoryContainer::constructUnique(sortedisa)), in(*Pin), n(getNumberOfEntries(in))
	{

	}

	WordPairAccessor(std::istream & rin) : in(rin), n(getNumberOfEntries(in))
	{

	}

	std::pair<uint64_t,uint64_t> get(uint64_t const i) const
	{
		if ( i >= n )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "WordPairAccessor::get(): index " << i << " is out of range [" << 0 << "," << i << ")" << std::endl;
			lme.finish();
			throw lme;
		}

		in.clear();
		in.seekg(i*(2*sizeof(uint64_t)));
		uint64_t const v0 = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		uint64_t const v1 = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		return std::pair<uint64_t,uint64_t>(v0,v1);
	}

	std::pair<uint64_t,uint64_t> operator[](uint64_t const i) const
	{
		return get(i);
	}

	uint64_t size() const
	{
		return n;
	}
};

struct WordPairAccessorFirstAdapter
{
	WordPairAccessor * acc;

	WordPairAccessorFirstAdapter(WordPairAccessor & racc) : acc(&racc) {}

	uint64_t get(uint64_t const i) const
	{
		return acc->get(i).first;
	}

	uint64_t operator[](uint64_t const i) const
	{
		return get(i);
	}

	uint64_t size() const
	{
		return acc->size();
	}
};

std::pair<uint64_t,uint64_t> getNextLarger(WordPairAccessor & W, uint64_t const p)
{
	WordPairAccessorFirstAdapter adp(W);
	libmaus2::util::ConstIterator<WordPairAccessorFirstAdapter,uint64_t> it(&adp);
	libmaus2::util::ConstIterator<WordPairAccessorFirstAdapter,uint64_t> ite = it;
	ite += adp.size();
	libmaus2::util::ConstIterator<WordPairAccessorFirstAdapter,uint64_t> itc =
		std::lower_bound(
			it,ite,p
		);

	if ( itc == ite )
		return W.get(0);
	else
		return W.get(itc-it);
}

std::pair<uint64_t,uint64_t> getNextLarger(std::string const & sortedisa, uint64_t const p)
{
	WordPairAccessor W(sortedisa);
	return getNextLarger(W,p);
}

struct SparseRank
{
	::libmaus2::huffman::IndexDecoderDataArray idda;
	libmaus2::aio::InputStreamInstance sparserankfd;
	int64_t const minsym;
	int64_t const maxsym;
	std::vector<uint64_t> const D;
	uint64_t const n;

	SparseRank(std::string const & bwt, std::string const & hist, int64_t const rminsym, int64_t const rmaxsym)
	: idda(std::vector<std::string>(1,bwt)), sparserankfd(bwt + ".sparserank"), minsym(rminsym), maxsym(rmaxsym), D(loadHMap(hist)),
	  n(std::accumulate(D.begin(),D.end(),0ull))
	{

	}


	static int64_t getMaxSym(std::string const & hist)
	{
		libmaus2::aio::InputStream::unique_ptr_type Phist(libmaus2::aio::InputStreamFactoryContainer::constructUnique(hist));
		std::map<int64_t,uint64_t> const hmap = ::libmaus2::util::NumberMapSerialisation::deserialiseMap<libmaus2::aio::InputStream,int64_t,uint64_t>(*Phist);

		if ( hmap.size() )
			return hmap.rbegin()->first;
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "LFStep::getMaxSym: no symbols" << std::endl;
			lme.finish();
			throw lme;

		}
	}

	static int64_t getMinSym(std::string const & hist)
	{
		libmaus2::aio::InputStream::unique_ptr_type Phist(libmaus2::aio::InputStreamFactoryContainer::constructUnique(hist));
		std::map<int64_t,uint64_t> const hmap = ::libmaus2::util::NumberMapSerialisation::deserialiseMap<libmaus2::aio::InputStream,int64_t,uint64_t>(*Phist);

		if ( hmap.size() )
			return hmap.begin()->first;
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "LFStep::getMinSym: no symbols" << std::endl;
			lme.finish();
			throw lme;

		}
	}

	uint64_t rankm(int64_t const sym, uint64_t const i)
	{
		libmaus2::huffman::FileBlockOffset const FBO = idda.findVBlock(i);
		assert ( FBO.fileptr == 0 );
		assert ( sym >= minsym );
		assert ( sym <= maxsym );
		sparserankfd.clear();
		sparserankfd.seekg(
			FBO.blockptr * (maxsym-minsym+1) * sizeof(uint64_t)
		);
		std::vector<uint64_t> H(maxsym-minsym+1,0);

		for ( uint64_t i = 0; i < H.size(); ++i )
			H[i] = libmaus2::util::NumberSerialisation::deserialiseNumber(sparserankfd);

		uint64_t offset = FBO.offset;
		::libmaus2::huffman::IndexDecoderData const & idata = idda.data[FBO.fileptr];
		::libmaus2::huffman::IndexEntry const data = idata.readEntry(FBO.blockptr);
		uint64_t const rlstart = data.vcnt;
		libmaus2::huffman::RLDecoder rldec(idda,rlstart);
		assert ( offset == i - rlstart );

		while ( offset )
		{
			std::pair<int64_t,uint64_t> const P = rldec.decodeRun();
			assert ( P.first != -1 );
			uint64_t touse = std::min(offset,P.second);
			H [ P.first - minsym ] += touse;
			offset -= touse;
		}

		return H[sym-minsym];
	}

	std::pair<int64_t,uint64_t> inverseSelect(uint64_t const i)
	{
		libmaus2::huffman::FileBlockOffset const FBO = idda.findVBlock(i);
		assert ( FBO.fileptr == 0 );
		sparserankfd.clear();
		sparserankfd.seekg(FBO.blockptr * (maxsym-minsym+1) * sizeof(uint64_t));
		std::vector<uint64_t> H(maxsym-minsym+1,0);

		for ( uint64_t i = 0; i < H.size(); ++i )
			H[i] = libmaus2::util::NumberSerialisation::deserialiseNumber(sparserankfd);

		uint64_t offset = FBO.offset;
		::libmaus2::huffman::IndexDecoderData const & idata = idda.data[FBO.fileptr];
		::libmaus2::huffman::IndexEntry const data = idata.readEntry(FBO.blockptr);
		uint64_t const rlstart = data.vcnt;
		libmaus2::huffman::RLDecoder rldec(idda,rlstart);
		assert ( offset == i - rlstart );

		while ( true )
		{
			std::pair<int64_t,uint64_t> const P = rldec.decodeRun();
			assert ( P.first != -1 );

			if ( offset >= P.second )
			{
				H [ P.first - minsym ] += P.second;
				offset -= P.second;
			}
			else
			{
				H [ P.first - minsym ] += offset;
				return std::pair<int64_t,uint64_t>(P.first,H [ P.first - minsym ]);
			}
		}
	}

	std::pair<int64_t,uint64_t> lf(uint64_t const i)
	{
		std::pair<int64_t,uint64_t> const P = inverseSelect(i);
		return std::pair<int64_t,uint64_t>(P.first,D [ P.first ] + P.second);
	}

	std::string decode(std::string const & sortedisa, uint64_t const low, uint64_t const len)
	{
		std::pair<uint64_t,uint64_t> const PP = getNextLarger(sortedisa,low+len);
		// std::cerr << "low+len=" << low+len << " PP.first=" << PP.first << " PP.second=" << PP.second << std::endl;
		uint64_t p = PP.first % n;
		uint64_t r = PP.second;

		while ( p != ((low+len)%n) )
		{
			// std::cerr << "p=" << p << " r=" << r << std::endl;
			r = lf(r).second;
			--p;
		}

		std::vector<char> V(len,0);
		for ( uint64_t i = 0; i < len; ++i )
		{
			std::pair<int64_t,uint64_t> const P = lf(r);

			// std::cerr << P.first << std::endl;

			if ( P.first )
				V[len-i-1] = libmaus2::fastx::remapChar(P.first-1);
			else
				V[len-i-1] = '\n';

			r = P.second;
		}

		return std::string(V.begin(),V.end());
	}
};

uint64_t computeNF(std::vector<uint64_t> const & H, std::vector<uint64_t> const & NZH, uint64_t const bm, std::vector<uint64_t> * HID = 0)
{
	uint64_t nf = 0;
	uint64_t s = 0;
	for ( uint64_t i = 0; i < NZH.size(); ++i )
	{
		if ( s && s >= bm )
		{
			s = 0;
			nf += 1;
		}

		uint64_t const sym = NZH[i];
		uint64_t const occ = H[sym+1]-H[sym];

		s += occ;

		// std::cerr << "sym=" << sym << " occ=" << occ << " s=" << s << std::endl;

		if ( HID )
			HID->push_back(nf);
	}

	if ( s )
	{
		s = 0;
		nf += 1;
	}

	return nf;
}

uint64_t computeSplit(uint64_t const tnumfiles, std::vector<uint64_t> const & H, std::vector<uint64_t> const & NZH, uint64_t const n)
{
	// look for smallest split value s.t. number of files is <= tnumfiles
	uint64_t bl = 0, bh = n;
	while ( bh - bl > 2 )
	{
		uint64_t const bm = (bh + bl)>>1;
		uint64_t const nf = computeNF(H,NZH,bm);

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
		if ( computeNF(H,NZH,i) <= tnumfiles )
		{
			split = i;
			break;
		}
	// std::cerr << "bl=" << bl << ",bh=" << bh << " split=" << split << " numfiles=" << computeNF(H,NZH,split) << std::endl;

	return split;
}

struct PreIsaInput
{
	typedef PreIsaInput this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus2::aio::ConcatInputStream CIS;
	libmaus2::aio::SynchronousGenericInput<uint64_t> SGI;

	PreIsaInput(std::vector<std::string> const & fn, uint64_t const roffset = 0)
	: CIS(fn), SGI(CIS,64*1024)
	{
		if ( roffset )
		{
			CIS.clear();
			CIS.seekg(roffset * 2 * sizeof(uint64_t));
		}
	}

	bool getNext(std::pair<uint64_t,uint64_t> & P)
	{
		bool const ok = SGI.getNext(P.first);

		if ( ok )
		{
			bool const ok2 = SGI.getNext(P.second);
			assert ( ok2 );
		}

		return ok;
	}
};

struct PreIsaAccessor
{
	typedef PreIsaAccessor this_type;

	mutable libmaus2::aio::ConcatInputStream CIS;
	// number of samples in file
	uint64_t n;

	typedef libmaus2::util::ConstIterator< this_type,std::pair<uint64_t,uint64_t> > iterator;

	PreIsaAccessor(std::vector<std::string> const & fn) : CIS(fn)
	{
		CIS.clear();
		CIS.seekg(0,std::ios::end);
		assert ( CIS.tellg() % (2*sizeof(uint64_t)) == 0 );
		n = CIS.tellg() / (2*sizeof(uint64_t));
		CIS.clear();
		CIS.seekg(0,std::ios::beg);
		// std::cerr << "constructed " << this << std::endl;
	}
	~PreIsaAccessor()
	{
		//std::cerr << "[V] destructing " << this << std::endl;
	}

	std::pair<uint64_t,uint64_t> get(uint64_t const i) const
	{
		//std::cerr << this << " " << "get(" << i << ")" << " n=" << n << std::endl;
		CIS.clear();
		CIS.seekg(i*2*sizeof(uint64_t));
		// std::cerr << CIS.tellg() << " " << i*2*sizeof(uint64_t) << std::endl;
		assert ( static_cast<std::streampos>(CIS.tellg()) == static_cast<std::streampos>(i*2*sizeof(uint64_t)) );
		uint64_t p0, p1;
		CIS.read(reinterpret_cast<char *>(&p0), sizeof(uint64_t));
		CIS.read(reinterpret_cast<char *>(&p1), sizeof(uint64_t));
		return std::pair<uint64_t,uint64_t>(p0,p1);
	}

	std::pair<uint64_t,uint64_t> operator[](uint64_t const i) const
	{
		return get(i);
	}

	iterator begin()
	{
		iterator it(this,0);
		return it;
	}

	iterator end()
	{
		iterator it(this,n);
		return it;
	}

	struct PFComp
	{
		bool operator()(std::pair<uint64_t,uint64_t> const & A, std::pair<uint64_t,uint64_t> const & B) const
		{
			return A.first < B.first;
		}
	};


	static std::vector< std::pair<uint64_t,uint64_t > > getOffsets(std::vector<std::string> const & fn, std::vector< uint64_t > const & splitblocks)
	{
		this_type acc(fn);

		assert ( splitblocks.size() );
		uint64_t const num = splitblocks.size()-1;

		std::vector< std::pair<uint64_t,uint64_t > >V(num+1);
		V[num] = std::pair<uint64_t,uint64_t >(splitblocks.back(), acc.n);

		for ( uint64_t i = 0; i < num; ++i )
		{
			uint64_t const rr = splitblocks[i];
			iterator it = std::lower_bound(acc.begin(),acc.end(),std::pair<uint64_t,uint64_t>(rr,0),PFComp());
			V[i] = std::pair<uint64_t,uint64_t > (rr,it-acc.begin());
		}

		return V;
	}
};

struct IsaTriple
{
	uint64_t s;
	uint64_t r;
	uint64_t p;

	IsaTriple() {}
	IsaTriple(uint64_t const & rs, uint64_t const & rr, uint64_t const & rp)
	: s(rs), r(rr), p(rp) {}

	bool operator<(IsaTriple const & rhs) const
	{
		if ( s != rhs.s )
			return s < rhs.s;
		else if ( r != rhs.r )
			return r < rhs.r;
		else
			return p < rhs.p;
	}

	bool operator==(IsaTriple const & rhs) const
	{
		return
			s == rhs.s
			&&
			r == rhs.r
			&&
			p == rhs.p;
	}

	bool operator!=(IsaTriple const & rhs) const
	{
		return !(*this == rhs);
	}
};

struct SaPair
{
	uint64_t r;
	uint64_t p;

	SaPair() {}
	SaPair(uint64_t const & rr, uint64_t const & rp)
	: r(rr), p(rp) {}

	bool operator<(SaPair const & rhs) const
	{
		if ( r != rhs.r )
			return r < rhs.r;
		else
			return p < rhs.p;
	}

	bool operator==(SaPair const & rhs) const
	{
		return
			r == rhs.r
			&&
			p == rhs.p;
	}

	bool operator!=(SaPair const & rhs) const
	{
		return !(*this == rhs);
	}
};

struct IsaPair
{
	uint64_t p;
	uint64_t r;

	IsaPair() {}
	IsaPair(uint64_t const & rp, uint64_t const & rr)
	: p(rp), r(rr) {}

	bool operator<(IsaPair const & rhs) const
	{
		if ( p != rhs.p )
			return p < rhs.p;
		else
			return r < rhs.r;
	}

	bool operator==(IsaPair const & rhs) const
	{
		return
			p == rhs.p
			&&
			r == rhs.r;
	}

	bool operator!=(IsaPair const & rhs) const
	{
		return !(*this == rhs);
	}
};

static std::vector< uint64_t > getSplitBlocks(uint64_t const tnum, uint64_t const bwtn)
{
	uint64_t const r = (bwtn + tnum-1)/tnum;
	uint64_t const num = (bwtn + r - 1)/r;
	std::vector<uint64_t> V(num+1);
	for ( uint64_t i = 0; i < num; ++i )
		V[i] = r*i;
	return V;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( arginfo.getNumRestArgs() < 4 )
		{
			std::cerr << "usage: " << argv[0] << " <bwt> <hist> <preisa> <preisasamplingrate>" << std::endl;
			return EXIT_FAILURE;
		}

		std::string bwt = arginfo.getUnparsedRestArg(0);
		std::string const origbwt = bwt;

		bool const copyinputtomemory = arginfo.getValue<uint64_t>("copyinputtomemory",false);

		if ( copyinputtomemory )
		{
			std::string const nbwt = "mem://copied_" + bwt;
			libmaus2::aio::InputStreamInstance ininst(bwt);
			libmaus2::aio::OutputStreamInstance outinst(nbwt);
			uint64_t const fs = libmaus2::util::GetFileSize::getFileSize(ininst);
			libmaus2::util::GetFileSize::copy(ininst,outinst,fs);
			std::cerr << "[V] copied " << bwt << " to " << nbwt << std::endl;
			bwt = nbwt;
		}


		std::string const histfn = arginfo.getUnparsedRestArg(1);
		std::string const preisa = arginfo.getUnparsedRestArg(2);
		uint64_t const preisasamplingrate = libmaus2::util::ArgInfo::parseValueUnsignedNumeric<uint64_t>("preisasamplingrate",arginfo.getUnparsedRestArg(3));
		uint64_t const sasamplingrate = arginfo.getValueUnsignedNumeric<uint64_t>("sasamplingrate",32);
		uint64_t const isasamplingrate = arginfo.getValueUnsignedNumeric<uint64_t>("isasamplingrate",32);
		std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpprefix",arginfo.getDefaultTmpFileName());

		libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type ref_isa;
		if ( arginfo.hasArg("ref_isa") )
		{
			libmaus2::aio::InputStreamInstance I(arginfo.getUnparsedValue("ref_isa",std::string()));
			uint64_t ref_isa_samplingrate;
			libmaus2::serialize::Serialize<uint64_t>::deserialize(I,&ref_isa_samplingrate);
			assert ( ref_isa_samplingrate == isasamplingrate );
			libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type Tptr(new libmaus2::autoarray::AutoArray<uint64_t>(I));
			ref_isa = UNIQUE_PTR_MOVE(Tptr);
		}
		libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type ref_sa;
		if ( arginfo.hasArg("ref_sa") )
		{
			libmaus2::aio::InputStreamInstance I(arginfo.getUnparsedValue("ref_sa",std::string()));
			uint64_t ref_sa_samplingrate;
			libmaus2::serialize::Serialize<uint64_t>::deserialize(I,&ref_sa_samplingrate);
			assert ( ref_sa_samplingrate == sasamplingrate );
			libmaus2::autoarray::AutoArray<uint64_t>::unique_ptr_type Tptr(new libmaus2::autoarray::AutoArray<uint64_t>(I));
			ref_sa = UNIQUE_PTR_MOVE(Tptr);
		}

		#if defined(_OPENMP)
		unsigned int const maxthreads = omp_get_max_threads();
		#else
		unsigned int const maxthreads = 1;
		#endif

		// get number of threads
		unsigned int const numthreads = arginfo.getValue<unsigned int>("threads", maxthreads);

		#if defined(_OPENMP)
		// set number of threads
		omp_set_num_threads(numthreads);
		#else
		if ( numthreads != 1 )
			std::cerr << "[V] warning, no OpenMP support, ignoring number of threads " << numthreads << std::endl;
		#endif

		// temp file for sa
		std::string const satmpfn = tmpfilenamebase + "_sa_tmp";
		libmaus2::util::TempFileRemovalContainer::addTempFile(satmpfn);
		// temp file for isa
		std::string const isatmpfn = tmpfilenamebase + "_isa_tmp";
		libmaus2::util::TempFileRemovalContainer::addTempFile(isatmpfn);

		// check sa sampling rate
		if ( libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(sasamplingrate) != 1 )
		{
			std::cerr << "sasamplingrate is not a power of 2" << std::endl;
			return EXIT_FAILURE;
		}
		// check isa sampling rate
		if ( libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(isasamplingrate) != 1 )
		{
			std::cerr << "isasamplingrate is not a power of 2" << std::endl;
			return EXIT_FAILURE;
		}

		// bit masks for sa and isa sampling rates
		uint64_t const sasamplingmask = sasamplingrate-1;
		uint64_t const isasamplingmask = isasamplingrate-1;

		// sorting output buffers for sa and isa
		libmaus2::sorting::SortingBufferedOutputFile<SaPair> sasorter(satmpfn,64*1024);
		libmaus2::sorting::SortingBufferedOutputFile<IsaPair> isasorter(isatmpfn,64*1024);

		// input file vector
		std::vector<std::string> Vbwtin(1,bwt);
		// length of input file in symbols
		uint64_t const n = libmaus2::huffman::RLDecoder::getLength(Vbwtin);

		std::vector<uint64_t> const H = loadHMap(histfn);

		uint64_t numnonzeroocc = 0;
		// sequence of non zero occ symbols
		std::vector<uint64_t> NZH;
		// map symbol to rank in sequence of non zero occ symbols
		std::vector<int64_t> NZHI;
		for ( uint64_t i = 1; i < H.size(); ++i )
			// if non zero occ count for symbol i-1
			if ( H[i] != H[i-1] )
			{
				numnonzeroocc++;
				// non zero occ symbol vector
				NZH.push_back(i-1);
				// rank mapping of non zero occ symbol
				NZHI.push_back(NZHI.size());
			}
			else
			{
				NZHI.push_back(-1);
			}

		// number of symbols with non zero occurence count
		std::cerr << "[D] numnonzeroocc=" << numnonzeroocc << std::endl;

		#if 0
		for ( uint64_t i = 0; i < NZHI.size(); ++i )
			std::cerr << i << " " << NZHI[i] << std::endl;
		#endif

		// target number of files
		uint64_t const tnumfiles = 256;
		// compute split
		uint64_t const split = computeSplit(tnumfiles, H, NZH, n);
		std::cerr << "[D] split=" << split << std::endl;

		// vector symbol to file for symbols in NZH
		std::vector<uint64_t> HID;
		computeNF(H, NZH, split, &HID);

		// map symbol -> file
		std::vector<int64_t> filemap;
		if ( NZH.size() )
			filemap.resize(NZH.back()+1);
		// initialise all to -1
		std::fill(filemap.begin(),filemap.end(),-1);

		// number of files used
		int64_t const numfiles = HID.size() ? static_cast<int64_t>(HID.back()+1) : 0;
		// bit vector marking files used by more than 1 symbol
		std::vector<bool> filemulti(numfiles);
		std::fill(filemulti.begin(),filemulti.end(),false);

		uint64_t hid_low = 0;
		while ( hid_low < HID.size() )
		{
			uint64_t hid_high = hid_low+1;
			while ( hid_high < HID.size() && HID[hid_high] == HID[hid_low] )
				++hid_high;

			if ( hid_high-hid_low > 1 )
				for ( uint64_t i = hid_low; i < hid_high; ++i )
					filemulti.at(HID[i]) = true;

			hid_low = hid_high;
		}

		for ( uint64_t i = 0; i < HID.size(); ++i )
		{
			int64_t const sym = NZH[i];
			uint64_t const mapping = HID[i];

			filemap .at ( sym ) = mapping;
			std::cerr << sym << "\t" << mapping << "\t" << filemulti[filemap[sym]] << "\t" << filemap[sym] << "\t" << H[sym] << std::endl;
		}

		std::cerr << "[V] numfiles=" << numfiles << std::endl;

		std::vector<std::string> isain(1,preisa);
		bool deletein = false;

		int64_t const maxsym = static_cast<int64_t>(H.size())-1;

		std::vector< std::vector<uint64_t> > Vblocksymhist;
		// block rank start vector
		std::vector< uint64_t > Vsplitblocks;

		if ( 1 /* maxsym > 256 */ )
		{
			std::string const rlhisttmp = tmpfilenamebase + "_rlhist_tmp";
			std::string const blockhisttmp = tmpfilenamebase + "_rlhist_final";
			libmaus2::util::TempFileRemovalContainer::addTempFile(blockhisttmp);
			libmaus2::huffman::RLDecoder::getFragmentSymHistograms(bwt,blockhisttmp,rlhisttmp,0,maxsym,numthreads);

			std::vector<uint64_t> rlblockhist(maxsym+1,0);
			std::vector<uint64_t> rlblockhistacc(maxsym+1,0);
			libmaus2::aio::InputStreamInstance::unique_ptr_type rlhistISI(new libmaus2::aio::InputStreamInstance(blockhisttmp));

			// decode block of 0 values
			for ( uint64_t i = 0; i < rlblockhist.size(); ++i )
				rlblockhist[i] = libmaus2::util::NumberSerialisation::deserialiseNumber(*rlhistISI) - rlblockhistacc[i];

			Vsplitblocks.push_back(std::accumulate(rlblockhistacc.begin(),rlblockhistacc.end(),0ull));
			assert ( Vsplitblocks.at(0) == 0ull );

			while ( rlhistISI->peek() != std::istream::traits_type::eof() )
			{
				for ( uint64_t i = 0; i < rlblockhist.size(); ++i )
					rlblockhist[i] = libmaus2::util::NumberSerialisation::deserialiseNumber(*rlhistISI) - rlblockhistacc[i];

				std::vector<uint64_t> bhist(maxsym+1,0);
				for ( uint64_t i = 0; i < rlblockhistacc.size(); ++i )
					bhist[i] = H[i] + rlblockhistacc[i];
				Vblocksymhist.push_back(bhist);

				for ( uint64_t i = 0; i < rlblockhist.size(); ++i )
					rlblockhistacc[i] += rlblockhist[i];

				Vsplitblocks.push_back(std::accumulate(rlblockhistacc.begin(),rlblockhistacc.end(),0ull));
			}

			#if 0
			std::vector<uint64_t> bhist(maxsym+1,0);
			for ( uint64_t i = 0; i < rlblockhistacc.size(); ++i )
				bhist[i] = H[i] + rlblockhistacc[i];
			Vblocksymhist.push_back(bhist);
			#endif

			assert ( Vsplitblocks.back() == n );
		}
		else
		{
			std::string const rlhisttmp = tmpfilenamebase + "_rlhist_tmp";
			std::string const blockhisttmp = tmpfilenamebase + "_rlhist_final";
			libmaus2::util::TempFileRemovalContainer::addTempFile(blockhisttmp);
			libmaus2::huffman::RLDecoder::getBlockSymHistograms(bwt,blockhisttmp,rlhisttmp,0,maxsym);

			libmaus2::aio::InputStreamInstance::unique_ptr_type rlhistISI(new libmaus2::aio::InputStreamInstance(blockhisttmp));
			// rlhistISI->seekg(0,std::ios::end);
			// std::cerr << rlhistISI->tellg() << std::endl;
			// rlhistISI->clear();
			// rlhistISI->seekg(0,std::ios::beg);

			uint64_t const tpacks = numthreads;
			uint64_t const packsize = (n + tpacks-1)/tpacks;

			std::vector<uint64_t> rlblockhistacc(maxsym+1,0);
			std::vector<uint64_t> rlblockhistprev(maxsym+1,0);
			std::vector<uint64_t> rlblockhist(maxsym+1,0);
			uint64_t ac = 0;
			uint64_t bc = 0;
			uint64_t maxbs = 0;
			std::vector<uint64_t> Vbs;
			while ( rlhistISI->peek() != std::istream::traits_type::eof() )
			{
				// compute block size
				uint64_t bs = 0;
				for ( uint64_t i = 0; i < rlblockhist.size(); ++i )
				{
					rlblockhist[i] = libmaus2::util::NumberSerialisation::deserialiseNumber((*rlhistISI));
					bs += (rlblockhist[i]-rlblockhistprev[i]);
				}
				maxbs = std::max(bs,maxbs);

				if ( ac > packsize )
				{
					assert ( std::accumulate(rlblockhistacc.begin(),rlblockhistacc.end(),0ull) == ac );

					Vbs.push_back(ac);
					Vblocksymhist.push_back(rlblockhistacc);
					// std::cerr << "bc=" << bc << std::endl;

					ac = 0;
					bc = 0;
					std::fill(rlblockhistacc.begin(),rlblockhistacc.end(),0ull);
				}

				for ( uint64_t i = 0; i < rlblockhist.size(); ++i )
				{
					rlblockhistacc[i] += rlblockhist[i]-rlblockhistprev[i];
					rlblockhistprev[i] = rlblockhist[i];
				}
				ac += bs;
				bc += 1;
			}

			assert ( std::accumulate(rlblockhistacc.begin(),rlblockhistacc.end(),0ull) == ac );

			if ( ac )
			{
				Vbs.push_back(ac);
				Vblocksymhist.push_back(rlblockhistacc);
				std::fill(rlblockhistacc.begin(),rlblockhistacc.end(),0ull);
			}

			Vblocksymhist.push_back(rlblockhistacc);

			#if 0
			std::cerr << "rlblockhist.size()=" << rlblockhist.size() << " " << Vbs.size() << " " << packsize << " maxbs=" << maxbs << std::endl;
			for ( uint64_t i = 0; i < Vblocksymhist.size(); ++i )
			{
				std::cerr << std::accumulate(Vblocksymhist[i].begin(),Vblocksymhist[i].end(),0ull) << std::endl;
			}

			exit(0);
			#endif

			Vsplitblocks = std::vector< uint64_t >(Vbs.size()+1);
			uint64_t as = 0;
			for ( uint64_t i = 0; i < Vsplitblocks.size(); ++i )
			{
				Vsplitblocks[i] = as;
				as += Vbs[i];
			}
			Vsplitblocks[Vbs.size()] = as;

			// iterate over symbols
			for ( uint64_t i = 0; i < rlblockhistacc.size(); ++i )
			{
				uint64_t c = H[i];

				// iterate over blocks
				for ( uint64_t j = 0; j < Vblocksymhist.size(); ++j )
				{
					uint64_t const t = Vblocksymhist[j][i];
					Vblocksymhist[j][i] = c;
					c += t;
				}
			}

			#if 0
			assert (
				std::accumulate(
					Vblocksymhist.back().begin(),
					Vblocksymhist.back().end(),
					0ull
				) == n
			);
			#endif
		}

		libmaus2::parallel::PosixSpinLock salock;
		libmaus2::parallel::PosixSpinLock isalock;

		// iterate pre isa sampling rate times (pre isa samples are at most at distance preisasamplingrate from each other)
		for ( uint64_t run = 0; run < preisasamplingrate; ++run )
		{
			std::cerr << "[V] entering run " << run << std::endl;

			{
				PreIsaInput PII(isain,0);
				std::pair<uint64_t,uint64_t> P, Pprev;
				bool haveprev = false;
				while ( PII.getNext(P) )
				{
					if ( haveprev )
					{
						bool ok = ( P.first >= Pprev.first );
						ok = ok && ( P.first != Pprev.first || P.second == Pprev.second );
						if ( ! ok )
						{
							std::cerr << "failed order " << Pprev.first << " " << P.first << std::endl;
						}
						assert ( ok );
					}
					Pprev = P;
					haveprev = true;
				}
			}

			std::vector< std::pair<uint64_t,uint64_t > > const piasplit = PreIsaAccessor::getOffsets(isain,Vsplitblocks);
			assert ( piasplit.size() == Vsplitblocks.size() );

			#if 0
			for ( uint64_t i = 0; i < piasplit.size(); ++i )
			{
				PreIsaAccessor PIA(isain);
				std::cerr << "[PIA] " << i << " " << piasplit[i].first << " " << piasplit[i].second << " " << PIA[piasplit[i].second].first << std::endl;
				assert ( i == 0 || PIA[piasplit[i].second-1].first < piasplit[i].first );
			}
			#endif

			std::vector<std::string> goutfilenames((piasplit.size()-1)*numfiles);
			libmaus2::parallel::PosixSpinLock goutfilenameslock;

			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( uint64_t t = 1; t < piasplit.size(); ++t )
			{
				uint64_t const rfrom = piasplit[t-1].first;
				uint64_t const rto = piasplit[t].first;
				uint64_t const preisafrom = piasplit[t-1].second;

				salock.lock();
				std::cerr << "[V] t=" << t << " rfrom=" << rfrom << " rto=" << rto << " preisafrom=" << preisafrom << std::endl;
				salock.unlock();

				std::ostringstream threadtmpprefixstr;
				threadtmpprefixstr << tmpfilenamebase << "_" << run << "_" << (t-1);
				std::string const threadtmpprefix = threadtmpprefixstr.str();

				// output file pointers
				libmaus2::autoarray::AutoArray<libmaus2::aio::OutputStreamInstance::unique_ptr_type > Aout(numfiles);
				// output file buffers
				libmaus2::autoarray::AutoArray<libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type > Sout(numfiles);
				// sorting output buffers (for files storing data for more than one symbol)
				libmaus2::autoarray::AutoArray<libmaus2::sorting::SortingBufferedOutputFile<IsaTriple>::unique_ptr_type > Tout(numfiles);
				// output file names
				std::vector<std::string> outfilenames(numfiles);
				// multi sort output file names
				std::vector<std::string> multisortfilenames(numfiles);

				//std::cerr << "[V+] " << numfiles << std::endl;

				// iterate over output files for this run
				for ( int64_t i = 0; i < numfiles; ++i )
				{
					std::ostringstream fnostr;
					fnostr << threadtmpprefix << "_" << std::setw(6) << std::setfill('0') << i << std::setw(0);
					std::string const fn = fnostr.str();
					outfilenames[i] = fn;

					libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
					libmaus2::aio::OutputStreamInstance::unique_ptr_type Tptr(new libmaus2::aio::OutputStreamInstance(fn));
					Aout[i] = UNIQUE_PTR_MOVE(Tptr);

					libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type Sptr(new libmaus2::aio::SynchronousGenericOutput<uint64_t>(*(Aout[i]),64*1024));
					Sout[i] = UNIQUE_PTR_MOVE(Sptr);

					multisortfilenames[i] = fn + ".multisort";
					if ( filemulti[i] )
					{
						libmaus2::util::TempFileRemovalContainer::addTempFile(multisortfilenames[i]);
						libmaus2::sorting::SortingBufferedOutputFile<IsaTriple>::unique_ptr_type Tptr(
							new libmaus2::sorting::SortingBufferedOutputFile<IsaTriple>(multisortfilenames[i],64*1024)
						);
						Tout[i] = UNIQUE_PTR_MOVE(Tptr);
					}
				}

				goutfilenameslock.lock();

				for ( int64_t i = 0; i < numfiles; ++i )
				{
					goutfilenames . at( i * (piasplit.size()-1) + (t-1) ) = outfilenames[i];
				}

				goutfilenameslock.unlock();

				// open decoder for RL
				libmaus2::huffman::RLDecoder rldec(std::vector<std::string>(1,bwt),rfrom);

				std::vector<uint64_t> HH = Vblocksymhist[t-1];

				PreIsaInput PII(isain,preisafrom);

				std::pair<uint64_t,uint64_t> P;
				#if 0
				uint64_t lc = rfrom;
				#endif
				uint64_t c = rfrom;

				while ( PII.getNext(P) && (P.first < rto) )
				{
					bool const ok = P.first >= rfrom;

					if ( ! ok )
					{
						std::cerr << "P.first=" << P.first << " rfrom=" << rfrom << std::endl;
					}

					assert ( P.first >= rfrom );

					if ( (P.first & sasamplingmask) == 0 )
					{
						salock.lock();
						sasorter.put(SaPair(P.first,P.second));
						salock.unlock();

						if ( ref_sa )
						{
							bool const ok = (*ref_sa)[P.first / sasamplingrate] == P.second;
							if ( ! ok )
							{
								std::cerr << "P.first=" << P.first << std::endl;
								std::cerr << "P.first/sasamplingrate=" << P.first/sasamplingrate << std::endl;
								std::cerr << "P.second=" << P.second << std::endl;
								std::cerr << "expected=" << (*ref_sa)[P.first / sasamplingrate] << std::endl;
							}
							assert ( ok );
							//std::cerr << "ok 1" << std::endl;
						}
					}
					if ( (P.second & isasamplingmask) == 0 )
					{
						isalock.lock();
						isasorter.put(IsaPair(P.second,P.first));
						isalock.unlock();

						if ( ref_isa )
						{
							assert ( (*ref_isa)[P.second / isasamplingrate] == P.first );
							//std::cerr << "ok 2" << std::endl;
						}
					}

					#if 0
					if ( (c >> 26) != (lc >> 26) )
					{
						std::cerr << "run=" << (run+1) << "/" << preisasamplingrate << " frac=" << static_cast<double>(c) / n << std::endl;
						lc = c;
					}
					#endif

					while ( c < P.first )
					{
						int64_t const sym = rldec.decode();
						assert ( sym >= 0 );
						assert ( sym < static_cast<int64_t>(HH.size()) );
						HH[sym] += 1;
						c += 1;
					}

					int64_t const sym = rldec.decode();
					assert ( sym >= 0 );
					assert ( sym < static_cast<int64_t>(HH.size()) );

					uint64_t const rout = HH[sym];
					uint64_t const pout = (P.second + n - 1)%n;
					assert ( filemap[sym] >= 0 );
					uint64_t const fileid = filemap[sym];

					// std::cerr << "sym=" << sym << " fileid=" << fileid << std::endl;

					bool const multi = filemulti[fileid];

					if ( multi )
					{
						Tout[fileid]->put(IsaTriple(sym,rout,pout));
					}
					else
					{
						Sout[fileid]->put(rout);
						Sout[fileid]->put(pout);
					}
					HH[sym] += 1;
					c += 1;
				}

				// post sort files storing data for more than one symbol
				for ( int64_t fileid = 0; fileid < numfiles; ++fileid )
					if ( filemulti[fileid] )
					{
						libmaus2::sorting::SortingBufferedOutputFile<IsaTriple>::merger_ptr_type merger(Tout[fileid]->getMerger());
						IsaTriple T;
						IsaTriple Tprev;
						bool haveprev = false;

						while ( merger->getNext(T) )
						{
							if ( haveprev )
								assert ( T.r >= Tprev.r );

							Sout[fileid]->put(T.r);
							Sout[fileid]->put(T.p);

							Tprev = T;
							haveprev = true;
						}
					}

				for ( int64_t i = 0; i < numfiles; ++i )
				{
					Sout[i]->flush();
					Sout[i].reset();
					Aout[i]->flush();
					Aout[i].reset();
				}

				for ( int64_t i = 0; i < numfiles; ++i )
					libmaus2::aio::FileRemoval::removeFile(multisortfilenames[i]);

			}

			for ( int64_t fileid = 0; fileid < numfiles; ++fileid )
				if ( filemulti[fileid] )
				{
					struct QueueElement
					{
						uint64_t r;
						uint64_t p;
						uint64_t id;

						QueueElement(uint64_t const rr = 0, uint64_t const rp = 0, uint64_t const rid = 0)
						: r(rr), p(rp), id(rid)
						{

						}

						bool operator<(QueueElement const & O) const
						{
							if ( r != O.r )
								return r > O.r;
							else if ( p != O.p )
								return p > O.p;
							else
								return id > O.id;
						}
					};

					std::priority_queue<QueueElement> Q;
					libmaus2::autoarray::AutoArray<PreIsaInput::unique_ptr_type> Ain(piasplit.size()-1);
					for ( uint64_t t = 1; t < piasplit.size(); ++t )
					{
						uint64_t const id = t-1;

						PreIsaInput::unique_ptr_type Tptr(new PreIsaInput(std::vector<std::string>(1,goutfilenames . at( fileid * (piasplit.size()-1) + (id) ))));
						Ain.at(id) = UNIQUE_PTR_MOVE(Tptr);

						std::pair<uint64_t,uint64_t> P;
						if ( Ain[id]->getNext(P) )
							Q.push(QueueElement(P.first,P.second,id));
					}

					std::ostringstream tmpfnstr;
					tmpfnstr << tmpfilenamebase << "_" << run << "_merge_" << fileid;
					std::string const tmpfn = tmpfnstr.str();
					libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfn);
					libmaus2::aio::OutputStreamInstance::unique_ptr_type out(new libmaus2::aio::OutputStreamInstance(tmpfn));
					libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type sout(new libmaus2::aio::SynchronousGenericOutput<uint64_t>(*out,8*1024));
					std::pair<uint64_t,uint64_t> Pprev;
					bool haveprev = false;
					while ( Q.size() )
					{
						QueueElement const E = Q.top();
						Q.pop();

						if ( haveprev )
						{
							assert ( E.r >= Pprev.first );
							assert ( E.r != Pprev.first || E.p == Pprev.second );
						}

						sout->put(E.r);
						sout->put(E.p);

						std::pair<uint64_t,uint64_t> P;
						if ( Ain.at(E.id)->getNext(P) )
							Q.push(QueueElement(P.first,P.second,E.id));

						Pprev.first = E.r;
						Pprev.second = E.p;
						haveprev = true;
					}
					sout->flush();
					sout.reset();
					out.reset();

					{
					PreIsaInput PII(std::vector<std::string>(1,tmpfn));
					std::pair<uint64_t,uint64_t> P, Pprev;
					bool haveprev = false;
					while ( PII.getNext(P) )
					{
						if ( haveprev )
						{
							assert ( P.first >= Pprev.first );
							assert ( P.first > Pprev.first || P.second == Pprev.second );
						}
						haveprev = true;
						Pprev = P;
					}
					}

					for ( uint64_t t = 1; t < piasplit.size(); ++t )
					{
						uint64_t const id = t-1;
						std::string & fn = goutfilenames . at( fileid * (piasplit.size()-1) + (id) );
						if ( id )
						{
							libmaus2::aio::OutputStreamInstance tout(fn);
						}
						else
						{
							libmaus2::aio::FileRemoval::removeFile(fn);
							fn = tmpfn;
						}
					}
				}

			if ( deletein )
				for ( uint64_t i = 0; i < isain.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile(isain[i]);

			isain = goutfilenames;
			deletein = true;
		}

		libmaus2::aio::OutputStreamInstance saout(origbwt + ".sa");
		::libmaus2::serialize::Serialize<uint64_t>::serialize(saout,sasamplingrate);
		::libmaus2::serialize::Serialize<uint64_t>::serialize(saout,(n+sasamplingrate-1)/sasamplingrate);
		libmaus2::sorting::SortingBufferedOutputFile<SaPair>::merger_ptr_type samerger(sasorter.getMerger());
		SaPair SaT;
		uint64_t prevr = std::numeric_limits<uint64_t>::max();
		uint64_t esa = 0;
		while ( samerger->getNext(SaT) )
			if ( SaT.r != prevr )
			{
				if ( ref_sa )
					assert ( (*ref_sa)[SaT.r/sasamplingrate] == SaT.p );

				assert ( SaT.r % sasamplingrate == 0 );
				bool const ok = (SaT.r/sasamplingrate) == (esa++);
				assert ( ok );

				::libmaus2::serialize::Serialize<uint64_t>::serialize(saout,SaT.p);
				prevr = SaT.r;
			}

		libmaus2::aio::OutputStreamInstance isaout(origbwt + ".isa");
		::libmaus2::serialize::Serialize<uint64_t>::serialize(isaout,isasamplingrate);
		::libmaus2::serialize::Serialize<uint64_t>::serialize(isaout,(n+isasamplingrate-1)/isasamplingrate);
		libmaus2::sorting::SortingBufferedOutputFile<IsaPair>::merger_ptr_type isamerger(isasorter.getMerger());
		IsaPair IsaT;
		uint64_t prevp = std::numeric_limits<uint64_t>::max();
		uint64_t eisa = 0;
		while ( isamerger->getNext(IsaT) )
			if ( IsaT.p != prevp )
			{
				assert ( IsaT.p % isasamplingrate == 0 );

				if ( ref_isa )
					assert ( (*ref_isa)[IsaT.p/isasamplingrate] == IsaT.r );

				bool const ok = (IsaT.p/isasamplingrate) == (eisa++);
				assert ( ok );

				::libmaus2::serialize::Serialize<uint64_t>::serialize(isaout,IsaT.r);
				prevp = IsaT.p;
			}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
