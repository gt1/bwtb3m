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

	static std::vector<uint64_t> loadHMap(std::string const & hist)
	{
		libmaus2::aio::InputStream::unique_ptr_type Phist(libmaus2::aio::InputStreamFactoryContainer::constructUnique(hist));
		std::map<int64_t,uint64_t> const hmap = ::libmaus2::util::NumberMapSerialisation::deserialiseMap<libmaus2::aio::InputStream,int64_t,uint64_t>(*Phist);
		std::vector<uint64_t> D;

		if ( hmap.size() && hmap.begin()->first >= 0 )
		{
			D.resize(hmap.rbegin()->first+1);
			uint64_t s = 0;
			for ( std::map<int64_t,uint64_t>::const_iterator ita = hmap.begin(); ita != hmap.end(); ++ita )
			{
				D [ ita->first ] = s;
				s += ita->second;
				// std::cerr << ita->first << "\t" << ita->second << std::endl;
			}
			#if 0
			for ( uint64_t i = 0; i < D.size(); ++i )
				std::cerr << i << "\t" << D[i] << std::endl;
			#endif
		}

		return D;
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

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		std::string const bwt = arginfo.getUnparsedRestArg(0);
		std::string const histfn = arginfo.getUnparsedRestArg(1);
		std::string const sortedisa = arginfo.getUnparsedRestArg(2);
		std::string const compact = arginfo.getUnparsedRestArg(3);
		uint64_t low = libmaus2::util::ArgInfo::parseValueUnsignedNumeric<uint64_t>("low",arginfo.getUnparsedRestArg(4));
		uint64_t len = libmaus2::util::ArgInfo::parseValueUnsignedNumeric<uint64_t>("len",arginfo.getUnparsedRestArg(5));

		libmaus2::aio::InputStream::unique_ptr_type Psortedisa(libmaus2::aio::InputStreamFactoryContainer::constructUnique(sortedisa));

		std::vector<std::string> Vbwtin(1,bwt);

		uint64_t const n = libmaus2::huffman::RLDecoder::getLength(Vbwtin);

		if ( low+len > n )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "low+len=" << low+len << " > n=" << n << std::endl;
			lme.finish();
			throw lme;
		}

		int64_t const maxsym = SparseRank::getMaxSym(histfn);
		if ( maxsym < 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "maxsym < 0" << std::endl;
			lme.finish();
			throw lme;
		}
		int64_t const minsym = SparseRank::getMinSym(histfn);

		if ( ! libmaus2::aio::InputStreamFactoryContainer::tryOpen(bwt + ".sparserank") )
			libmaus2::huffman::RLDecoder::getBlockSymHistograms(
				bwt,
				bwt + ".sparserank",
				bwt + "_prerank",
				minsym,
				maxsym,
				&std::cerr
			);

		SparseRank SR(bwt,histfn,minsym,maxsym);

		uint64_t const bs = 2048;
		while ( len )
		{
			uint64_t const todecode = std::min(len,bs);
			std::cerr << "[V] decoding [" << low << "," << low+todecode << ")" << std::endl;
			std::string const block = SR.decode(sortedisa,low,todecode);
			std::cout << block;

			libmaus2::bitio::CompactDecoderWrapper wrap(compact);
			wrap.clear();
			wrap.seekg(low);
			for ( uint64_t i = 0; i < len; ++i )
			{
				char const c = wrap.get();
				if ( c == 0 )
				{
					assert ( block[i] == '\n' );
				}
				else
				{
					assert ( block[i] == libmaus2::fastx::remapChar(c-1) );
				}
			}

			len -= todecode;
			low += todecode;
		}

		#if 0
		libmaus2::huffman::RLDecoder rldec(std::vector<std::string>(1,bwt));
		std::vector<uint64_t> H(maxsym-minsym+1,0);
		for ( uint64_t i = 0; i < 128000; ++i )
		{
			int64_t const sym = rldec.decode();

			for ( int64_t j = minsym; j <= maxsym; ++j )
			{
				uint64_t const r = SR.rankm(j,i);
				std::cerr << i << " " << j << " " << r << " " << H[j-minsym] << std::endl;
				assert ( r == H[j-minsym] );
			}

			assert ( sym == SR.inverseSelect(i).first );
			assert ( H[sym-minsym] == SR.inverseSelect(i).second );

			H [ sym - minsym ] += 1;
		}
		#endif



		#if 0
		for ( uint64_t i = 0; i < 1024; ++i )
		{
			std::pair<uint64_t,uint64_t> const P = getNextLarger(sortedisa,i);
			std::cerr << i << "\t" << P.first << "," << P.second << std::endl;
		}

		if ( n >= 1024 )
			for ( uint64_t i = n-1024; i <= n; ++i )
			{
				std::pair<uint64_t,uint64_t> const P = getNextLarger(sortedisa,i);
				std::cerr << i << "\t" << P.first << "," << P.second << std::endl;
			}
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
