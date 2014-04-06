/**
    bwtb3m
    Copyright (C) 2009-2014 German Tischler
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

#include <libmaus/wavelet/ImpExternalWaveletGeneratorCompactHuffman.hpp>
#include <libmaus/wavelet/ImpExternalWaveletGeneratorCompactHuffmanParallel.hpp>

#include <iostream>
#include <iomanip>

// #define HUFGAP
#define HUFRL

#include <libmaus/aio/CheckedInputStream.hpp>
#include <libmaus/aio/CheckedOutputStream.hpp>
#include <libmaus/aio/CircularWrapper.hpp>
#include <libmaus/aio/FileFragment.hpp>
#include <libmaus/aio/ReorderConcatGenericInput.hpp>

#include <libmaus/autoarray/AutoArray.hpp>

#include <libmaus/bitio/BitStreamFileDecoder.hpp>
#include <libmaus/bitio/CompactDecoderBuffer.hpp>
#include <libmaus/bitio/BitVectorInput.hpp>
#include <libmaus/bitio/BitVectorOutput.hpp>

#include <libmaus/gamma/GammaGapEncoder.hpp>
#include <libmaus/gamma/GammaGapDecoder.hpp>
#include <libmaus/gamma/SparseGammaGapFileSet.hpp>
#include <libmaus/gamma/SparseGammaGapFileLevelSet.hpp>
#include <libmaus/gamma/SparseGammaGapMultiFileLevelSet.hpp>
#include <libmaus/gamma/SparseGammaGapEncoder.hpp>
#include <libmaus/gamma/SparseGammaGapDecoder.hpp>

#include <libmaus/huffman/GapDecoder.hpp>
#include <libmaus/huffman/GapEncoder.hpp>
#include <libmaus/huffman/HuffmanTree.hpp>

#include <libmaus/gamma/GammaGapEncoder.hpp>
#include <libmaus/gamma/GammaGapDecoder.hpp>
#include <libmaus/gamma/GammaRLEncoder.hpp>
#include <libmaus/gamma/GammaRLDecoder.hpp>

#include <libmaus/huffman/huffman.hpp>
#include <libmaus/huffman/HuffmanEncoderFile.hpp>
#include <libmaus/huffman/RLEncoder.hpp>
#include <libmaus/huffman/RLDecoder.hpp>

#include <libmaus/lf/DArray.hpp>
#include <libmaus/lf/LF.hpp>

#include <libmaus/math/numbits.hpp>

#include <libmaus/parallel/SynchronousCounter.hpp>
#include <libmaus/parallel/LockedBool.hpp>

#include <libmaus/sorting/PairFileSorting.hpp>

#include <libmaus/suffixsort/BwtMergeBlockSortResult.hpp>
#include <libmaus/suffixsort/BwtMergeTempFileNameSet.hpp>
#include <libmaus/suffixsort/BwtMergeTempFileNameSetVector.hpp>
#include <libmaus/suffixsort/BwtMergeZBlock.hpp>
#include <libmaus/suffixsort/BwtMergeZBlockRequest.hpp>
#include <libmaus/suffixsort/BwtMergeZBlockRequestVector.hpp>
#include <libmaus/suffixsort/ByteInputTypes.hpp>
#include <libmaus/suffixsort/CircularBwt.hpp>
#include <libmaus/suffixsort/CircularSuffixComparator.hpp>
#include <libmaus/suffixsort/CompactInputTypes.hpp>
#include <libmaus/suffixsort/divsufsort.hpp>
#include <libmaus/suffixsort/GapArrayByte.hpp>
#include <libmaus/suffixsort/GapMergePacket.hpp>
#include <libmaus/suffixsort/PacInputTypes.hpp>
#include <libmaus/suffixsort/PacTermInputTypes.hpp>
#include <libmaus/suffixsort/Utf8InputTypes.hpp>
#include <libmaus/suffixsort/Lz4InputTypes.hpp>

#include <libmaus/wavelet/ImpExternalWaveletGeneratorHuffman.hpp>
#include <libmaus/wavelet/ImpExternalWaveletGeneratorHuffmanParallel.hpp>
#include <libmaus/wavelet/ImpHuffmanWaveletTree.hpp>
#include <libmaus/wavelet/Utf8ToImpCompactHuffmanWaveletTree.hpp>
#include <libmaus/wavelet/Utf8ToImpHuffmanWaveletTree.hpp>

#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/FileTempFileContainer.hpp>
#include <libmaus/util/GetFileSize.hpp>
#include <libmaus/util/Histogram.hpp>
#include <libmaus/util/HistogramSet.hpp>
#include <libmaus/util/KMP.hpp>
#include <libmaus/util/MemUsage.hpp>
#include <libmaus/util/NumberMapSerialisation.hpp>
#include <libmaus/util/OctetString.hpp>
#include <libmaus/util/OutputFileNameTools.hpp>
#include <libmaus/util/StringSerialisation.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/util/Utf8String.hpp>
#include <libmaus/util/SimpleCountingHash.hpp>

#if defined(HUFRL)
typedef ::libmaus::huffman::RLDecoder rl_decoder;
#else
typedef ::libmaus::gamma::GammaRLDecoder rl_decoder;
#endif

libmaus::parallel::OMPLock gcerrlock;
libmaus::parallel::OMPLock glock;

struct WtTodo
{
	uint64_t low;
	uint64_t high;
	unsigned int depth;
	uint32_t node;
	
	WtTodo() : low(0), high(0), depth(0), node(0) {}
	WtTodo(
		uint64_t const rlow,
		uint64_t const rhigh,
		unsigned int const rdepth,
		uint32_t const rnode
	) : low(rlow), high(rhigh), depth(rdepth), node(rnode) {}
};

std::ostream & operator<<(std::ostream & out, WtTodo const & O)
{
	out << "WtTodo(low=" << O.low << ",high=" << O.high << ",depth=" << O.depth << ",node=" << O.node << ")";
	return out;
}

template<bool _utf8_input_type>
struct RlToHwtBase
{
	static bool utf8Wavelet()
	{
		return _utf8_input_type;
	}
	
	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type loadWaveletTree(std::string const & hwt)
	{
		libmaus::aio::CheckedInputStream CIS(hwt);
		libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type pICHWT(new libmaus::wavelet::ImpCompactHuffmanWaveletTree(CIS));
		return UNIQUE_PTR_MOVE(pICHWT);
	}
	
	struct RlDecoderInfoObject
	{
		rl_decoder * dec;
		uint64_t block;
		uint64_t blocksleft;
		
		uint64_t low;
		uint64_t end;
		uint64_t blocksize;
		uint64_t blockoffset;
		
		RlDecoderInfoObject() : dec(0), block(0), blocksleft(0), low(0), end(0), blocksize(0), blockoffset(0) {}
		RlDecoderInfoObject(
			rl_decoder * rdec,
			uint64_t const rblock,
			uint64_t const rblocksleft,
			uint64_t const rlow,
			uint64_t const rend,
			uint64_t const rblocksize,
			uint64_t const rblockoffset
		) : dec(rdec), block(rblock), blocksleft(rblocksleft), low(rlow), end(rend), blocksize(rblocksize), blockoffset(rblockoffset)
		{
		
		}
		
		bool hasNextBlock() const
		{
			return blocksleft > 1;
		}
		
		RlDecoderInfoObject nextBlock() const
		{
			assert ( hasNextBlock() );
			return RlDecoderInfoObject(dec,block+1,blocksleft-1,low+blocksize,end,blocksize,blockoffset);
		}
		
		uint64_t getLow() const
		{
			return low;	
		}
		
		uint64_t getHigh() const
		{
			return std::min(low+blocksize,end);
		}
		
		uint64_t getAbsoluteBlock() const
		{
			return block + blockoffset;
		}
	};
	
	template<typename _object_type>
	struct Todo
	{
		typedef _object_type object_type;
		std::stack<object_type> todo;
		libmaus::parallel::OMPLock lock;
		
		Todo()
		{
		
		}
		
		void push(object_type const & o)
		{
			libmaus::parallel::ScopeLock sl(lock);
			todo.push(o);
		}
		
		bool pop(object_type & o)
		{
			libmaus::parallel::ScopeLock sl(lock);
			if ( todo.size() )
			{
				o = todo.top();
				todo.pop();
				return true;
			}
			else
			{
				return false;
			}	
		}
	};

	/**
	 * specialised version for small alphabets
	 **/
	template<typename entity_type = uint8_t>
	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type rlToHwtTermSmallAlphabet(
		std::vector<std::string> const & bwt, 
		std::string const & huftreefilename,
		uint64_t const bwtterm,
		uint64_t const p0r
		)
	{
		// std::cerr << "(" << libmaus::util::Demangle::demangle<entity_type>() << ")";
	
		// load the huffman tree
		::libmaus::huffman::HuffmanTree::unique_ptr_type UH = loadCompactHuffmanTree(huftreefilename);
		::libmaus::huffman::HuffmanTree & H = *UH;
		// check depth is low
		assert ( H.maxDepth() <= 8*sizeof(entity_type) );
		// get encoding table
		::libmaus::huffman::HuffmanTree::EncodeTable const E(H);
		// get symbol array
		libmaus::autoarray::AutoArray<int64_t> const symbols = H.symbolArray();
		// get maximum symbol
		int64_t const maxsym = symbols.size() ? symbols[symbols.size()-1] : -1;
		// check it is in range
		assert ( 
			maxsym < 0
			||
			static_cast<uint64_t>(maxsym) <= static_cast<uint64_t>(std::numeric_limits<entity_type>::max()) 
		);
		// size of symbol table
		uint64_t const tablesize = maxsym+1;
		// number of inner nodes in Huffman tree
		uint64_t const inner = H.inner();

		::libmaus::huffman::IndexDecoderDataArray IDD(bwt);
		::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV = ::libmaus::huffman::IndexLoader::loadAccIndex(bwt);

		// compute symbol to node mapping
		uint64_t symtonodesvecsize = 0;
		// depth <= 16, tablesize <= 64k
		libmaus::autoarray::AutoArray<uint32_t> symtonodevecoffsets(tablesize,false);
		for ( uint64_t i = 0; i < symbols.size(); ++i )
		{
			uint64_t const sym = symbols[i];
			symtonodevecoffsets[sym] = symtonodesvecsize;
			assert ( symtonodesvecsize <= std::numeric_limits<uint32_t>::max() );
			symtonodesvecsize += E.getCodeLength(sym);
		}
		
		libmaus::autoarray::AutoArray<uint32_t> symtonodes(symtonodesvecsize,false);
		uint32_t * symtonodesp = symtonodes.begin();
		for ( uint64_t i = 0; i < symbols.size(); ++i )
		{
			// symbol
			uint64_t const sym = symbols[i];
			// node
			uint64_t node = H.root();
			
			assert ( symtonodesp-symtonodes.begin() == symtonodevecoffsets[sym] );
			
			for ( uint64_t j = 0; j < E.getCodeLength(sym); ++j )
			{
				// add symbol to inner node
				//symtonodes [ sym ] . push_back( node - H.leafs() );
				*(symtonodesp++) = node - H.leafs();
				
				// follow tree link according to code of symbol
				if ( E.getBitFromTop(sym,j) )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
		}
		assert ( symtonodesp = symtonodes.end() );
				
		// total size
		uint64_t const n = rl_decoder::getLength(bwt);
		// before
		uint64_t const n_0 = p0r;
		// terminator
		uint64_t const n_1 = 1;
		// after
		uint64_t const n_2 = n-(n_0+n_1);
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		assert ( numthreads );
		
		uint64_t const blocksperthread = 4;
		uint64_t const targetblocks = numthreads * blocksperthread;
		// maximum internal memory for blocks
		uint64_t const blockmemthres = 1024ull*1024ull;
		// per thread
		uint64_t const blocksizethres = (blockmemthres*sizeof(entity_type) + numthreads-1)/numthreads;
		
		// block sizes
		uint64_t const t_0 = std::min(blocksizethres,(n_0 + (targetblocks-1))/targetblocks);
		uint64_t const t_1 = std::min(blocksizethres,(n_1 + (targetblocks-1))/targetblocks);
		uint64_t const t_2 = std::min(blocksizethres,(n_2 + (targetblocks-1))/targetblocks);
		uint64_t const t_g = std::max(std::max(t_0,t_1),t_2);
		// actual number of blocks
		uint64_t const b_0 = t_0 ? ((n_0 + (t_0-1))/t_0) : 0;
		uint64_t const b_1 = t_1 ? ((n_1 + (t_1-1))/t_1) : 0;
		uint64_t const b_2 = t_2 ? ((n_2 + (t_2-1))/t_2) : 0;
		uint64_t const b_g = b_0 + b_1 + b_2;
		// blocks per thread
		uint64_t const blocks_per_thread_0 = (b_0 + numthreads-1)/numthreads;
		uint64_t const blocks_per_thread_2 = (b_2 + numthreads-1)/numthreads;
		
		// local character histograms
		libmaus::autoarray::AutoArray<uint64_t> localhist(numthreads * tablesize,false);
		// global node sizes
		libmaus::autoarray::AutoArray2d<uint64_t> localnodehist(inner,b_g+1,true);

		#if 0
		std::cerr << "(" << symtonodevecoffsets.byteSize() << "," << symtonodes.byteSize() << "," << localhist.byteSize() << "," 
			<< (inner * (b_g+1) * sizeof(uint64_t)) << ")";
		#endif
		
		libmaus::parallel::OMPLock cerrlock;
		
		libmaus::autoarray::AutoArray<rl_decoder::unique_ptr_type> rldecs(2*numthreads);
		Todo<RlDecoderInfoObject> todostack;
		// set up block information for symbols after terminator
		for ( uint64_t ii = 0; ii < numthreads; ++ii )
		{
			uint64_t const i = numthreads-ii-1;
			uint64_t const baseblock = i*blocks_per_thread_2;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_2, b_2);
			uint64_t const low = n_0+n_1+baseblock * t_2;
						
			if ( low < n )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif

				rldecs[numthreads+i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(RlDecoderInfoObject(
					rldecs[numthreads+i].get(),baseblock,highblock-baseblock,
					low,n_0+n_1+n_2,t_2,b_0+b_1
					)
				);
			}
		}
		// set up block information for symbols before terminator
		for ( uint64_t ii = 0; ii < numthreads; ++ii )
		{
			uint64_t const i = numthreads-ii-1;
			uint64_t const baseblock = i*blocks_per_thread_0;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_0, b_0);
			uint64_t const low = baseblock * t_0;
						
			if ( low < n_0 )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif
				rldecs[i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(RlDecoderInfoObject(
					rldecs[i].get(),baseblock,highblock-baseblock,
					low,n_0,t_0,0
					)
				);				
			}
		}
		
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			RlDecoderInfoObject dio;
			
			while ( todostack.pop(dio) )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif
			
				// pointer to histogram memory for this block
				uint64_t * const hist = localhist.begin() + tid*tablesize;
				// erase histogram
				std::fill(hist,hist+tablesize,0ull);

				// lower and upper bound of block in symbols
				uint64_t const low  = dio.getLow();
				uint64_t const high = dio.getHigh();
				
				assert ( high > low );
				
				rl_decoder & dec = *(dio.dec);

				// symbols to be processed
				uint64_t todo = high-low;
				std::pair<int64_t,uint64_t> R(-1,0);
				uint64_t toadd = 0;
					
				while ( todo )
				{
					// decode a run
					R = dec.decodeRun();
					// this should not be an end of file marker
					// assert ( R.first >= 0 );
					// count to be added
					toadd = std::min(todo,R.second);
					// add it
					hist [ R.first ] += toadd;
					// reduce todo
					todo -= toadd;					
				}
				
				if ( R.first != -1 && toadd != R.second )
					dec.putBack(std::pair<int64_t,uint64_t>(R.first,R.second-toadd));
				
				if ( dio.hasNextBlock() )
					todostack.push(dio.nextBlock());

				uint64_t const blockid = dio.getAbsoluteBlock();
				// compute bits per node in this block
				for ( uint64_t sym = 0; sym < tablesize; ++sym )
					if ( E.hasSymbol(sym) )
					{
						uint32_t * symtonodesp = symtonodes.begin() + symtonodevecoffsets[sym];
					
						for ( uint64_t i = 0; i < E.getCodeLength(sym); ++i )
							localnodehist(*(symtonodesp++) /* symtonodes[sym][i] */,blockid) += hist[sym];
							// localnodehist(symtonodes[sym][i],blockid) += hist[sym];
					}
			}
		}
		
		for ( uint64_t i = 0; i < rldecs.size(); ++i )
			rldecs[i].reset();
		
		// terminator
		for ( uint64_t i = 0; i < E.getCodeLength(bwtterm); ++i )
			// localnodehist ( symtonodes[bwtterm][i], b_0 ) += 1;
			localnodehist ( symtonodes[symtonodevecoffsets[bwtterm] + i], b_0 ) += 1;

		// compute prefix sums for each node
		localnodehist.prefixSums();
		
		#if 0
		for ( uint64_t node = 0; node < inner; ++node )
			for ( uint64_t b = 0; b < b_g+1; ++b )
				std::cerr << "localnodehist(" << node << "," << b << ")=" << localnodehist(node,b) << std::endl;
		#endif

		typedef libmaus::rank::ImpCacheLineRank rank_type;
		typedef rank_type::unique_ptr_type rank_ptr_type;
		libmaus::autoarray::AutoArray<rank_ptr_type> R(inner);
		libmaus::autoarray::AutoArray<uint64_t *> P(inner);
		libmaus::autoarray::AutoArray<entity_type> U(numthreads*t_g*2,false);

		for ( uint64_t node = 0; node < inner; ++node )
		{
			uint64_t const nodebits = localnodehist(node,b_g) + 1;
			uint64_t const nodedatawords = (nodebits+63)/64;
				
			rank_ptr_type tRnode(new rank_type(nodebits));
			R[node] = UNIQUE_PTR_MOVE(tRnode);
			P[node] = R[node]->A.end() - nodedatawords;
			
			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t i = 0; i < static_cast<int64_t>(nodedatawords); ++i )
				P[node][i] = 0;
		}
		
		libmaus::parallel::OMPLock nodelock;

		// set up block information for symbols after terminator
		for ( uint64_t ii = 0; ii < numthreads; ++ii )
		{
			uint64_t const i = numthreads-ii-1;
			uint64_t const baseblock = i*blocks_per_thread_2;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_2, b_2);
			uint64_t const low = n_0+n_1+baseblock * t_2;
						
			if ( low < n )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif

				rldecs[numthreads+i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(RlDecoderInfoObject(
					rldecs[numthreads+i].get(),baseblock,highblock-baseblock,
					low,n_0+n_1+n_2,t_2,b_0+b_1
					)
				);
			}
		}
		// set up block information for symbols before terminator
		for ( uint64_t ii = 0; ii < numthreads; ++ii )
		{
			uint64_t const i = numthreads-ii-1;
			uint64_t const baseblock = i*blocks_per_thread_0;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_0, b_0);
			uint64_t const low = baseblock * t_0;
						
			if ( low < n_0 )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif
				rldecs[i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(RlDecoderInfoObject(
					rldecs[i].get(),baseblock,highblock-baseblock,
					low,n_0,t_0,0
					)
				);				
			}
		}

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			RlDecoderInfoObject dio;
			
			while ( todostack.pop(dio) )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif

				#if 0
				nodelock.lock();
				std::cerr << "[" << tid << "] processing block " << dio.getAbsoluteBlock() << std::endl;
				nodelock.unlock();
				#endif
				
				// lower and upper bound of block in symbols
				uint64_t const low  = dio.getLow();
				uint64_t const high = dio.getHigh();
				
				assert ( high > low );
							
				entity_type * const block = U.begin() + (tid*t_g*2);
				entity_type * const tblock = block + t_g;
				
				// open decoder
				// rl_decoder dec(bwt,low);
				
				// symbols to be processed
				uint64_t todo = high-low;
				
				entity_type * tptr = block;

				rl_decoder & dec = *(dio.dec);

				// symbols to be processed
				std::pair<int64_t,uint64_t> R(-1,0);
				uint64_t toadd = 0;
				
				// load data
				while ( todo )
				{
					// decode a run
					R = dec.decodeRun();
					// this should not be an end of file marker
					// assert ( R.first >= 0 );
					// count to be added
					toadd = std::min(todo,R.second);
					// reduce todo
					todo -= toadd;
					
					for ( uint64_t i = 0; i < toadd; ++i )
						*(tptr++) = R.first;
				}

				if ( R.first != -1 && toadd != R.second )
					dec.putBack(std::pair<int64_t,uint64_t>(R.first,R.second-toadd));
				
				if ( dio.hasNextBlock() )
					todostack.push(dio.nextBlock());
				
				for ( uint64_t i = 0; i < (high-low); ++i )
				{
					uint64_t const sym = block[i];
					block[i] = E.getCode(sym) << (8*sizeof(entity_type)-E.getCodeLength(sym));
				}
				
				std::stack<WtTodo> todostack;
				todostack.push(WtTodo(0,high-low,0,H.root()));
				
				while ( todostack.size() )
				{
					WtTodo const wtcur = todostack.top();
					todostack.pop();
					
					#if 0
					std::cerr << wtcur << std::endl;
					#endif
					
					// inner node id
					uint64_t const unode = wtcur.node - H.root();
					// shift for relevant bit
					unsigned int const shift = (8*sizeof(entity_type)-1)-wtcur.depth;
					
					// bit offset
					uint64_t bitoff = localnodehist(unode,dio.getAbsoluteBlock());
					uint64_t towrite = wtcur.high-wtcur.low;
					entity_type * ptr = block+wtcur.low;
					
					entity_type * optrs[2] = { tblock, tblock+t_g-1 };
					static int64_t const inc[2] = { 1, -1 };
					uint64_t bcnt[2]; bcnt[0] = 0; bcnt[1] = 0;
				
					// align to word boundary	
					nodelock.lock();
					while ( towrite && ((bitoff%64)!=0) )
					{
						entity_type const sym = *(ptr++);
						uint64_t const bit = (sym >> shift) & 1;
						// assert ( bit < 2 );
						bcnt[bit]++;
						(*optrs[bit]) = sym;
						optrs[bit] += inc[bit];
					
						libmaus::bitio::putBit(P[unode],bitoff,bit);
					
						towrite--;
						bitoff++;
					}
					nodelock.unlock();
					
					// write full words
					uint64_t * PP = P[unode] + (bitoff/64);
					while ( towrite >= 64 )
					{
						uint64_t v = 0;
						for ( unsigned int i = 0; i < 64; ++i )
						{
							entity_type const sym = *(ptr++);
							uint64_t const bit = (sym >> shift) & 1;
							// assert ( bit < 2 );
							bcnt[bit]++;
							(*optrs[bit]) = sym;
							optrs[bit] += inc[bit];

							v <<= 1;
							v |= bit;
						}
					
						*(PP++) = v;
						towrite -= 64;
						bitoff += 64;
					}
					
					// write rest
					// assert ( towrite < 64 );
					nodelock.lock();
					while ( towrite )
					{
						entity_type const sym = *(ptr++);
						uint64_t const bit = (sym >> shift) & 1;
						// assert ( bit < 2 );
						bcnt[bit]++;
						(*optrs[bit]) = sym;
						optrs[bit] += inc[bit];
	
						libmaus::bitio::putBit(P[unode],bitoff,bit);
					
						towrite--;
						bitoff++;	
					}
					nodelock.unlock();					

					/* write data back */	
					ptr = block + wtcur.low;
					entity_type * inptr = tblock;
					for ( uint64_t i = 0; i < bcnt[0]; ++i )
						*(ptr++) = *(inptr++);

					inptr = tblock + t_g;
					for ( uint64_t i = 0; i < bcnt[1]; ++i )
						*(ptr++) = *(--inptr);
						
					// recursion
					if ( ( ! H.isLeaf( H.rightChild(wtcur.node) ) ) && (bcnt[1] != 0) )
						todostack.push(WtTodo(wtcur.high-bcnt[1],wtcur.high,wtcur.depth+1,H.rightChild(wtcur.node)));
					if ( (! H.isLeaf( H.leftChild(wtcur.node) )) && (bcnt[0] != 0) )
						todostack.push(WtTodo(wtcur.low,wtcur.low+bcnt[0],wtcur.depth+1,H.leftChild(wtcur.node)));
				}
			}
		}

		// reset decoders
		for ( uint64_t i = 0; i < numthreads; ++i )
			rldecs[i].reset();

		// terminator
		{
			uint64_t node = H.root();

			for ( uint64_t i = 0; i < E.getCodeLength(bwtterm); ++i )
			{
				uint64_t const bit = E.getBitFromTop(bwtterm,i);
				uint64_t const unode = node - H.root();				
				libmaus::bitio::putBit ( P[unode] , localnodehist ( unode , b_0 ), bit );
				
				if ( bit )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
		}

		// set up rank dictionaries
		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( uint64_t node = 0; node < inner; ++node )
		{
			uint64_t const nodebits = localnodehist(node,b_g) + 1;
			uint64_t const nodedatawords = (nodebits+63)/64;
			
			uint64_t * inptr = P[node];
			uint64_t * outptr = R[node]->A.begin();	
			
			uint64_t accb = 0;
			
			uint64_t todo = nodedatawords;
			
			while ( todo )
			{
				uint64_t const toproc = std::min(todo,static_cast<uint64_t>(6));
				
				uint64_t miniword = 0;
				uint64_t miniacc = 0;
				for ( uint64_t i = 0; i < toproc; ++i )
				{
					miniword |= miniacc << (i*9);
					miniacc  += libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(inptr[i]);
				}
				miniword |= (miniacc << (toproc*9));
				
				for ( uint64_t i = 0; i < toproc; ++i )
					outptr[2+i] = inptr[i];
				outptr[0] = accb;
				outptr[1] = miniword;
				
				accb += miniacc;
				inptr += toproc;
				outptr += 2+toproc;
				todo -= toproc;
			}
		}
		
		libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type pICHWT(
			new libmaus::wavelet::ImpCompactHuffmanWaveletTree(n,H,R)
		);
		
		#if 0
		libmaus::wavelet::ImpCompactHuffmanWaveletTree const & ICHWT = *pICHWT;
		// open decoder
		rl_decoder dec(bwt,0);
		libmaus::autoarray::AutoArray<uint64_t> ranktable(tablesize);
		for ( uint64_t i = 0; i < ICHWT.size(); ++i )
		{
			int64_t const sym = dec.decode();
			#if 0
			std::cerr << ICHWT[i] << "(" << sym << ")" << ";";
			#endif
						
			if ( i == p0r )
				assert ( ICHWT[i] == bwtterm );
			else
			{
				assert ( sym == ICHWT[i] );
				assert ( ranktable[sym] == ICHWT.rankm(sym,i) );
				assert ( ICHWT.select(sym,ranktable[sym]) == i );
				ranktable[sym]++;
				assert ( ranktable[sym] == ICHWT.rank(sym,i) );				
			}
		}
		#if 0
		std::cerr << std::endl;
		#endif
		#endif
		
		return UNIQUE_PTR_MOVE(pICHWT);
	}

	/**
	 * specialised version for small alphabets
	 **/
	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type rlToHwtSmallAlphabet(
		std::string const & bwt, 
		std::string const & huftreefilename
	)
	{
		// load the huffman tree
		::libmaus::huffman::HuffmanTree::unique_ptr_type UH = loadCompactHuffmanTree(huftreefilename);
		::libmaus::huffman::HuffmanTree & H = *UH;
		libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(rlToHwtSmallAlphabet(bwt,H));		
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	/**
	 * specialised version for small alphabets
	 **/
	template<typename entity_type>
	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type rlToHwtSmallAlphabet(
		std::string const & bwt, 
		::libmaus::huffman::HuffmanTree const & H
	)
	{
		// check depth is low
		assert ( H.maxDepth() <= 8*sizeof(entity_type) );
		// get encoding table
		::libmaus::huffman::HuffmanTree::EncodeTable const E(H);
		// get symbol array
		libmaus::autoarray::AutoArray<int64_t> const symbols = H.symbolArray();
		// get maximum symbol
		int64_t const maxsym = symbols.size() ? symbols[symbols.size()-1] : -1;
		// check it is in range
		assert (
			(maxsym < 0)
			|| 
			static_cast<uint64_t>(maxsym) <= static_cast<uint64_t>(std::numeric_limits<entity_type>::max()) 
		);
		// size of symbol table
		uint64_t const tablesize = maxsym+1;
		// number of inner nodes in Huffman tree
		uint64_t const inner = H.inner();

		// compute symbol to node mapping
		uint64_t symtonodesvecsize = 0;
		// depth <= 16, tablesize <= 64k
		libmaus::autoarray::AutoArray<uint32_t> symtonodevecoffsets(tablesize,false);
		for ( uint64_t i = 0; i < symbols.size(); ++i )
		{
			uint64_t const sym = symbols[i];
			symtonodevecoffsets[sym] = symtonodesvecsize;
			assert ( symtonodesvecsize <= std::numeric_limits<uint32_t>::max() );
			symtonodesvecsize += E.getCodeLength(sym);
		}
		
		libmaus::autoarray::AutoArray<uint32_t> symtonodes(symtonodesvecsize,false);
		uint32_t * symtonodesp = symtonodes.begin();
		for ( uint64_t i = 0; i < symbols.size(); ++i )
		{
			// symbol
			uint64_t const sym = symbols[i];
			// node
			uint64_t node = H.root();
			
			assert ( symtonodesp-symtonodes.begin() == symtonodevecoffsets[sym] );
			
			for ( uint64_t j = 0; j < E.getCodeLength(sym); ++j )
			{
				// add symbol to inner node
				//symtonodes [ sym ] . push_back( node - H.leafs() );
				*(symtonodesp++) = node - H.leafs();
				
				// follow tree link according to code of symbol
				if ( E.getBitFromTop(sym,j) )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
		}
		assert ( symtonodesp = symtonodes.end() );



		#if 0		
		// compute symbol to node mapping
		std::vector < std::vector < uint32_t > > symtonodes(tablesize);
		for ( uint64_t i = 0; i < symbols.size(); ++i )
		{
			// symbol
			uint64_t const sym = symbols[i];
			// node
			uint64_t node = H.root();
			
			for ( uint64_t j = 0; j < E.getCodeLength(sym); ++j )
			{
				// add symbol to inner node
				symtonodes [ sym ] . push_back( node - H.leafs() );
				
				// follow tree link according to code of symbol
				if ( E.getBitFromTop(sym,j) )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
		}
		#endif
		
		// total size
		uint64_t const n = rl_decoder::getLength(bwt);

		::libmaus::huffman::IndexDecoderDataArray IDD(std::vector<std::string>(1,bwt));
		::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV = 
			::libmaus::huffman::IndexLoader::loadAccIndex(std::vector<std::string>(1,bwt));
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		assert ( numthreads );
		
		uint64_t const blocksperthread = 4;
		uint64_t const targetblocks = numthreads * blocksperthread;
		// maximum internal memory for blocks
		uint64_t const blockmemthres = 1024ull*1024ull;
		// per thread
		uint64_t const blocksizethres = (blockmemthres*sizeof(entity_type) + numthreads-1)/numthreads;
		
		// block size
		uint64_t const t_g = std::min(blocksizethres,(n + (targetblocks-1))/targetblocks);
		// actual number of blocks
		uint64_t const b_g = t_g ? ((n + (t_g-1))/t_g) : 0;
		// blocks per thread
		uint64_t const blocks_per_thread_g = (b_g + numthreads-1)/numthreads;

		// local character histograms
		libmaus::autoarray::AutoArray<uint64_t> localhist(numthreads * tablesize,false);
		// global node sizes
		libmaus::autoarray::AutoArray2d<uint64_t> localnodehist(inner,b_g+1,true);

		#if 0
		std::cerr << "(" << symtonodevecoffsets.byteSize() << "," << symtonodes.byteSize() << "," << localhist.byteSize() << "," 
			<< (inner * (b_g+1) * sizeof(uint64_t)) << ")";
		#endif
		
		libmaus::parallel::OMPLock cerrlock;

		libmaus::autoarray::AutoArray<rl_decoder::unique_ptr_type> rldecs(numthreads);
		Todo<RlDecoderInfoObject> todostack;
		// set up block information
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			uint64_t const baseblock = i*blocks_per_thread_g;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_g, b_g);
			uint64_t const low = baseblock * t_g;
						
			if ( low < n )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif

				rldecs[i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(
					RlDecoderInfoObject(
						rldecs[i].get(),baseblock,highblock-baseblock,
						low,n,t_g,0
					)
				);
			}
		}
	
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			RlDecoderInfoObject dio;
			
			while ( todostack.pop(dio) )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif
			
				// pointer to histogram memory for this block
				uint64_t * const hist = localhist.begin() + tid*tablesize;
				// erase histogram
				std::fill(hist,hist+tablesize,0ull);

				// lower and upper bound of block in symbols
				uint64_t const low  = dio.getLow();
				uint64_t const high = dio.getHigh();
				
				assert ( high > low );

				rl_decoder & dec = *(dio.dec);

				// symbols to be processed
				uint64_t todo = high-low;
				std::pair<int64_t,uint64_t> R(-1,0);
				uint64_t toadd = 0;
				
				while ( todo )
				{
					// decode a run
					R = dec.decodeRun();
					// this should not be an end of file marker
					// assert ( R.first >= 0 );
					// count to be added
					toadd = std::min(todo,R.second);
					// add it
					hist [ R.first ] += toadd;
					// reduce todo
					todo -= toadd;					
				}

				if ( R.first != -1 && toadd != R.second )
					dec.putBack(std::pair<int64_t,uint64_t>(R.first,R.second-toadd));
				
				if ( dio.hasNextBlock() )
					todostack.push(dio.nextBlock());

				// compute bits per node in this block
				uint64_t const blockid = dio.getAbsoluteBlock();
				for ( uint64_t sym = 0; sym < tablesize; ++sym )
					if ( E.hasSymbol(sym) )
					{
						uint32_t * symtonodesp = symtonodes.begin() + symtonodevecoffsets[sym];
						
						for ( uint64_t i = 0; i < E.getCodeLength(sym) /* symtonodes[sym].size() */; ++i )
							localnodehist(*(symtonodesp++) /* symtonodes[sym][i] */,blockid) += hist[sym];
							// localnodehist(symtonodes[sym][i],blockid) += hist[sym];
					}
			}
		}

		for ( uint64_t i = 0; i < numthreads; ++i )
			rldecs[i].reset();
		
		// compute prefix sums for each node
		localnodehist.prefixSums();
		
		#if 0
		for ( uint64_t node = 0; node < inner; ++node )
			for ( uint64_t b = 0; b < b_g+1; ++b )
				std::cerr << "localnodehist(" << node << "," << b << ")=" << localnodehist(node,b) << std::endl;
		#endif

		typedef libmaus::rank::ImpCacheLineRank rank_type;
		typedef rank_type::unique_ptr_type rank_ptr_type;
		libmaus::autoarray::AutoArray<rank_ptr_type> R(inner);
		libmaus::autoarray::AutoArray<uint64_t *> P(inner);
		libmaus::autoarray::AutoArray<entity_type> U(numthreads*t_g*2,false);

		for ( uint64_t node = 0; node < inner; ++node )
		{
			uint64_t const nodebits = localnodehist(node,b_g) + 1;
			uint64_t const nodedatawords = (nodebits+63)/64;

			#if 0
			std::cerr << "node " << node << " bits " << localnodehist(node,b_g) << " words " <<
				(localnodehist(node,b_g)+63)/64 << std::endl;
			#endif
				
			rank_ptr_type tRnode(new rank_type(nodebits));
			R[node] = UNIQUE_PTR_MOVE(tRnode);
			P[node] = R[node]->A.end() - nodedatawords;
			
			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t i = 0; i < static_cast<int64_t>(nodedatawords); ++i )
				P[node][i] = 0;
		}
		
		libmaus::parallel::OMPLock nodelock;

		// set up block information
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			uint64_t const baseblock = i*blocks_per_thread_g;
			uint64_t const highblock = std::min(baseblock + blocks_per_thread_g, b_g);
			uint64_t const low = baseblock * t_g;
						
			if ( low < n )
			{
				#if defined(HUFRL)
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,low));
				#else
				rl_decoder::unique_ptr_type tptr(new rl_decoder(IDD,IECV.get(),low));
				#endif

				rldecs[i] = UNIQUE_PTR_MOVE(tptr);
				todostack.push(
					RlDecoderInfoObject(
						rldecs[i].get(),baseblock,highblock-baseblock,
						low,n,t_g,0
					)
				);
			}
		}

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			RlDecoderInfoObject dio;
			
			while ( todostack.pop(dio) )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif

				#if 0
				nodelock.lock();
				std::cerr << "[" << tid << "] processing block " << b << std::endl;
				nodelock.unlock();
				#endif
				
				// lower and upper bound of block in symbols
				uint64_t const low  = dio.getLow();
				uint64_t const high = dio.getHigh();
				
				assert ( high > low );

				entity_type * const block = U.begin() + (tid*t_g*2);
				entity_type * const tblock = block + t_g;
				
				entity_type * tptr = block;

				rl_decoder & dec = *(dio.dec);

				// symbols to be processed
				uint64_t todo = high-low;
				std::pair<int64_t,uint64_t> R(-1,0);
				uint64_t toadd = 0;

				// load data
				while ( todo )
				{
					// decode a run
					R = dec.decodeRun();
					// this should not be an end of file marker
					// assert ( R.first >= 0 );
					// count to be added
					toadd = std::min(todo,R.second);
					// reduce todo
					todo -= toadd;
					
					for ( uint64_t i = 0; i < toadd; ++i )
						*(tptr++) = R.first;
				}

				if ( R.first != -1 && toadd != R.second )
					dec.putBack(std::pair<int64_t,uint64_t>(R.first,R.second-toadd));
				
				if ( dio.hasNextBlock() )
					todostack.push(dio.nextBlock());
				
				for ( uint64_t i = 0; i < (high-low); ++i )
				{
					uint64_t const sym = block[i];
					block[i] = E.getCode(sym) << (8*sizeof(entity_type)-E.getCodeLength(sym));
				}
				
				std::stack<WtTodo> todostack;
				todostack.push(WtTodo(0,high-low,0,H.root()));
				
				while ( todostack.size() )
				{
					WtTodo const wtcur = todostack.top();
					todostack.pop();
					
					#if 0
					std::cerr << wtcur << std::endl;
					#endif
					
					// inner node id
					uint64_t const unode = wtcur.node - H.root();
					// shift for relevant bit
					unsigned int const shift = (8*sizeof(entity_type)-1)-wtcur.depth;
					
					// bit offset
					uint64_t bitoff = localnodehist(unode,dio.getAbsoluteBlock());
					uint64_t towrite = wtcur.high-wtcur.low;
					entity_type * ptr = block+wtcur.low;
					
					entity_type * optrs[2] = { tblock, tblock+t_g-1 };
					static int64_t const inc[2] = { 1, -1 };
					uint64_t bcnt[2]; bcnt[0] = 0; bcnt[1] = 0;
				
					// align to word boundary	
					nodelock.lock();
					while ( towrite && ((bitoff%64)!=0) )
					{
						entity_type const sym = *(ptr++);
						uint64_t const bit = (sym >> shift) & 1;
						// assert ( bit < 2 );
						bcnt[bit]++;
						(*optrs[bit]) = sym;
						optrs[bit] += inc[bit];
					
						libmaus::bitio::putBit(P[unode],bitoff,bit);
					
						towrite--;
						bitoff++;
					}
					nodelock.unlock();
					
					// write full words
					uint64_t * PP = P[unode] + (bitoff/64);
					while ( towrite >= 64 )
					{
						uint64_t v = 0;
						for ( unsigned int i = 0; i < 64; ++i )
						{
							entity_type const sym = *(ptr++);
							uint64_t const bit = (sym >> shift) & 1;
							// assert ( bit < 2 );
							bcnt[bit]++;
							(*optrs[bit]) = sym;
							optrs[bit] += inc[bit];

							v <<= 1;
							v |= bit;
						}
					
						*(PP++) = v;
						towrite -= 64;
						bitoff += 64;
					}
					
					// write rest
					// assert ( towrite < 64 );
					nodelock.lock();
					while ( towrite )
					{
						entity_type const sym = *(ptr++);
						uint64_t const bit = (sym >> shift) & 1;
						// assert ( bit < 2 );
						bcnt[bit]++;
						(*optrs[bit]) = sym;
						optrs[bit] += inc[bit];

						libmaus::bitio::putBit(P[unode],bitoff,bit);
					
						towrite--;
						bitoff++;	
					}
					nodelock.unlock();					

					/* write data back */	
					ptr = block + wtcur.low;
					entity_type * inptr = tblock;
					for ( uint64_t i = 0; i < bcnt[0]; ++i )
						*(ptr++) = *(inptr++);

					inptr = tblock + t_g;
					for ( uint64_t i = 0; i < bcnt[1]; ++i )
						*(ptr++) = *(--inptr);
						
					// recursion
					if ( ( ! H.isLeaf( H.rightChild(wtcur.node) ) ) && (bcnt[1] != 0) )
						todostack.push(WtTodo(wtcur.high-bcnt[1],wtcur.high,wtcur.depth+1,H.rightChild(wtcur.node)));
					if ( (! H.isLeaf( H.leftChild(wtcur.node) )) && (bcnt[0] != 0) )
						todostack.push(WtTodo(wtcur.low,wtcur.low+bcnt[0],wtcur.depth+1,H.leftChild(wtcur.node)));
				}
			}
		}

		for ( uint64_t i = 0; i < numthreads; ++i )
			rldecs[i].reset();

		// set up rank dictionaries
		for ( uint64_t node = 0; node < inner; ++node )
		{
			uint64_t const nodebits = localnodehist(node,b_g) + 1;
			uint64_t const nodedatawords = (nodebits+63)/64;
			
			uint64_t * inptr = P[node];
			uint64_t * outptr = R[node]->A.begin();	
			
			uint64_t accb = 0;
			
			uint64_t todo = nodedatawords;
			
			while ( todo )
			{
				uint64_t const toproc = std::min(todo,static_cast<uint64_t>(6));
				
				uint64_t miniword = 0;
				uint64_t miniacc = 0;
				for ( uint64_t i = 0; i < toproc; ++i )
				{
					miniword |= miniacc << (i*9);
					miniacc  += libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(inptr[i]);
				}
				miniword |= (miniacc << (toproc*9));
				
				for ( uint64_t i = 0; i < toproc; ++i )
					outptr[2+i] = inptr[i];
				outptr[0] = accb;
				outptr[1] = miniword;
				
				accb += miniacc;
				inptr += toproc;
				outptr += 2+toproc;
				todo -= toproc;
			}
		}
		
		libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type pICHWT(
			new libmaus::wavelet::ImpCompactHuffmanWaveletTree(n,H,R)
		);
		
		#if 0
		libmaus::wavelet::ImpCompactHuffmanWaveletTree const & ICHWT = *pICHWT;
		// open decoder
		rl_decoder dec(std::vector<std::string>(1,bwt),0);
		libmaus::autoarray::AutoArray<uint64_t> ranktable(tablesize);
		for ( uint64_t i = 0; i < ICHWT.size(); ++i )
		{
			int64_t const sym = dec.decode();
			std::cerr << ICHWT[i] << "(" << sym << ")" << ";";
						
			assert ( sym == ICHWT[i] );
			assert ( ranktable[sym] == ICHWT.rankm(sym,i) );
			assert ( ICHWT.select(sym,ranktable[sym]) == i );
			ranktable[sym]++;
			assert ( ranktable[sym] == ICHWT.rank(sym,i) );				
		}

		std::cerr << std::endl;

		#endif
		
		return UNIQUE_PTR_MOVE(pICHWT);
	}

	static ::libmaus::util::Histogram::unique_ptr_type computeRlSymHist(std::string const & bwt)
	{
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif

		uint64_t const n = rl_decoder::getLength(bwt);
		
		uint64_t const numpacks = 4*numthreads;
		uint64_t const packsize = (n + numpacks - 1)/numpacks;
		::libmaus::parallel::OMPLock histlock;
		::libmaus::util::Histogram::unique_ptr_type mhist(new ::libmaus::util::Histogram);

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t t = 0; t < static_cast<int64_t>(numpacks); ++t )
		{
			uint64_t const low = std::min(t*packsize,n);
			uint64_t const high = std::min(low+packsize,n);
			::libmaus::util::Histogram lhist;
			
			if ( high-low )
			{
				rl_decoder dec(std::vector<std::string>(1,bwt),low);
				
				uint64_t todec = high-low;
				std::pair<int64_t,uint64_t> P;
				
				while ( todec )
				{
					P = dec.decodeRun();
					assert ( P.first >= 0 );
					P.second = std::min(P.second,todec);
					lhist.add(P.first,P.second);
					
					todec -= P.second;
				}
			}
			
			histlock.lock();
			mhist->merge(lhist);
			histlock.unlock();
		}

		return UNIQUE_PTR_MOVE(mhist);
	}
	
	static unsigned int utf8WaveletMaxThreads()
	{
		return 24;
	}
	
	static uint64_t utf8WaveletMaxPartMem()
	{
		return 64ull*1024ull*1024ull;
	}

	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type rlToHwt(
		std::string const & bwt, 
		std::string const & hwt, 
		std::string const tmpprefix
	)
	{
		::libmaus::util::Histogram::unique_ptr_type mhist = UNIQUE_PTR_MOVE(computeRlSymHist(bwt));
		::std::map<int64_t,uint64_t> const chist = mhist->getByType<int64_t>();

		::libmaus::huffman::HuffmanTree H ( chist.begin(), chist.size(), false, true, true );
		
		if ( utf8Wavelet() )
		{
			::libmaus::wavelet::Utf8ToImpCompactHuffmanWaveletTree::constructWaveletTreeFromRl<rl_decoder,true /* radix sort */>(
				bwt,hwt,tmpprefix,H,
				utf8WaveletMaxPartMem() /* part size maximum */,
				utf8WaveletMaxThreads());

			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type IHWT(libmaus::wavelet::ImpCompactHuffmanWaveletTree::load(hwt));
			return UNIQUE_PTR_MOVE(IHWT);

			#if 0
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type IHWT(libmaus::wavelet::ImpCompactHuffmanWaveletTree::load(hwt));
			libmaus::lf::ImpCompactHuffmanWaveletLF IHWL(hwt);
			rl_decoder rldec(std::vector<std::string>(1,bwt));
			std::cerr << "Checking output bwt of length " << IHWT->size() << "...";
			for ( uint64_t i = 0; i < IHWT->size(); ++i )
			{
				uint64_t const sym = rldec.get();
				assert ( sym == (*IHWT)[i] );
				assert ( sym == IHWL[i] );
			}
			std::cerr << "done." << std::endl;
			
			std::cerr << "Running inverse select loop...";
			uint64_t r = 0;
			for ( uint64_t i = 0; i < IHWT->size(); ++i )
			{
				IHWT->inverseSelect(i);
				// r = IHWL(r);
			}
			std::cerr << "done." << std::endl;
			#endif
		}
		else
		{
			// special case for very small alphabets
			if ( H.maxDepth() <= 8*sizeof(uint8_t) && H.maxSymbol() <= std::numeric_limits<uint8_t>::max() )
			{
				libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(rlToHwtSmallAlphabet<uint8_t>(bwt,H));
				libmaus::aio::CheckedOutputStream COS(hwt);
				ptr->serialise(COS);
				COS.flush();
				COS.close();
				return UNIQUE_PTR_MOVE(ptr);
			}
			else if ( H.maxDepth() <= 8*sizeof(uint16_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
			{
				libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(rlToHwtSmallAlphabet<uint16_t>(bwt,H));
				libmaus::aio::CheckedOutputStream COS(hwt);
				ptr->serialise(COS);
				COS.flush();
				COS.close();
				return UNIQUE_PTR_MOVE(ptr);
			}
			else if ( H.maxDepth() <= 8*sizeof(uint32_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
			{
				libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(rlToHwtSmallAlphabet<uint32_t>(bwt,H));
				libmaus::aio::CheckedOutputStream COS(hwt);
				ptr->serialise(COS);
				COS.flush();
				COS.close();
				return UNIQUE_PTR_MOVE(ptr);
			}
			else if ( H.maxDepth() <= 8*sizeof(uint64_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
			{
				libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(rlToHwtSmallAlphabet<uint64_t>(bwt,H));
				libmaus::aio::CheckedOutputStream COS(hwt);
				ptr->serialise(COS);
				COS.flush();
				COS.close();
				return UNIQUE_PTR_MOVE(ptr);
			}
			else
			{
				::libmaus::util::TempFileNameGenerator tmpgen(tmpprefix,3);

				#if defined(_OPENMP)
				uint64_t const numthreads = omp_get_max_threads();
				#else
				uint64_t const numthreads = 1;
				#endif

				uint64_t const n = rl_decoder::getLength(bwt);
				uint64_t const packsize = (n + numthreads - 1)/numthreads;
				uint64_t const numpacks = (n + packsize-1)/packsize;

				::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffmanParallel IEWGH(H,tmpgen,numthreads);

				#if defined(_OPENMP)
				#pragma omp parallel for
				#endif
				for ( int64_t t = 0; t < static_cast<int64_t>(numpacks); ++t )
				{
					assert ( t < static_cast<int64_t>(numthreads) );
					uint64_t const low = std::min(t*packsize,n);
					uint64_t const high = std::min(low+packsize,n);
					
					::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffmanParallel::BufferType & BTS = IEWGH[t];

					if ( high-low )
					{
						rl_decoder dec(std::vector<std::string>(1,bwt),low);

						uint64_t todec = high-low;
						std::pair<int64_t,uint64_t> P;
						
						while ( todec )
						{
							P = dec.decodeRun();
							assert ( P.first >= 0 );
							P.second = std::min(P.second,todec);

							for ( uint64_t i = 0; i < P.second; ++i )
								BTS.putSymbol(P.first);
							
							todec -= P.second;
						}
					}
				}
				
				IEWGH.createFinalStream(hwt);

				libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type IHWT(libmaus::wavelet::ImpCompactHuffmanWaveletTree::load(hwt));
				return UNIQUE_PTR_MOVE(IHWT);
			}
		}
	}
	
	static void rlToHwtTerm(
		std::vector<std::string> const & bwt, 
		std::string const & hwt, 
		std::string const tmpprefix,
		::libmaus::huffman::HuffmanTree & H,
		uint64_t const bwtterm,
		uint64_t const p0r
		)
	{
		if ( utf8Wavelet() )
		{
			::libmaus::wavelet::Utf8ToImpCompactHuffmanWaveletTree::constructWaveletTreeFromRlWithTerm<rl_decoder,true /* radix sort */>(
				bwt,hwt,tmpprefix,H,p0r,bwtterm,
				utf8WaveletMaxPartMem() /* maximum part size */,
				utf8WaveletMaxThreads()
			);
			
			#if 0
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type IHWT(libmaus::wavelet::ImpCompactHuffmanWaveletTree::load(hwt));
			rl_decoder rldec(bwt);

			std::cerr << "Checking output bwt of length " << IHWT->size() << "...";
			for ( uint64_t i = 0; i < IHWT->size(); ++i )
			{
				uint64_t const sym = rldec.get();
				assert ( sym == (*IHWT)[i] );
			}
			std::cerr << "done." << std::endl;
			#endif
		}
		else
		{			
			::libmaus::util::TempFileNameGenerator tmpgen(tmpprefix,3);

			#if defined(_OPENMP)
			uint64_t const numthreads = omp_get_max_threads();
			#else
			uint64_t const numthreads = 1;
			#endif

			uint64_t const n = rl_decoder::getLength(bwt);

			assert ( p0r < n );
			uint64_t const nlow = p0r;
			uint64_t const nhigh = n - (nlow + 1);
			
			uint64_t const packsizelow  = (nlow + numthreads - 1)/numthreads;
			uint64_t const numpackslow  = packsizelow ? ( (nlow + packsizelow-1)/packsizelow ) : 0;
			uint64_t const packsizehigh = (nhigh + numthreads - 1)/numthreads;
			uint64_t const numpackshigh = packsizehigh ? ( (nhigh + packsizehigh-1)/packsizehigh ) : 0;

			::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffmanParallel IEWGH(H,tmpgen,2*numthreads+1);

			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t t = 0; t < static_cast<int64_t>(numpackslow); ++t )
			{
				assert ( t < static_cast<int64_t>(numthreads) );
				uint64_t const low  = std::min(t*packsizelow,nlow);
				uint64_t const high = std::min(low+packsizelow,nlow);
				
				::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffmanParallel::BufferType & BTS = IEWGH[t];

				if ( high-low )
				{
					rl_decoder dec(bwt,low);

					uint64_t todec = high-low;
					std::pair<int64_t,uint64_t> P;
					
					while ( todec )
					{
						P = dec.decodeRun();
						assert ( P.first >= 0 );
						P.second = std::min(P.second,todec);

						for ( uint64_t i = 0; i < P.second; ++i )
							BTS.putSymbol(P.first);
						
						todec -= P.second;
					}
				}
			}
			
			IEWGH[numthreads].putSymbol(bwtterm);

			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t t = 0; t < static_cast<int64_t>(numpackshigh); ++t )
			{
				assert ( t < static_cast<int64_t>(numthreads) );
				uint64_t const low  = std::min(nlow + 1 + t*packsizehigh,n);
				uint64_t const high = std::min(low+packsizehigh,n);
				
				::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffmanParallel::BufferType & BTS = IEWGH[numthreads+1+t];

				if ( high-low )
				{
					rl_decoder dec(bwt,low);

					uint64_t todec = high-low;
					std::pair<int64_t,uint64_t> P;
					
					while ( todec )
					{
						P = dec.decodeRun();
						assert ( P.first >= 0 );
						P.second = std::min(P.second,todec);

						for ( uint64_t i = 0; i < P.second; ++i )
							BTS.putSymbol(P.first);
						
						todec -= P.second;
					}
				}
			}
			
			IEWGH.createFinalStream(hwt);
		}
	}
	
	static void rlToHwtTerm(
		std::vector<std::string> const & bwt, 
		std::string const & hwt, 
		std::string const tmpprefix,
		::std::map<int64_t,uint64_t> const & chist,
		uint64_t const bwtterm,
		uint64_t const p0r
		)
	{
		::libmaus::huffman::HuffmanTree H(chist.begin(),chist.size(),false,true,true);
		rlToHwtTerm(bwt,hwt,tmpprefix,H,bwtterm,p0r);
	}

	static ::libmaus::huffman::HuffmanTree::unique_ptr_type loadCompactHuffmanTree(std::string const & huftreefilename)
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type CIN(new libmaus::aio::CheckedInputStream(huftreefilename));
		::libmaus::huffman::HuffmanTree::unique_ptr_type tH(new ::libmaus::huffman::HuffmanTree(*CIN));
		CIN->close();
		CIN.reset();
		
		return UNIQUE_PTR_MOVE(tH);
	}

	static ::libmaus::huffman::HuffmanTreeNode::shared_ptr_type loadHuffmanTree(std::string const & huftreefilename)
	{
		// deserialise symbol frequences
		libmaus::aio::CheckedInputStream::unique_ptr_type chistCIN(new libmaus::aio::CheckedInputStream(huftreefilename));
		::libmaus::huffman::HuffmanTreeNode::shared_ptr_type shnode = 
			::libmaus::huffman::HuffmanTreeNode::deserialize(*chistCIN);
		chistCIN->close();
		chistCIN.reset();
		
		return shnode;
	}

	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type rlToHwtTerm(
		std::vector<std::string> const & bwt, 
		std::string const & hwt, 
		std::string const tmpprefix,
		std::string const huftreefilename,
		uint64_t const bwtterm,
		uint64_t const p0r
		)
	{
		::libmaus::huffman::HuffmanTree::unique_ptr_type UH = loadCompactHuffmanTree(huftreefilename);
		::libmaus::huffman::HuffmanTree & H = *UH;

		// std::cerr << "(maxdepth=" << H.maxDepth() << ",maxSymbol=" << H.maxSymbol() << ")";

		if ( H.maxDepth() <= 8*sizeof(uint8_t) && H.maxSymbol() <= std::numeric_limits<uint8_t>::max() )
		{
			// std::cerr << "(small)";
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(rlToHwtTermSmallAlphabet<uint8_t>(bwt,huftreefilename,bwtterm,p0r));			
			return UNIQUE_PTR_MOVE(tICHWT);
		}
		else if ( H.maxDepth() <= 8*sizeof(uint16_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
		{
			// std::cerr << "(small)";
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(rlToHwtTermSmallAlphabet<uint16_t>(bwt,huftreefilename,bwtterm,p0r));			
			return UNIQUE_PTR_MOVE(tICHWT);
		}
		else if ( H.maxDepth() <= 8*sizeof(uint32_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
		{
			// std::cerr << "(small)";
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(rlToHwtTermSmallAlphabet<uint32_t>(bwt,huftreefilename,bwtterm,p0r));			
			return UNIQUE_PTR_MOVE(tICHWT);
		}
		else if ( H.maxDepth() <= 8*sizeof(uint64_t) && H.maxSymbol() <= std::numeric_limits<uint16_t>::max() )
		{
			// std::cerr << "(small)";
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(rlToHwtTermSmallAlphabet<uint64_t>(bwt,huftreefilename,bwtterm,p0r));			
			return UNIQUE_PTR_MOVE(tICHWT);
		}
		else
		{
		
			// std::cerr << "(large)";
			rlToHwtTerm(bwt,hwt,tmpprefix,H,bwtterm,p0r);
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(loadWaveletTree(hwt));
			return UNIQUE_PTR_MOVE(tICHWT);
		}
	}

};

struct RlToHwtTermRequest
{
	typedef RlToHwtTermRequest this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;
	
	std::vector<std::string> bwt;
	std::string hwt;
	std::string tmpprefix;
	std::string huftreefilename;
	uint64_t bwtterm;
	uint64_t p0r;
	bool utf8;
	
	RlToHwtTermRequest() {}
	RlToHwtTermRequest(
		std::vector<std::string> const & rbwt,
		std::string const & rhwt,
		std::string const & rtmpprefix,
		std::string const & rhuftreefilename,
		uint64_t const rbwtterm,
		uint64_t const rp0r,
		bool const rutf8
	) : bwt(rbwt), hwt(rhwt), tmpprefix(rtmpprefix), huftreefilename(rhuftreefilename), bwtterm(rbwtterm), p0r(rp0r), utf8(rutf8) {}
	
	RlToHwtTermRequest(std::istream & in)
	:
		bwt(libmaus::util::StringSerialisation::deserialiseStringVector(in)),
		hwt(libmaus::util::StringSerialisation::deserialiseString(in)),
		tmpprefix(libmaus::util::StringSerialisation::deserialiseString(in)),
		huftreefilename(libmaus::util::StringSerialisation::deserialiseString(in)),
		bwtterm(libmaus::util::NumberSerialisation::deserialiseSignedNumber(in)),
		p0r(libmaus::util::NumberSerialisation::deserialiseNumber(in)),
		utf8(libmaus::util::NumberSerialisation::deserialiseNumber(in))
	{
	
	}
	
	static unique_ptr_type load(std::string const & filename)
	{
		libmaus::aio::CheckedInputStream CIS(filename);
		unique_ptr_type ptr(new this_type(CIS));
		return UNIQUE_PTR_MOVE(ptr);
	}
	
	std::ostream & serialise(std::ostream & out) const
	{
		libmaus::util::StringSerialisation::serialiseStringVector(out,bwt);
		libmaus::util::StringSerialisation::serialiseString(out,hwt);
		libmaus::util::StringSerialisation::serialiseString(out,tmpprefix);
		libmaus::util::StringSerialisation::serialiseString(out,huftreefilename);
		libmaus::util::NumberSerialisation::serialiseSignedNumber(out,bwtterm);
		libmaus::util::NumberSerialisation::serialiseNumber(out,p0r);
		libmaus::util::NumberSerialisation::serialiseNumber(out,utf8);
		return out;
	}

	static std::ostream & serialise(
		std::ostream & out,
		std::vector<std::string> const & bwt,
		std::string const & hwt,
		std::string const & tmpprefix,
		std::string const & huftreefilename,
		uint64_t const bwtterm,
		uint64_t const p0r,
		bool const utf8
	)
	{
		libmaus::util::StringSerialisation::serialiseStringVector(out,bwt);
		libmaus::util::StringSerialisation::serialiseString(out,hwt);
		libmaus::util::StringSerialisation::serialiseString(out,tmpprefix);
		libmaus::util::StringSerialisation::serialiseString(out,huftreefilename);
		libmaus::util::NumberSerialisation::serialiseSignedNumber(out,bwtterm);
		libmaus::util::NumberSerialisation::serialiseNumber(out,p0r);
		libmaus::util::NumberSerialisation::serialiseNumber(out,utf8);
		return out;
	}

	libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type dispatch()
	{
		if ( utf8 )
		{
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tptr(
				RlToHwtBase<true>::rlToHwtTerm(bwt,hwt,tmpprefix,huftreefilename,bwtterm,p0r)
			);
			return UNIQUE_PTR_MOVE(tptr);
		}
		else
		{
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tptr(
				RlToHwtBase<false>::rlToHwtTerm(bwt,hwt,tmpprefix,huftreefilename,bwtterm,p0r)
			);
			return UNIQUE_PTR_MOVE(tptr);
		}
	}
};

struct BwtMergeBlockSortRequest : libmaus::suffixsort::BwtMergeEnumBase
{
	bwt_merge_sort_input_type inputtype;
	std::string fn; // file name of complete file
	uint64_t fs; // size of complete file
	std::string chistfilename; // file name of global character histogram
	std::string huftreefilename; // file name of global huffman tree
	uint64_t bwtterm; // bwt term symbol (!= any other symbol)
	uint64_t maxsym; // maximal appearing symbol
	std::string tmpfilenamesser; // temp file names for this block
	std::string tmpfilenamebase; // temp file name base for files temporary to the computation of this block
	uint64_t rlencoderblocksize; // block size for run length encoder
	uint64_t isasamplingrate; // sampling rate for inverse suffix array
	uint64_t blockstart; // start of this block (in symbols)
	uint64_t cblocksize; // size of this block (in symbols)
	::libmaus::suffixsort::BwtMergeZBlockRequestVector zreqvec; // vector of positions in file where rank in this block is requested
	bool computeTermSymbolHwt;
	
	static bwt_merge_sort_input_type decodeInputType(uint64_t const i)
	{
		switch ( i )
		{
			case bwt_merge_input_type_byte:
				return bwt_merge_input_type_byte;
			case bwt_merge_input_type_compact:
				return bwt_merge_input_type_compact;
			case bwt_merge_input_type_pac:
				return bwt_merge_input_type_pac;
			case bwt_merge_input_type_pac_term:
				return bwt_merge_input_type_pac_term;
			case bwt_merge_input_type_utf8:
				return bwt_merge_input_type_utf8;
			default:
			{
				::libmaus::exception::LibMausException ex;
				ex.getStream() << "Number " << i << " is not a valid input type designator." << std::endl;
				ex.finish();
				throw ex;
			}
		}
	}
	
	template<typename stream_type>
	void serialise(stream_type & stream) const
	{
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,static_cast<int>(inputtype));		
		::libmaus::util::StringSerialisation::serialiseString(stream,fn);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,fs);
		::libmaus::util::StringSerialisation::serialiseString(stream,chistfilename);
		::libmaus::util::StringSerialisation::serialiseString(stream,huftreefilename);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,bwtterm);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,maxsym);
		::libmaus::util::StringSerialisation::serialiseString(stream,tmpfilenamesser);
		::libmaus::util::StringSerialisation::serialiseString(stream,tmpfilenamebase);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,rlencoderblocksize);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,isasamplingrate);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,blockstart);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,cblocksize);
		zreqvec.serialise(stream);
		::libmaus::util::NumberSerialisation::serialiseNumber(stream,computeTermSymbolHwt);
	}
	
	std::string serialise() const
	{
		std::ostringstream ostr;
		serialise(ostr);
		return ostr.str();
	}
	
	BwtMergeBlockSortRequest()
	{
	}

	template<typename stream_type>	
	BwtMergeBlockSortRequest(stream_type & stream)
	:
		inputtype(decodeInputType(::libmaus::util::NumberSerialisation::deserialiseNumber(stream))),
		fn(::libmaus::util::StringSerialisation::deserialiseString(stream)),
		fs(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		chistfilename(::libmaus::util::StringSerialisation::deserialiseString(stream)),
		huftreefilename(::libmaus::util::StringSerialisation::deserialiseString(stream)),
		bwtterm(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		maxsym(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		tmpfilenamesser(::libmaus::util::StringSerialisation::deserialiseString(stream)),
		tmpfilenamebase(::libmaus::util::StringSerialisation::deserialiseString(stream)),
		rlencoderblocksize(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		isasamplingrate(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		blockstart(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		cblocksize(::libmaus::util::NumberSerialisation::deserialiseNumber(stream)),
		zreqvec(stream),
		computeTermSymbolHwt(::libmaus::util::NumberSerialisation::deserialiseNumber(stream))
	{
	}

	BwtMergeBlockSortRequest(
		bwt_merge_sort_input_type rinputtype,
		std::string rfn,
		uint64_t rfs,
		std::string rchistfilename,
		std::string rhuftreefilename,
		uint64_t rbwtterm,
		uint64_t rmaxsym,
		std::string rtmpfilenamesser,
		std::string rtmpfilenamebase,
		uint64_t rrlencoderblocksize,
		uint64_t risasamplingrate,
		uint64_t rblockstart,
		uint64_t rcblocksize,
		::libmaus::suffixsort::BwtMergeZBlockRequestVector const & rzreqvec,
		bool const rcomputeTermSymbolHwt
	)
	: 
		inputtype(rinputtype),
		fn(rfn),
		fs(rfs),
		chistfilename(rchistfilename),
		huftreefilename(rhuftreefilename),
		bwtterm(rbwtterm),
		maxsym(rmaxsym),
		tmpfilenamesser(rtmpfilenamesser),
		tmpfilenamebase(rtmpfilenamebase),
		rlencoderblocksize(rrlencoderblocksize),
		isasamplingrate(risasamplingrate),
		blockstart(rblockstart),
		cblocksize(rcblocksize),
		zreqvec(rzreqvec),
		computeTermSymbolHwt(rcomputeTermSymbolHwt)
	{
	}
	
	static BwtMergeBlockSortRequest load(std::string const & s)
	{
		std::istringstream istr(s);
		return BwtMergeBlockSortRequest(istr);
	}

	template<typename input_types_type>
	static uint64_t findSplitCommon(
		std::string const & fn,
		// position of textblock
		uint64_t const t,
		// length of textblock
		uint64_t const n,
		// position of pattern
		uint64_t const p,
		// length of file
		uint64_t const m
	)
	{
		typedef typename input_types_type::base_input_stream base_input_stream;
		typedef typename input_types_type::circular_wrapper circular_wrapper;
		
		circular_wrapper textstr(fn,t);
		circular_wrapper patstr(fn,p);
		
		// dynamically growing best prefix table
		::libmaus::util::KMP::BestPrefix<base_input_stream> BP(patstr,m);
		// adapter for accessing pattern in BP
		typename ::libmaus::util::KMP::BestPrefix<base_input_stream>::BestPrefixXAdapter xadapter = BP.getXAdapter();
		// call KMP adaption
		std::pair<uint64_t, uint64_t> Q = ::libmaus::util::KMP::PREFIX_SEARCH_INTERNAL_RESTRICTED(
			// pattern
			xadapter,m,BP,
			// text
			textstr,m,
			// restriction for position
			n
		);
		
		return Q.second;
	}
	
	static ::libmaus::huffman::HuffmanTree::unique_ptr_type loadCompactHuffmanTree(std::string const & huftreefilename)
	{
		libmaus::aio::CheckedInputStream::unique_ptr_type CIN(new libmaus::aio::CheckedInputStream(huftreefilename));
		::libmaus::huffman::HuffmanTree::unique_ptr_type tH(new ::libmaus::huffman::HuffmanTree(*CIN));
		CIN->close();
		CIN.reset();
		
		return UNIQUE_PTR_MOVE(tH);
	}

	static ::libmaus::huffman::HuffmanTreeNode::shared_ptr_type loadHuffmanTree(std::string const & huftreefilename)
	{
		// deserialise symbol frequences
		libmaus::aio::CheckedInputStream::unique_ptr_type chistCIN(new libmaus::aio::CheckedInputStream(huftreefilename));
		::libmaus::huffman::HuffmanTreeNode::shared_ptr_type shnode = 
			::libmaus::huffman::HuffmanTreeNode::deserialize(*chistCIN);
		chistCIN->close();
		chistCIN.reset();
		
		return shnode;
	}

	template<typename input_types_type>
	::libmaus::suffixsort::BwtMergeBlockSortResult sortBlock() const
	{
		// typedef typename input_types_type::base_input_stream base_input_stream;
		// typedef typename base_input_stream::char_type char_type;
		// typedef typename ::libmaus::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
	
		// glock.lock();
	
		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB1] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
	
		std::ostringstream tmpfilenamedirstr;
		tmpfilenamedirstr 
			<< tmpfilenamebase << "_sortblock_" 
			<< std::setw(10) << std::setfill('0') << blockstart
			<< std::setw(0) << "_"
			<< std::setw(10) << std::setfill('0') << cblocksize
			;
		std::string const tmpfilenamedir = tmpfilenamedirstr.str();
		::libmaus::util::TempFileNameGenerator tmpgen(tmpfilenamedir,3);

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB2] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		::libmaus::suffixsort::BwtMergeTempFileNameSet const tmpfilenames = 
			::libmaus::suffixsort::BwtMergeTempFileNameSet::load(tmpfilenamesser);

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB3] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// result
		::libmaus::suffixsort::BwtMergeBlockSortResult result;
		// copy request values
		result.setBlockStart( blockstart );
		result.setCBlockSize( cblocksize );
		result.setTempFileSet( tmpfilenames );

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB4] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// set up huffman tree
		unsigned int const albits = maxsym ? (8*sizeof(uint64_t) - ::libmaus::bitio::Clz::clz(maxsym)) : 0;
		
		// symbol before block
		int64_t const presym = input_types_type::linear_wrapper::getSymbolAtPosition(fn,blockstart ? (blockstart-1) : (fs-1));

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB5] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// start of next block
		uint64_t const nextblockstart = (blockstart + cblocksize) % fs;
		
		// find lcp between this block and start of next
		uint64_t const blcp = findSplitCommon<input_types_type>(fn,blockstart,cblocksize,nextblockstart,fs);

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB6] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// size of input string we need to read
		uint64_t const readsize = (cblocksize + blcp + 1);
		
		// 
		typedef typename input_types_type::string_type string_type;
		typedef typename input_types_type::circular_wrapper circular_wrapper;
		typedef typename circular_wrapper::unique_ptr_type circular_wrapper_ptr_type;
		circular_wrapper_ptr_type cwptr(new circular_wrapper(fn,blockstart));
		uint64_t const octetlength = string_type::computeOctetLength(*cwptr,readsize);
		cwptr.reset();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB7] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// set up circular reader
		typename input_types_type::circular_wrapper circ(fn,blockstart);
		cwptr = UNIQUE_PTR_MOVE(circular_wrapper_ptr_type(new circular_wrapper(fn,blockstart)));
		// construct string (read text and preprocess it for random symbol access)
		typename string_type::unique_ptr_type PT(new string_type(*cwptr, octetlength, readsize));
		string_type & T = *PT;

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB8] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		typedef typename string_type::saidx_t saidx_t;
		::libmaus::autoarray::AutoArray<saidx_t, ::libmaus::autoarray::alloc_type_c> SA =
			T.computeSuffixArray32();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB9] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// compute character histogram
		::libmaus::util::Histogram hist;
		for ( uint64_t i = 0; i < cblocksize; ++i )
			hist ( T[i] );

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB10] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		std::map<int64_t,uint64_t> const histm = hist.getByType<int64_t>();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB11] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		::libmaus::lf::DArray D(histm,bwtterm);
		D.serialise(tmpfilenames.getHist());

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB12] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// check whether first suffix of next block is smaller or larger than first suffix of this block
		bool gtlast = false;
		for ( uint64_t i = 0; i < SA.size(); ++i )
			// position 0 comes first, first suffix of next block is larger than first suffix of this block
			if ( !SA[i] )
			{
				gtlast = true;
				break;
			}
			// first suffix of next block comes first
			else if ( SA[i] == static_cast<saidx_t>(cblocksize) )
			{
				gtlast = false;
				break;
			}

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB13] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// remove terminator symbols from suffix array
		saidx_t * out = SA.begin();
		for ( saidx_t * in = SA.begin(); in != SA.end(); ++in )
			if ( *in < static_cast<saidx_t>(cblocksize) )
				*(out++) = *in;
		assert ( out-SA.begin() == static_cast<ptrdiff_t>(cblocksize) );

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB14] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// search for rank of first position in block
		for ( saidx_t * in = SA.begin(); in != out; ++in )
			if ( ! *in )
				result.setBlockP0Rank( (in-SA.begin()) );

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB15] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// check if we find the same via binary search		
		assert ( result.getBlockP0Rank() == input_types_type::circular_suffix_comparator::suffixSearch(SA.begin(), cblocksize, blockstart /* offset */, blockstart, fn, fs) );

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB16] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// search for rank of first position in complete file
		// result.absp0rank = ::libmaus::suffixsort::CircularSuffixComparator::suffixSearch(SA.begin(), cblocksize, blockstart, 0, fn, fs);
		
		// store sampled inverse suffix array
		assert ( ::libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(isasamplingrate) == 1 );
		uint64_t const isasamplingmask = isasamplingrate-1;
		::libmaus::aio::SynchronousGenericOutput<uint64_t> SGOISA(tmpfilenames.getSampledISA(),16*1024);
		for ( uint64_t r = 0; r < cblocksize; ++r )
		{
			uint64_t const p = SA[r];
			
			if ( ! (p & isasamplingmask) )
			{
				SGOISA.put(r);
				SGOISA.put(p + blockstart);
			}
		}
		SGOISA.flush();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB17] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		/**
		 * compute ranks for lf mapping blocks
		 **/
		// std::cerr << "[V] searching for " << zreqvec.size() << " suffixes...";
		::libmaus::timing::RealTimeClock sufsertc; sufsertc.start();
		result.resizeZBlocks(zreqvec.size());
		#if defined(_OPENMP)
		// #pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( uint64_t z = 0; z < zreqvec.size(); ++z )
		{
			uint64_t const zabspos = zreqvec[z]; // .zabspos;

			uint64_t const zrank = input_types_type::circular_suffix_comparator::suffixSearchTryInternal(
				SA.begin(), T.begin(), T.end(), cblocksize,
				blockstart, zabspos%fs, 
				fn, fs
			);
			
			#if 0
			uint64_t const zrankext = input_types_type::circular_suffix_comparator::suffixSearch(
				SA.begin(), cblocksize, 
				blockstart, zabspos%fs, 
				fn, fs
			);
			assert ( zrankext == zrank );
			#endif

			::libmaus::suffixsort::BwtMergeZBlock zblock(zabspos,zrank);
			result.setZBlock(z,zblock);
		}
		// std::cerr << "done, time " << sufsertc.getElapsedSeconds() << std::endl;

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB18] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// compute BWT		
		::libmaus::bitio::BitVector::unique_ptr_type pGT(new ::libmaus::bitio::BitVector(cblocksize+1));
		::libmaus::bitio::BitVector & GT = *pGT;

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB19] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		bool gtflag = false;
		uint64_t const outcnt = out-SA.begin();
		uint64_t r0 = outcnt;
		
		// construct modified bwt
		for ( saidx_t * in = SA.begin(); in != out; ++in )
		{
			saidx_t const saval = *in;
		
			GT [ saval ] = gtflag;
		
			if ( saval )
			{
				*in = T[saval-1];
			}
			else
			{
				*in = bwtterm;
				// update gt flag
				gtflag = true;
				// set rank of position 0
				r0 = in-SA.begin();
			}
		}

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB20] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		// deallocate text
		PT.reset();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB21] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		GT [ cblocksize ] = gtlast;

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB22] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB23] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		// save gt array
		#if 0
		::libmaus::huffman::HuffmanEncoderFileStd GTHEF(tmpfilenames.getGT());
		for ( int64_t i = cblocksize; i > 0; --i )
			GTHEF.writeBit(GT[i]);
		GTHEF.flush();
		#endif
		
		for ( uint64_t j = 0; j < tmpfilenames.getGT().size(); ++j )
		{
			uint64_t const gtpartsize = (cblocksize + tmpfilenames.getGT().size() - 1)/tmpfilenames.getGT().size();
			uint64_t const low = std::min(j * gtpartsize,cblocksize);
			uint64_t const high = std::min(low + gtpartsize,cblocksize);
			libmaus::bitio::BitVectorOutput BVO(tmpfilenames.getGT()[j]);
			for ( int64_t i = cblocksize-low; i > static_cast<int64_t>(cblocksize-high); --i )
				BVO.writeBit(GT[i]);
			BVO.flush();
		}
		
		#if 0
		libmaus::bitio::BitVectorOutput BVO(tmpfilenames.getGT());
		for ( int64_t i = cblocksize; i > 0; --i )
			BVO.writeBit(GT[i]);
		BVO.flush();
		#endif
		
		pGT.reset();

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB24] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif
		
		uint64_t const targetbwtfilesize = (outcnt + tmpfilenames.getBWT().size() - 1) / tmpfilenames.getBWT().size();

		for ( uint64_t b = 0; b < tmpfilenames.getBWT().size(); ++b )
		{
			uint64_t const low  = std::min(  b * targetbwtfilesize, outcnt);
			uint64_t const high = std::min(low + targetbwtfilesize, outcnt);
			
			#if defined(HUFRL)
			::libmaus::huffman::RLEncoderStd bwtenc(tmpfilenames.getBWT()[b],albits                   ,high-low,rlencoderblocksize);
			#else
			::libmaus::gamma::GammaRLEncoder bwtenc(tmpfilenames.getBWT()[b],albits/* alphabet bits */,high-low,rlencoderblocksize);
			#endif

			if ( low <= r0 && r0 < high )
			{
				for ( uint64_t i = low; i < r0; ++i )
					bwtenc.encode(SA[i]);
				bwtenc.encode(presym);
				for ( uint64_t i = r0+1; i < high; ++i )
					bwtenc.encode(SA[i]);
			}
			else
			{
				for ( uint64_t i = low; i < high; ++i )
					bwtenc.encode(SA[i]);
			}

			#if 0
			// run-length coding for bwt
			for ( uint64_t i = 0; i < r0; ++i )
				bwtenc.encode(SA[i]);
			bwtenc.encode(presym);
			for ( uint64_t i = r0+1; i < outcnt; ++i )
				bwtenc.encode(SA[i]);
			#endif
			
			bwtenc.flush();
			
			#if defined(FERAMANZGEN_DEBUG)
			gcerrlock.lock();
			std::cerr << "[D] generated " << tmpfilenames.getBWT()[b] << " with size " << high-low << std::endl;
			gcerrlock.unlock();
			#endif
		}

		#if defined(FERAMANZGEN_DEBUG)
		gcerrlock.lock();
		std::cerr << "[SB25] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		gcerrlock.unlock();
		#endif

		if ( computeTermSymbolHwt )
		{
			if ( input_types_type::utf8Wavelet() )
			{
				std::string utftmp = tmpfilenames.getHWT() + ".utf8tmp";
				libmaus::util::TempFileRemovalContainer::addTempFile(utftmp);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB26] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
				
				// std::cerr << "writing " << utftmp << std::endl;
				
				::libmaus::util::CountPutObject CPO;
				for ( uint64_t i = 0; i < outcnt; ++i )
					::libmaus::util::UTF8::encodeUTF8(SA[i],CPO);
				uint64_t const ucnt = CPO.c;

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB27] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				libmaus::aio::CheckedOutputStream::unique_ptr_type utfCOS(new libmaus::aio::CheckedOutputStream(utftmp));
				for ( uint64_t i = 0; i < outcnt; ++i )
					::libmaus::util::UTF8::encodeUTF8(SA[i],*utfCOS);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB28] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				utfCOS->flush();
				utfCOS->close();
				utfCOS.reset();

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB29] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				#if 0
				::libmaus::autoarray::AutoArray<uint8_t> UT(ucnt,false);
				::libmaus::util::PutObject<uint8_t *> P(UT.begin());
				for ( uint64_t i = 0; i < outcnt; ++i )
					::libmaus::util::UTF8::encodeUTF8(SA[i],P);
				#endif
					
				SA.release();

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB30] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
				
				::libmaus::autoarray::AutoArray<uint8_t> UT(ucnt,false);
				libmaus::aio::CheckedInputStream::unique_ptr_type utfCIS(new libmaus::aio::CheckedInputStream(utftmp));
				utfCIS->read(reinterpret_cast<char *>(UT.begin()),ucnt,64*1024);
				utfCIS->close();
				utfCIS.reset();

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB31] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				::libmaus::huffman::HuffmanTree::unique_ptr_type uhnode = loadCompactHuffmanTree(huftreefilename);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB32] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
				
				std::string const tmpfileprefix = tmpfilenamedir + "_wt";
				::libmaus::wavelet::Utf8ToImpCompactHuffmanWaveletTree::constructWaveletTree<true>(
					UT,tmpfilenames.getHWT(),tmpfileprefix,uhnode.get(),
					1/* num threads */
				);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB33] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
				
				remove(utftmp.c_str());

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB34] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
			}
			else
			{

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB35] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				::libmaus::huffman::HuffmanTree::unique_ptr_type uhnode = loadCompactHuffmanTree(huftreefilename);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB36] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
				
				// construct huffman shaped wavelet tree
				libmaus::util::FileTempFileContainer FTFC(tmpgen);

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB37] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				::libmaus::wavelet::ImpExternalWaveletGeneratorCompactHuffman IEWGH(*uhnode,FTFC);

				for ( uint64_t i = 0; i < r0; ++i )
					IEWGH.putSymbol(SA[i]);
				IEWGH.putSymbol(bwtterm);
				for ( uint64_t i = r0+1; i < outcnt; ++i )
					IEWGH.putSymbol(SA[i]);		

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB38] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif

				// create final stream for huffman coded wavelet tree
				::libmaus::aio::CheckedOutputStream HCOS(tmpfilenames.getHWT());
				IEWGH.createFinalStream(HCOS);
				HCOS.flush();

				#if defined(FERAMANZGEN_DEBUG)
				gcerrlock.lock();
				std::cerr << "[SB39] " << libmaus::util::MemUsage() << "," << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
				gcerrlock.unlock();
				#endif
			}
		}
		else
		{
			libmaus::aio::CheckedOutputStream hwtReqCOS(tmpfilenames.getHWTReq());
			RlToHwtTermRequest::serialise(hwtReqCOS,
				tmpfilenames.getBWT(),
				tmpfilenames.getHWT(),
				tmpfilenamedir,
				huftreefilename,
				bwtterm,
				r0,
				input_types_type::utf8Wavelet()
			);
			hwtReqCOS.flush();
			hwtReqCOS.close();
		}
		
		// glock.unlock();
		
		return result;
	}

	
	std::string dispatch() const
	{
		::libmaus::suffixsort::BwtMergeBlockSortResult result;
		
		switch (inputtype)
		{
			case bwt_merge_input_type_byte:
				result = sortBlock<libmaus::suffixsort::ByteInputTypes>();
				break;
			case bwt_merge_input_type_compact:
				result = sortBlock<libmaus::suffixsort::CompactInputTypes>();
				break;
			case bwt_merge_input_type_pac:
				result = sortBlock<libmaus::suffixsort::PacInputTypes>();
				break;
			case bwt_merge_input_type_pac_term:
				result = sortBlock<libmaus::suffixsort::PacTermInputTypes>();
				break;
			case bwt_merge_input_type_lz4:
				result = sortBlock<libmaus::suffixsort::Lz4InputTypes>();
				break;
			case bwt_merge_input_type_utf8:
				result = sortBlock<libmaus::suffixsort::Utf8InputTypes>();
				break;
			default:
			{
				::libmaus::exception::LibMausException ex;
				ex.getStream() << "Number " << inputtype << " is not a valid input type designator." << std::endl;
				ex.finish();
				throw ex;				
			}
		}
		return result.serialise();
	}
};


std::ostream & operator<<(std::ostream & out, BwtMergeBlockSortRequest const & o)
{
	out << "BwtMergeBlockSortRequest(";
	out << o.inputtype;
	out << ",";
	out << o.fn;
	out << ",";
	out << o.fs;
	out << ",";
	out << o.rlencoderblocksize;
	out << ",";
	out << o.isasamplingrate;
	out << ",";
	out << o.blockstart;
	out << ",";
	out << o.cblocksize;
	out << ",";
	out << "{";
	for ( uint64_t i = 0; i < o.zreqvec.size(); ++i )
	{
		out << o.zreqvec[i].getZAbsPos();
		if ( i+1 < o.zreqvec.size() )
			out << ";";
	}
	out << "}";
	out << ",";
	out << o.computeTermSymbolHwt;
	out << ")";
	
	return out;
}

// forward declaration
struct MergeStrategyMergeGapRequest;

struct MergeStrategyMergeGapRequestQueryObject
{
	uint64_t p;
	uint64_t r;
	MergeStrategyMergeGapRequest * o;
	
	MergeStrategyMergeGapRequestQueryObject()
	: p(0), r(0), o(0)
	{
	
	}
	
	MergeStrategyMergeGapRequestQueryObject(uint64_t const rp)
	: p(rp), r(0), o(0)
	{
	
	}

	MergeStrategyMergeGapRequestQueryObject(uint64_t const rp, uint64_t const rr, MergeStrategyMergeGapRequest * ro)
	: p(rp), r(rr), o(ro)
	{
	
	}
	
	bool operator<(MergeStrategyMergeGapRequestQueryObject const & o) const
	{
		return p < o.p;
	}
};

struct MergeStrategyBlock
{
	typedef MergeStrategyBlock this_type;
	typedef libmaus::util::shared_ptr<MergeStrategyBlock>::type shared_ptr_type;

	//! low symbol offset in file
	uint64_t low;
	//! high symbol offset in file
	uint64_t high;
	//! length of the source block (on disk) in bits
	uint64_t sourcelengthbits;
	//! length of the source block in bytes (for byte oriented suffix sorting)
	uint64_t sourcelengthbytes;
	//! length of sequence in Huffman code
	uint64_t codedlength;
	//! length of index for random access in text
	uint64_t sourcetextindexbits;
	//! node id
	uint64_t nodeid;
	//! node depth
	uint64_t nodedepth;
	//! sorting result
	::libmaus::suffixsort::BwtMergeBlockSortResult sortresult;
	//! parent node
	MergeStrategyBlock * parent;

	MergeStrategyBlock()
	: low(0), high(0), sourcelengthbits(0), sourcelengthbytes(0), codedlength(0), sourcetextindexbits(0), nodeid(0), nodedepth(0), parent(0) {}
	MergeStrategyBlock(uint64_t const rlow, uint64_t const rhigh, uint64_t const rsourcelengthbits, uint64_t const rsourcelengthbytes, uint64_t const rsourcetextindexbits)
	: low(rlow), high(rhigh), sourcelengthbits(rsourcelengthbits), sourcelengthbytes(rsourcelengthbytes), codedlength(0), sourcetextindexbits(rsourcetextindexbits), nodeid(0), nodedepth(0), parent(0) {}
	
	virtual ~MergeStrategyBlock() {}
	virtual std::ostream & print(std::ostream & out, uint64_t const indent) const = 0;
	
	/**
	 * get estimate for space used by Huffman shaped wavelet tree in bits
	 **/
	uint64_t getIHWTSpaceBits() const
	{
		return (codedlength * 4+2)/3;
	}

	/**
	 * get estimate for space used by Huffman shaped wavelet tree in bytes
	 **/
	uint64_t getIHWTSpaceBytes() const
	{
		return (getIHWTSpaceBits()+7)/8;
	}
	
	/**
	 * get estimate for number of bytes required to merge into this block using an internal memory gap array
	 **/
	uint64_t getMergeSpaceInternalGap() const
	{
		return 
			getIHWTSpaceBytes() + (high-low) * sizeof(uint32_t);
	}
	/**
	 * get estimate for number of bytes required to merge into this block using an internal memory small gap array (1 byte + overflow)
	 **/
	uint64_t getMergeSpaceInternalSmallGap() const
	{
		return 
			getIHWTSpaceBytes() + (high-low) * sizeof(uint8_t);
	}

	/**
	 * get estimate for number of bytes required to merge into this block using an external memory gap array
	 **/
	uint64_t getMergeSpaceExternalGap(uint64_t const threads, uint64_t wordsperthread) const
	{
		return 
			getIHWTSpaceBytes() + threads * wordsperthread * sizeof(uint64_t);
	}
	
	/**
	 * @return space required in bytes for sorting block in internal memory using divsufsort
	 **/
	uint64_t directSortSpace() const
	{
		if ( (high-low) < (1ull << 31) )
		{
			return 
				(high-low)*(sizeof(uint32_t)) + // suffix array
				(high-low+7)/8 + // GT bit vector
				(sourcetextindexbits+7)/8 +
				sourcelengthbytes;
		}
		else
		{
			return (high-low)*(sizeof(uint64_t)) + // suffix array
				(high-low+7)/8 + // GT bit vector
				(sourcetextindexbits+7)/8 +
				sourcelengthbytes;
		}
	}

	std::ostream & printBase(std::ostream & out) const
	{
		out 
			<< "id=" << nodeid
			<< ",d=" << nodedepth
			<< ",[" << low << "," << high << "),sourcelengthbits=" << sourcelengthbits << ",sourcelengthbytes=" << sourcelengthbytes << ",codedlength=" << codedlength << ",directsortspace=" << directSortSpace() << ",sourcetextindexbits=" << sourcetextindexbits;
		return out;
	}
	
	virtual uint64_t fillNodeId(uint64_t i) = 0;
	virtual void fillNodeDepth(uint64_t const i) = 0;
	// register the given query positions in this leaf or all leafs under this inner node
	virtual void registerQueryPositions(std::vector<uint64_t> const & V) = 0;
	// fill query positions for t thread for each gap request in the tree
	virtual void fillQueryPositions(uint64_t const t) = 0;

	// fill query objects for inner nodes given finished leaf node information
	virtual void fillQueryObjects(libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV) = 0;
	// fill zblock ranks from finished base blocks for t threads
	virtual void fillGapRequestObjects(uint64_t const t) = 0;
	
	virtual bool isLeaf() const = 0;
	virtual void setParent(MergeStrategyBlock * rparent) = 0;
	virtual bool childFinished() = 0;
};

struct MergeStrategyBaseBlock : public MergeStrategyBlock
{
	BwtMergeBlockSortRequest sortreq;
	std::vector<uint64_t> querypos;

	MergeStrategyBaseBlock() : MergeStrategyBlock() {}

	MergeStrategyBaseBlock(
		uint64_t const rlow, uint64_t const rhigh, 
		::libmaus::huffman::HuffmanTree::EncodeTable & EC,
		std::map<int64_t,uint64_t> const & blockhist,
		uint64_t const rsourcelengthbits,
		uint64_t const rsourcelengthbytes,
		uint64_t const rsourcetextindexbits
	)
	: MergeStrategyBlock(rlow,rhigh,rsourcelengthbits,rsourcelengthbytes,rsourcetextindexbits)
	{
		for (
			std::map<int64_t,uint64_t>::const_iterator ita = blockhist.begin();
			ita != blockhist.end();
			++ita )
		{
			MergeStrategyBlock::codedlength += (EC).getCodeLength(ita->first) * ita->second;
		}
	}
	
	std::ostream & print(std::ostream & out, uint64_t const indent) const
	{
		out << "[V]" << std::string(indent,' ') << "MergeStrategyBaseBlock(";
		printBase(out);
		out << ")" << std::endl;
		out << "[V]" << std::string(indent+1,' ') << "qp={";
		for ( uint64_t i = 0; i < querypos.size(); ++i )
		{
			out << querypos[i];
			if ( i+1 < querypos.size() )
				out << ";";
		}
		out << "}" << std::endl;
		out << "[V]" << std::string(indent+1,' ') << "req=";
		out << sortreq;
		out << std::endl;
		return out;
	}
	
	static MergeStrategyBlock::shared_ptr_type construct(
		uint64_t const rlow, uint64_t const rhigh, 
		::libmaus::huffman::HuffmanTree::EncodeTable & EC,
		std::map<int64_t,uint64_t> const & blockhist,
		uint64_t const rsourcelengthbits,
		uint64_t const rsourcelengthbytes,
		uint64_t const rsourcetextindexbits
	)
	{
		return MergeStrategyBlock::shared_ptr_type(new MergeStrategyBaseBlock(rlow,rhigh,EC,blockhist,rsourcelengthbits,rsourcelengthbytes,rsourcetextindexbits));
	}

	uint64_t fillNodeId(uint64_t i)
	{
		MergeStrategyBlock::nodeid = i++;
		return i;
	}

	void fillNodeDepth(uint64_t const i)
	{
		MergeStrategyBlock::nodedepth = i;
	}

	virtual void registerQueryPositions(std::vector<uint64_t> const & V)
	{
		std::copy(V.begin(),V.end(),std::back_insert_iterator< std::vector<uint64_t> >(querypos));
	}
	
	void fillQueryPositions(uint64_t const /* t */)
	{
		::libmaus::suffixsort::BwtMergeZBlockRequestVector zreqvec;
		zreqvec.resize(querypos.size());
		for ( uint64_t i = 0; i < querypos.size(); ++i )
			zreqvec[i] = ::libmaus::suffixsort::BwtMergeZBlockRequest(querypos[i]);
		sortreq.zreqvec = zreqvec;
	}

	void fillQueryObjects(libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV)
	{
		libmaus::autoarray::AutoArray < ::libmaus::suffixsort::BwtMergeZBlock > const & Z = sortresult.getZBlocks();
		
		for ( uint64_t i = 0; i < Z.size(); ++i )
		{
			uint64_t const p = Z[i].getZAbsPos();
			uint64_t const r = Z[i].getZRank();
			
			// search objects with matching position
			typedef libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject>::iterator it;
			std::pair<it,it> range = std::equal_range(VV.begin(),VV.end(),MergeStrategyMergeGapRequestQueryObject(p));
			
			// add rank r to each object for position p
			for ( it ita = range.first; ita != range.second; ++ita )
				ita->r += r;
		}
	}
	void fillGapRequestObjects(uint64_t const)
	{
	
	}
	bool isLeaf() const
	{
		return true;
	}
	void setParent(MergeStrategyBlock * rparent)
	{
		parent = rparent;
	}
	bool childFinished()
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "childFinished called on an object of struct MergeStrategyBaseBlock" << std::endl;
		se.finish();
		throw se;
	}
};

/**
 * sorting thread for base blocks
 **/
struct BaseBlockSortThread : public libmaus::parallel::PosixThread
{
	typedef BaseBlockSortThread this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	//! thread id
	uint64_t tid;

	/**
	 * semaphore. delivers a message whenever there is sufficient
	 * free space to process the next element
	 **/
	libmaus::parallel::PosixSemaphore & P;
	//! vector of blocks to be processed
	std::vector < MergeStrategyBlock::shared_ptr_type > & V;

	//! next package to be processed
	volatile uint64_t & next;
	//! amount of free memory
	volatile uint64_t & freemem;
	//! number of finished threads
	volatile uint64_t & finished;
	//! lock for the above
	libmaus::parallel::PosixMutex & freememlock;
	//! inner node queue
	std::deque<MergeStrategyBlock *> & itodo;
	
	BaseBlockSortThread(
		uint64_t rtid,
		libmaus::parallel::PosixSemaphore & rP,
		std::vector < MergeStrategyBlock::shared_ptr_type > & rV,
		uint64_t & rnext,
		uint64_t & rfreemem,
		uint64_t & rfinished,
		libmaus::parallel::PosixMutex & rfreememlock,
		std::deque<MergeStrategyBlock *> & ritodo
	) : tid(rtid), P(rP), V(rV), next(rnext), freemem(rfreemem), finished(rfinished), freememlock(rfreememlock),
	    itodo(ritodo)
	{
	
	}
	
	void * run()
	{
		bool running = true;
		
		while ( running )
		{
			// wait until sufficient memory is free
			P.wait();
						
			// get package id
			uint64_t pack = 0;

			{
				// get lock
				libmaus::parallel::ScopePosixMutex scopelock(freememlock);

				if ( next < V.size() )
				{
					assert ( V[next]->directSortSpace() <= freemem );
					freemem -= V[next]->directSortSpace();
					pack = next++;				
					
					if ( next >= V.size() || freemem >= V[next]->directSortSpace() )
						P.post();			
				}
				else
				{
					running = false;
					P.post();
				}				
			}

			if ( running )
			{
				try
				{
					// perform sorting
					MergeStrategyBaseBlock * block = dynamic_cast<MergeStrategyBaseBlock *>(V[pack].get());					
					block->sortresult = ::libmaus::suffixsort::BwtMergeBlockSortResult::load(block->sortreq.dispatch());
				}
				catch(std::exception const & ex)
				{
					libmaus::parallel::ScopePosixMutex scopelock(freememlock);
					std::cerr << tid << " failed " << pack << " " << ex.what() << std::endl;
				}

				{
					// get lock
					libmaus::parallel::ScopePosixMutex scopelock(freememlock);
					
					std::cerr << "[V] [" << tid << "] sorted block " << pack << std::endl;
				
					// "free" memory
					freemem += V[pack]->directSortSpace();

					if ( V[pack]->parent )
					{
						bool const pfinished = V[pack]->parent->childFinished();
					
						if ( pfinished )
						{
							itodo.push_back(V[pack]->parent);
						}
					}
					
					// post if there is room for another active sorting thread
					if ( next == V.size() || freemem >= V[next]->directSortSpace() )
						P.post();
				}
			}
		}
		
		{
		libmaus::parallel::ScopePosixMutex scopelock(freememlock);
		finished++;
		}
		
		// quit
		return 0;
	}
};

/**
 * a set of thread for block sorting
 **/
struct BaseBlockSorting
{
	typedef BaseBlockSorting this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;

	std::vector < MergeStrategyBlock::shared_ptr_type > & V;
	
	libmaus::parallel::PosixSemaphore P;
	uint64_t next;
	uint64_t freemem;
	uint64_t finished;
	libmaus::parallel::PosixMutex freememlock;
	//! inner node queue
	std::deque<MergeStrategyBlock *> & itodo;

	libmaus::autoarray::AutoArray<BaseBlockSortThread::unique_ptr_type> threads;

	BaseBlockSorting(
		std::vector < MergeStrategyBlock::shared_ptr_type > & rV,
		uint64_t const rfreemem,
		uint64_t const numthreads,
		//! inner node queue
		std::deque<MergeStrategyBlock *> & ritodo
	)
	: V(rV), P(), next(0), freemem(rfreemem), finished(0), freememlock(), itodo(ritodo), threads(numthreads)
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i]->directSortSpace() > freemem )
			{
				libmaus::exception::LibMausException se;
				se.getStream() << "Memory provided is " << freemem << " but " 
					<< V[i]->directSortSpace() << " are required for sorting block " << i << std::endl;
				se.finish();
				throw se;
			}
	
		for ( uint64_t i = 0; i < numthreads; ++i )
			threads[i] =
				UNIQUE_PTR_MOVE(
					BaseBlockSortThread::unique_ptr_type(
						new BaseBlockSortThread(
							i,P,V,next,freemem,finished,freememlock,itodo
						)
					)
				);
	}
	
	void start(uint64_t const stacksize)
	{
		for ( uint64_t i = 0; i < threads.size(); ++i )
			threads[i]->startStack(stacksize);
			
		P.post();
	}

	void start()
	{
		for ( uint64_t i = 0; i < threads.size(); ++i )
			threads[i]->start();
			
		P.post();
	}

	void join()
	{
		try
		{
			for ( uint64_t i = 0; i < threads.size(); ++i )
				threads[i]->join();	
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
		}
	}
};

struct MergeStrategyMergeGapRequest
{
	typedef MergeStrategyMergeGapRequest this_type;
	typedef libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus::util::shared_ptr<this_type>::type shared_ptr_type;

	std::vector<MergeStrategyBlock::shared_ptr_type> const * pchildren;
	uint64_t into;
	std::vector < ::libmaus::suffixsort::BwtMergeZBlock > zblocks;
	
	MergeStrategyMergeGapRequest() : pchildren(0), into(0), zblocks() {}
	MergeStrategyMergeGapRequest(
		std::vector<MergeStrategyBlock::shared_ptr_type> const * rpchildren, 
		uint64_t const rinto)
	: pchildren(rpchildren), into(rinto), zblocks()
	{
	}
	
	std::vector<uint64_t> getQueryPositions(uint64_t const t) const
	{
		assert ( pchildren );
		assert ( pchildren->size() );
		assert ( into != pchildren->size()-1 );
		assert ( t );
		
		std::vector<MergeStrategyBlock::shared_ptr_type> const & children = *pchildren;
		std::vector<uint64_t> Q;
		
		// absolute low pos of inserted blocks
		uint64_t const abslow = children[into+1]->low;
		// absolute high pos of inserted blocks
		uint64_t const abshigh = children[children.size()-1]->high;
		// acc size of inserted blocks
		uint64_t const abssize = abshigh-abslow;
		assert ( abssize );
		// distance between query points
		uint64_t const absdif = (abssize + t - 1)/t;
		
		// count query positions
		uint64_t n = 0;
		for ( uint64_t i = 0; i < t; ++i )
			if ( i * absdif < abssize )
				++n;
		
		// store query points (falling positions)
		Q.resize(n);
		for ( uint64_t i = 0; i < t; ++i )
			if ( i * absdif < abssize )
				Q[i] = abshigh - i * absdif;

		// check query points
		for ( uint64_t i = 0; i < Q.size(); ++i )
			assert ( Q[i] > abslow && Q[i] <= abshigh );
		
		return Q;
	}
	
	void getQueryPositionObjects(
		std::vector < MergeStrategyMergeGapRequestQueryObject > & VV,
		uint64_t const t
	)
	{
		std::vector<uint64_t> const Q = getQueryPositions(t);
		
		for ( uint64_t i = 0; i < Q.size(); ++i )
			VV.push_back(MergeStrategyMergeGapRequestQueryObject(Q[i],0,this));
	}

	libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> getQueryPositionObjects(uint64_t const t)
	{
		std::vector<uint64_t> const Q = getQueryPositions(t);
		libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> VV(Q.size());
		
		for ( uint64_t i = 0; i < Q.size(); ++i )
			VV[i] = MergeStrategyMergeGapRequestQueryObject(Q[i],0,this);
			
		return VV;
	}
};

std::ostream & operator<<(std::ostream & out, MergeStrategyMergeGapRequest const & G)
{
	out << "MergeStrategyMergeGapRequest(into=" << G.into << ",";
	for ( uint64_t i = 0; i < G.zblocks.size(); ++i )
		out << G.zblocks[i] << ";";
	out << ")";
	
	return out;
}

struct MergeStrategyMergeBlock : public MergeStrategyBlock
{
	std::vector<MergeStrategyBlock::shared_ptr_type> children;
	std::vector<MergeStrategyMergeGapRequest::shared_ptr_type> gaprequests;
	uint64_t unfinishedChildren;

	MergeStrategyMergeBlock() 
	: MergeStrategyBlock()
	{
	}
	
	void releaseChildren()
	{
		gaprequests.resize(0);
		children.resize(0);
	}

	virtual std::ostream & print(std::ostream & out, uint64_t const indent) const = 0;

	//void fillQueryObjects(std::vector<MergeStrategyMergeGapRequestQueryObject> & VV)
	void fillQueryObjects(libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV)
	{
		for ( uint64_t i = 0; i < children.size(); ++i )
			children[i]->fillQueryObjects(VV);
	}

	void fillGapRequestObjects(uint64_t const t)
	{
		for ( uint64_t i = 0; i < gaprequests.size(); ++i )
		{
			// std::vector < MergeStrategyMergeGapRequestQueryObject > VV;
			// get gap query objects
			libmaus::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> VV = gaprequests[i]->getQueryPositionObjects(/* VV, */t);
			std::sort(VV.begin(),VV.end());
			// fill the objects (add up ranks for each query position)
			children[gaprequests[i]->into]->fillQueryObjects(VV);
			// size of vector
			uint64_t const vn = VV.size();

			// push z blocks back to front
			for ( uint64_t i = 0; i < VV.size(); ++i )
			{
				uint64_t const ii = vn-i-1;
				VV[i].o->zblocks.push_back(::libmaus::suffixsort::BwtMergeZBlock(VV[ii].p,VV[ii].r));
			}
		}
				
		for ( uint64_t i = 0; i < children.size(); ++i )
			children[i]->fillGapRequestObjects(t);			
	}
	
	std::ostream & print(std::ostream & out, uint64_t const indent, std::string const & name) const
	{
		out << "[V]" << std::string(indent,' ')<< name << "(";
		printBase(out);
		out << ")" << std::endl;

		for ( uint64_t i = 0; i < gaprequests.size(); ++i )
			out << "[V]" << std::string(indent+1,' ') << (*(gaprequests[i])) << std::endl;

		for ( uint64_t i = 0; i < children.size(); ++i )
			children[i]->print(out,indent+1);

		return out;
	}
	
	static void construct(MergeStrategyMergeBlock * pobj, std::vector<MergeStrategyBlock::shared_ptr_type> const children)
	{
		pobj->children = children;

		uint64_t low = children.size() ? children[0]->low : 0;
		uint64_t high = children.size() ? children[0]->high : 0;
		uint64_t sourcelengthbits = 0;
		uint64_t sourcelengthbytes = 0;
		uint64_t codedlength = 0;
		uint64_t sourcetextindexbits = 0;
				
		for ( uint64_t i = 0; i < children.size(); ++i )
		{
			low = std::min(low,children[i]->low);
			high = std::max(high,children[i]->high);
			sourcelengthbits += children[i]->sourcelengthbits;
			sourcelengthbytes += children[i]->sourcelengthbytes;
			codedlength += children[i]->codedlength;
			sourcetextindexbits += children[i]->sourcetextindexbits;
		}
		
		pobj->low = low;
		pobj->high = high;
		pobj->sourcelengthbits = sourcelengthbits;
		pobj->sourcelengthbytes = sourcelengthbytes;
		pobj->codedlength = codedlength;
		pobj->sourcetextindexbits = sourcetextindexbits;
		pobj->unfinishedChildren = children.size();
		
		for ( uint64_t i = 0; i+1 < children.size(); ++i )
		{
			pobj->gaprequests.push_back(
				MergeStrategyMergeGapRequest::shared_ptr_type(
					new MergeStrategyMergeGapRequest(&(pobj->children),i))
			);
		}
	}

	uint64_t fillNodeId(uint64_t i)
	{
		MergeStrategyBlock::nodeid = i++;
		
		for ( uint64_t j = 0; j < children.size(); ++j )
			i = children[j]->fillNodeId(i);
			
		return i;
	}

	void fillNodeDepth(uint64_t const i)
	{
		MergeStrategyBlock::nodedepth = i;

		for ( uint64_t j = 0; j < children.size(); ++j )
			children[j]->fillNodeDepth(i+1);
	}

	virtual void registerQueryPositions(std::vector<uint64_t> const & V)
	{
		for ( uint64_t i = 0; i < children.size(); ++ i )
			children[i]->registerQueryPositions(V);
	}

	virtual void fillQueryPositions(uint64_t const t)
	{
		for ( uint64_t i = 0; i < gaprequests.size(); ++i )
		{
			// get query position from gap request
			std::vector<uint64_t> const Q = gaprequests[i]->getQueryPositions(t);
			// register these positions in the leafs
			children[i]->registerQueryPositions(Q);
		}
		// recursive call for children
		for ( uint64_t i = 0; i < children.size(); ++i )
			children[i]->fillQueryPositions(t);
	}
	virtual bool isLeaf() const
	{
		return false;
	}

	void setParent(MergeStrategyBlock * rparent)
	{
		parent = rparent;
		
		for ( uint64_t i = 0; i < children.size(); ++i )
			children[i]->setParent(this);
	}
	bool childFinished()
	{
		return --unfinishedChildren == 0;
	}
};

struct MergeStrategyMergeInternalBlock : public MergeStrategyMergeBlock
{
	// std::vector<MergeStrategyBlock::shared_ptr_type> children;

	MergeStrategyMergeInternalBlock() 
	: MergeStrategyMergeBlock()
	{
	}

	std::ostream & print(std::ostream & out, uint64_t const indent) const
	{
		return MergeStrategyMergeBlock::print(out,indent,"MergeStrategyMergeInternalBlock");
	}
	
	static MergeStrategyBlock::shared_ptr_type construct(
		std::vector<MergeStrategyBlock::shared_ptr_type> const children
	)
	{
		MergeStrategyMergeInternalBlock * pobj = new MergeStrategyMergeInternalBlock();
		MergeStrategyBlock::shared_ptr_type sobj(pobj);
		MergeStrategyMergeBlock::construct(pobj,children);
		return sobj;
	}
};

struct MergeStrategyMergeInternalSmallBlock : public MergeStrategyMergeBlock
{
	// std::vector<MergeStrategyBlock::shared_ptr_type> children;

	MergeStrategyMergeInternalSmallBlock() 
	: MergeStrategyMergeBlock()
	{
	}

	std::ostream & print(std::ostream & out, uint64_t const indent) const
	{
		return MergeStrategyMergeBlock::print(out,indent,"MergeStrategyMergeInternalSmallBlock");
	}
	
	static MergeStrategyBlock::shared_ptr_type construct(
		std::vector<MergeStrategyBlock::shared_ptr_type> const children
	)
	{
		MergeStrategyMergeInternalSmallBlock * pobj = new MergeStrategyMergeInternalSmallBlock();
		MergeStrategyBlock::shared_ptr_type sobj(pobj);
		MergeStrategyMergeBlock::construct(pobj,children);
		return sobj;
	}
};

struct MergeStrategyMergeExternalBlock : public MergeStrategyMergeBlock
{
	// std::vector<MergeStrategyBlock::shared_ptr_type> children;

	MergeStrategyMergeExternalBlock() 
	: MergeStrategyMergeBlock()
	{
	}

	std::ostream & print(std::ostream & out, uint64_t const indent) const
	{
		return MergeStrategyMergeBlock::print(out,indent,"MergeStrategyMergeExternalBlock");
	}
	
	static MergeStrategyBlock::shared_ptr_type construct(
		std::vector<MergeStrategyBlock::shared_ptr_type> const children
	)
	{
		MergeStrategyMergeExternalBlock * pobj = new MergeStrategyMergeExternalBlock();
		MergeStrategyBlock::shared_ptr_type sobj(pobj);
		MergeStrategyMergeBlock::construct(pobj,children);
		return sobj;
	}
};

std::ostream & operator<<(std::ostream & out, MergeStrategyBlock const & MSB)
{
	return MSB.print(out,0);
}

MergeStrategyBlock::shared_ptr_type constructMergeTreeRec(
	std::vector<MergeStrategyBlock::shared_ptr_type> nodes,
	uint64_t const mem,
	uint64_t const numthreads,
	uint64_t const exwordsperthread
)
{
	if ( !nodes.size() )
		return MergeStrategyBlock::shared_ptr_type();
	else if ( nodes.size() == 1 )
		return nodes.at(0);
	
	/**
	 * see if we can do a pairwise merge
	 **/
	std::vector<MergeStrategyBlock::shared_ptr_type> outnodes;
	uint64_t const outnodessize = (nodes.size()+1)/2;
	bool pairok = true;
		
	for ( uint64_t i = 0; pairok && i < outnodessize; ++i )
	{
		uint64_t const low = 2*i;
		uint64_t const high = std::min(low+2,static_cast<uint64_t>(nodes.size()));
		
		MergeStrategyBlock::shared_ptr_type node;
		
		if ( high-low < 2 )
		{
			node = nodes[low];
		}
		if ( (!node) && (nodes[low]->getMergeSpaceInternalGap() <= mem) )
		{
			node = MergeStrategyMergeInternalBlock::construct(std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin()+low,nodes.begin()+high));
		}
		if ( (!node) && (nodes[low]->getMergeSpaceInternalSmallGap() <= mem) )
		{
			node = MergeStrategyMergeInternalSmallBlock::construct(std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin()+low,nodes.begin()+high));
		}
		if ( (!node) && (nodes[low]->getMergeSpaceExternalGap(numthreads,exwordsperthread) <= mem) )
		{
			node = MergeStrategyMergeExternalBlock::construct(std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin()+low,nodes.begin()+high));
		}
		
		if ( node )
			outnodes.push_back(node);
		else
			pairok = false;
	}
	
	MergeStrategyBlock::shared_ptr_type node;
	
	if ( pairok )
		node = constructMergeTreeRec(outnodes,mem,numthreads,exwordsperthread);

	/*
	 * if pairwise merge does not work, try multi-way merge without recursion
	 */ 
	if ( ! node )
	{
		bool internalok = true;
		for ( uint64_t i = 0; internalok && i+1 < nodes.size(); ++i )
			if ( nodes.at(i)->getMergeSpaceInternalGap() > mem )
				internalok = false;

		if ( internalok )		
			node = MergeStrategyMergeInternalBlock::construct(
				std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin(),nodes.end())
			);
	}

	if ( ! node )
	{
		bool internalsmallok = true;
		for ( uint64_t i = 0; internalsmallok && i+1 < nodes.size(); ++i )
			if ( nodes.at(i)->getMergeSpaceInternalSmallGap() > mem )
				internalsmallok = false;

		if ( internalsmallok )		
			node = MergeStrategyMergeInternalSmallBlock::construct(
				std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin(),nodes.end())
			);
	}
	
	if ( ! node )
	{
		bool externalok = true;
		for ( uint64_t i = 0; externalok && i+1 < nodes.size(); ++i )
			if ( nodes.at(i)->getMergeSpaceExternalGap(numthreads,exwordsperthread) > mem )
				externalok = false;
		
		if ( externalok )
			node = MergeStrategyMergeExternalBlock::construct(
				std::vector<MergeStrategyBlock::shared_ptr_type>(nodes.begin(),nodes.end())
			);
	}
				
	return node;
}

MergeStrategyBlock::shared_ptr_type constructMergeTree(
	std::vector<MergeStrategyBlock::shared_ptr_type> nodes,
	uint64_t const mem,
	uint64_t const numthreads,
	uint64_t const exwordsperthread
)
{
	std::cerr << "[V] Number of input leaf nodes for merge tree construction is " << nodes.size() << std::endl;

	MergeStrategyBlock::shared_ptr_type node = constructMergeTreeRec(nodes,mem,numthreads,exwordsperthread);

	if ( ! node )
	{
		libmaus::exception::LibMausException se;
		se.getStream() << "There is no way to merge with these settings." << std::endl;
		se.finish();
		throw se;
	}
	
	node->fillNodeId(0);
	node->fillNodeDepth(0);
	node->fillQueryPositions(numthreads);
	node->setParent(0);
	
	return node;
}


template<typename input_types_type>
struct BwtMergeSort
{
	static ::std::map<int64_t,uint64_t> computeSymbolHistogram(std::string const & fn)
	{
		typedef typename input_types_type::linear_wrapper stream_type;
		typedef typename input_types_type::base_input_stream::char_type char_type;
		typedef typename ::libmaus::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
		uint64_t const fs = stream_type::getFileSize(fn);
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		uint64_t const symsperfrag = (fs + numthreads - 1)/numthreads;
		uint64_t const numfrags = (fs + symsperfrag - 1)/symsperfrag;

		libmaus::util::HistogramSet HS(numfrags,256);
		
		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif
		for ( int64_t t = 0; t < static_cast<int64_t>(numfrags); ++t )
		{
			uint64_t const low = t * symsperfrag;
			uint64_t const high = std::min(low+symsperfrag,fs);
			uint64_t const size = high-low;
			uint64_t todo = size;

			stream_type CIS(fn);
			CIS.seekg(low);

			::libmaus::autoarray::AutoArray<unsigned_char_type> B(16*1024,false);
			libmaus::util::Histogram & H = HS[t];
			
			while ( todo )
			{
				uint64_t const toread = std::min(todo,B.size());
				CIS.read(reinterpret_cast<char_type *>(B.begin()),toread);
				assert ( CIS.gcount() == static_cast<int64_t>(toread) );
				for ( uint64_t i = 0; i < toread; ++i )
					H (B[i]);
				todo -= toread;
			}
		}

		::libmaus::util::Histogram::unique_ptr_type PH(UNIQUE_PTR_MOVE(HS.merge()));
		::libmaus::util::Histogram & H(*PH);

		return H.getByType<int64_t>();
	}

	static ::std::map<int64_t,uint64_t> computeSymbolHistogramPlusTerm(std::string const & fn, int64_t & term)
	{
		::std::map<int64_t,uint64_t> M = computeSymbolHistogram(fn);

		term = M.size() ? (M.rbegin()->first+1) : 0;
		M [ term ] ++;
		
		return M;
	}

	template<typename hef_type>
	static void appendBitVectorToFile(
		std::string const & fn,
		uint64_t const n,
		hef_type & HEF
	)
	{
		typedef ::libmaus::huffman::BitInputBuffer4 sbis_type;			
		::libmaus::aio::CheckedInputStream::unique_ptr_type istr;
		sbis_type::raw_input_ptr_type ript;
		sbis_type::unique_ptr_type SBIS;
		
		istr = UNIQUE_PTR_MOVE(::libmaus::aio::CheckedInputStream::unique_ptr_type(new ::libmaus::aio::CheckedInputStream(fn)));
		ript = UNIQUE_PTR_MOVE(sbis_type::raw_input_ptr_type(new sbis_type::raw_input_type(*istr)));
		SBIS = UNIQUE_PTR_MOVE(sbis_type::unique_ptr_type(new sbis_type(ript,static_cast<uint64_t>(64*1024))));
		
		uint64_t const fullbytes = n / 8;
		for ( uint64_t i = 0; i < fullbytes; ++i )
			HEF.write(SBIS->read(8),8);
		
		for ( uint64_t i = fullbytes*8; i < n; ++i )
		{
			HEF.writeBit(SBIS->readBit());
		}
	}

	static void checkBwtBlockDecode(
		std::pair<uint64_t,uint64_t> const & isai,
		std::pair<uint64_t,uint64_t> const & isapre,
		std::string const & fn,
		uint64_t const fs,
		::libmaus::lf::ImpCompactHuffmanWaveletLF const & IHWT,
		::libmaus::aio::SynchronousGenericOutput<uint64_t> & SGO,
		::libmaus::aio::SynchronousGenericOutput<uint64_t> & ISGO,
		uint64_t const sasamplingrate,
		uint64_t const isasamplingrate,
		int64_t const ibs = -1
	)
	{
		assert ( ::libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(sasamplingrate) == 1 );
		assert ( ::libmaus::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(isasamplingrate) == 1 );
		
		uint64_t r = isai.first;
		uint64_t p = isai.second;
		uint64_t const sasamplingmask = sasamplingrate-1;
		uint64_t const isasamplingmask = isasamplingrate-1;
		// size of block
		uint64_t const bs = ibs >= 0 ? ibs : ((p >  isapre.second) ? (p-isapre.second) : (fs - isapre.second));

		// std::cerr << "bs=" << bs << std::endl;

		typename input_types_type::circular_reverse_wrapper CRWR(fn,p);
		
		if ( ! p )
		{
			for ( uint64_t j = 0; j < bs; ++j )
			{
				if ( !(r & sasamplingmask) )
				{
					SGO.put(r);
					SGO.put(p);
				}
				if ( !(p & isasamplingmask) )
				{
					ISGO.put(p);
					ISGO.put(r);
				}
			
				int64_t syma = CRWR.get();
				int64_t symb = IHWT[r];
				assert ( syma == symb );
				// std::cerr << "(" << syma << "," << symb << ")";
				
				r = IHWT(r);
				
				if ( ! p )
					p = fs;
				--p;
			}
		}
		else
		{
			for ( uint64_t j = 0; j < bs; ++j )
			{
				if ( !(r & sasamplingmask) )
				{
					SGO.put(r);
					SGO.put(p);
				}
				if ( !(p & isasamplingmask) )
				{
					ISGO.put(p);
					ISGO.put(r);
				}
			
				int64_t syma = CRWR.get();
				int64_t symb = IHWT[r];
				assert ( syma == symb );
				// std::cerr << "(" << syma << "," << symb << ")";
				
				r = IHWT(r);
				
				--p;
			}
		}
		
		assert ( r == isapre.first );
	}

	static uint64_t getRlLength(std::vector < std::vector < std::string > > const & bwtfilenames)
	{
		uint64_t fs = 0;
		for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
			fs += rl_decoder::getLength(bwtfilenames[i]);
		return fs;
	}

	static std::vector< std::vector<std::string> > stringVectorPack(std::vector<std::string> const & Vin)
	{
		std::vector< std::vector<std::string> > Vout;
		for ( uint64_t i = 0; i < Vin.size(); ++i )
			Vout.push_back(std::vector<std::string>(1,Vin[i]));
		return Vout;
	}

	static std::vector<std::string> parallelGapFragMerge(
		std::vector < std::vector < std::string > > const & bwtfilenames,
		std::vector < std::vector < std::string > > const & gapfilenames,
		std::string const & outputfilenameprefix,
		// std::string const & tempfilenameprefix,
		uint64_t const numthreads,
		uint64_t const lfblockmult,
		uint64_t const rlencoderblocksize
	)
	{
		// no bwt input files, create empty bwt file
		if ( ! bwtfilenames.size() )
		{
			std::string const outputfilename = outputfilenameprefix + ".bwt";
		
			#if defined(HUFRL)
			::libmaus::huffman::RLEncoderStd rlenc(outputfilename,0,               0,rlencoderblocksize);
			#else
			::libmaus::gamma::GammaRLEncoder rlenc(outputfilename,0 /* alphabet */,0,rlencoderblocksize);
			#endif
			
			rlenc.flush();
			
			return std::vector<std::string>(1,outputfilename);
		}
		// no gap files, rename input
		else if ( ! gapfilenames.size() )
		{
			assert ( bwtfilenames.size() == 1 );

			std::vector<std::string> outputfilenames;
			
			for ( uint64_t i = 0; i < bwtfilenames[0].size(); ++i )
			{
				std::ostringstream outputfilenamestr;
				outputfilenamestr << outputfilenameprefix << '_'
					<< std::setw(4) << std::setfill('0') << i << std::setw(0)
					<< ".bwt";
				std::string const outputfilename = outputfilenamestr.str();
				
				// ::libmaus::util::GetFileSize::copy(bwtfilenames[0][i],outputfilename);
				rename ( bwtfilenames[0][i].c_str(), outputfilename.c_str() );
				
				outputfilenames.push_back(outputfilename);
			}
			
			return outputfilenames;
		}
		// at least one gap file, merge
		else
		{
			#if defined(HUFGAP)
			typedef ::libmaus::huffman::GapDecoder gapfile_decoder_type;
			#else
			typedef ::libmaus::gamma::GammaGapDecoder gapfile_decoder_type;
			#endif

			#if !defined(HUFRL)
			unsigned int const albits = rl_decoder::getAlBits(bwtfilenames.front());
			#endif
			
			uint64_t const firstblockgapfilesize = gapfilenames.size() ? gapfile_decoder_type::getLength(gapfilenames[0]) : 0;
			assert ( firstblockgapfilesize );

			// uint64_t const fs = rl_decoder::getLength(bwtfilenames);
			uint64_t const fs = getRlLength(bwtfilenames);

			// first gap file meta information
			::libmaus::huffman::IndexDecoderDataArray::unique_ptr_type Pfgapidda(new ::libmaus::huffman::IndexDecoderDataArray(gapfilenames[0]));
			// first gap file
			gapfile_decoder_type::unique_ptr_type fgap(new gapfile_decoder_type(*Pfgapidda /* gapfilenames[0] */));
			// low marker
			uint64_t hlow = 0;
			// target g parts
			uint64_t const tgparts = numthreads*lfblockmult;
			// size of parts
			uint64_t const gpartsize = (fs + tgparts - 1) / tgparts;
			// maximum parts
			uint64_t const maxgparts = (fs + gpartsize - 1) / gpartsize;
			::libmaus::autoarray::AutoArray< ::libmaus::suffixsort::GapMergePacket> gmergepackets(maxgparts);
			// actual merge parts
			uint64_t actgparts = 0;
			uint64_t gs = 0;
			
			// std::cerr << "fs=" << fs << " tgparts=" << tgparts << " gpartsize=" << gpartsize << std::endl;
						
			// compute number of suffixes per gblock
			while ( hlow != firstblockgapfilesize )
			{
				uint64_t const gsrest = (fs-gs);
				uint64_t const gskip = std::min(gsrest,gpartsize);
				libmaus::huffman::KvInitResult kvinit;
				gapfile_decoder_type lgapdec(*Pfgapidda,gs + gskip, kvinit);
				
				// avoid small rest
				if ( fs - kvinit.kvoffset <= 1024 )
				{
					kvinit.kvoffset = fs;
					kvinit.koffset  = firstblockgapfilesize;
					// std::cerr << "+++ end" << std::endl;
				}
				// we did not end up on a gap array element, go to the next one
				else if ( kvinit.kvtarget )
				{
					kvinit.koffset += 1;
					kvinit.voffset += kvinit.kvtarget + lgapdec.peek();
					kvinit.kvoffset += (1 + kvinit.kvtarget + lgapdec.peek());
					
					// there is no new suffix after the last gap element, so deduct 1 if we reached the end of the array
					if ( kvinit.koffset == firstblockgapfilesize )
						kvinit.kvoffset -= 1;
					
					// std::cerr << "intermediate partial" << std::endl;
				}
				
				#if 1
				uint64_t const s = kvinit.kvoffset-gs;
				uint64_t const hhigh = kvinit.koffset;
				#else
				uint64_t hhigh = hlow;
				uint64_t s = 0;
				
				// sum up until we have reached at least gpartsize or the end of the file
				while ( hhigh != firstblockgapfilesize && s < gpartsize )
				{
					uint64_t const d = fgap->decode();
					s += (d + 1);
					hhigh += 1;
				}
								
				// last gap file value has no following suffix in block merged into
				if ( hhigh == firstblockgapfilesize )
					s -= 1;
				
				// if there is only one suffix left, then include it
				if ( hhigh+1 == firstblockgapfilesize && fgap->peek() == 0 )
				{
					fgap->decode();
					hhigh += 1;
				}

				std::cerr << "*** hhigh: " << hhigh << " kvinit.koffset: " << kvinit.koffset << " s=" << s << " kvinit.kvoffset-gs=" << kvinit.kvoffset-gs << std::endl;
				// kvinit.koffset, kvinit.kvoffset, rest: kvinit.kvtarget

				if ( actgparts >= maxgparts )
				{
					::libmaus::exception::LibMausException se;
					se.getStream() << "acgtgparts=" << actgparts << " >= maxgparts=" << maxgparts << std::endl;
					se.finish();
					throw se;
					// assert ( actgparts < maxgparts );
				}
				#endif
				
				// save interval on first gap array and number of suffixes on this interval
				gmergepackets[actgparts++] = ::libmaus::suffixsort::GapMergePacket(hlow,hhigh,s);
				
				// std::cerr << "got packet " << gmergepackets[actgparts-1] << " firstblockgapfilesize=" << firstblockgapfilesize << std::endl;
				
				// add suffixes in this block to global count
				gs += s;				
				
				// set new low marker	
				hlow = hhigh;
			}
			
			// we should have seen all the suffixes
			assert ( gs == fs );
			
			#if 0
			std::cerr << "actgparts=" << actgparts << std::endl;
			
			for ( uint64_t i = 0; i < actgparts; ++i )
				std::cerr << gmergepackets[i] << std::endl;
			#endif
			
			// compute prefix sums over number of suffixes per gblock
			::libmaus::autoarray::AutoArray<uint64_t> spref(actgparts+1,false);
			for ( uint64_t i = 0; i < actgparts; ++i )
				spref[i] = gmergepackets[i].s;
			spref.prefixSums();
			
			// compute prefix sums over number of suffixes per block used for each gpart
			::libmaus::autoarray::AutoArray < uint64_t > bwtusedcntsacc( (actgparts+1)*(gapfilenames.size()+1), false );
			
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for ( int64_t z = 0; z < static_cast<int64_t>(spref.size()); ++z )
			{
				#if 0
				std::cerr << "proc first: " << gmergepackets[z] << ",spref=" << spref[z] << std::endl;
				#endif
				
				// offset in first gap array
				uint64_t lspref = spref[z];
				
				// array of decoders
				::libmaus::autoarray::AutoArray < gapfile_decoder_type::unique_ptr_type > gapdecs(gapfilenames.size());

				for ( uint64_t j = 0; j < gapfilenames.size(); ++j )
				{
					::libmaus::huffman::KvInitResult kvinitresult;
					gapdecs[j] = UNIQUE_PTR_MOVE(
						gapfile_decoder_type::unique_ptr_type(
							new gapfile_decoder_type(
								gapfilenames[j],
								lspref,kvinitresult
							)
						)
					);
					
					// key offset for block z and file j
					bwtusedcntsacc [ j * (actgparts+1) + z ] = kvinitresult.koffset;

					#if 0
					std::cerr << "lspref=" << lspref << "," << kvinitresult << std::endl;				
					#endif

					// we should be on a key if j is the first file
					if ( j == 0 )
					{
						if ( kvinitresult.kvtarget != 0 )
						{
							std::cerr << "j=0 " << " z=" << z << " lspref=" << lspref <<
								" kvinitresult.koffset=" << kvinitresult.koffset <<
								" kvinitresult.voffset=" << kvinitresult.voffset <<
								" kvinitresult.kvoffset=" << kvinitresult.kvoffset <<
								" kvinitresult.kvtarget=" << kvinitresult.kvtarget << std::endl;
						}
						assert ( kvinitresult.kvtarget == 0 );
					}

					// offset for next gap file:
					// sum of values up to key lspref in this file + number of values used for next key
					lspref = kvinitresult.voffset + kvinitresult.kvtarget;

				}
				
				#if 0
				std::cerr << "lspref=" << lspref << std::endl;
				#endif
				
				// set end pointer
				bwtusedcntsacc [ gapfilenames.size() * (actgparts+1) + z ] = lspref;
			}
			
			// how many suffixes of each block do we use in each gpart?
			// (turn prefix sums in count per block)
			::libmaus::autoarray::AutoArray < uint64_t > bwtusedcnts( (gapfilenames.size()+1) * actgparts, false );
			for ( uint64_t block = 0; block < gapfilenames.size()+1; ++block )
			{
				uint64_t const * lbwtusedcntsacc = bwtusedcntsacc.begin() + block*(actgparts+1);

				for ( uint64_t i = 0; i < actgparts; ++i )
					bwtusedcnts [ block * actgparts + i ] = lbwtusedcntsacc[i+1]-lbwtusedcntsacc[i];
				
				assert ( lbwtusedcntsacc [ actgparts ] == rl_decoder::getLength(bwtfilenames[block]) );
			}

			// vector of output file names
			std::vector < std::string > gpartfrags(actgparts);

			// now use the block counts computed above
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for ( int64_t z = 0; z < static_cast<int64_t>(actgparts); ++z )
			{
				std::ostringstream ostr;
				ostr << outputfilenameprefix << "_" << std::setw(4) << std::setfill('0') << z << std::setw(0) << ".bwt";
				std::string const gpartfrag = ostr.str();
				::libmaus::util::TempFileRemovalContainer::addTempFile(gpartfrag);
				gpartfrags[z] = gpartfrag;
			
				#if 0
				std::cerr << gmergepackets[z] << ",spref=" << spref[z] << std::endl;
				#endif
				uint64_t lspref = spref[z];
				::libmaus::autoarray::AutoArray < gapfile_decoder_type::unique_ptr_type > gapdecoders(gapfilenames.size());
				::libmaus::autoarray::AutoArray< uint64_t > gapcur(gapfilenames.size());

				// set up gap file decoders at the proper offsets
				for ( uint64_t j = 0; j < gapfilenames.size(); ++j )
				{
					// sum up number of suffixes in later blocks for this gpart
					uint64_t suflat = 0;
					for ( uint64_t k = j+1; k < bwtfilenames.size(); ++k )
						suflat += bwtusedcnts [ k*actgparts + z ];
				
					::libmaus::huffman::KvInitResult kvinitresult;
					gapdecoders[j] = UNIQUE_PTR_MOVE(
						gapfile_decoder_type::unique_ptr_type(
							new gapfile_decoder_type(
								gapfilenames[j],
								lspref,kvinitresult
							)
						)
					);
					if ( suflat )
						gapcur[j] = gapdecoders[j]->decode();
					else
						gapcur[j] = 0;
						
					lspref = kvinitresult.voffset + kvinitresult.kvtarget;

					if ( j == 0 )
						assert ( kvinitresult.kvtarget == 0 );
				}				

				::libmaus::autoarray::AutoArray < uint64_t > bwttowrite(bwtfilenames.size(),false);
				::libmaus::autoarray::AutoArray < rl_decoder::unique_ptr_type > bwtdecoders(bwtfilenames.size());
				
				for ( uint64_t j = 0; j < bwtfilenames.size(); ++j )
				{
					uint64_t const bwtoffset = bwtusedcntsacc [ j * (actgparts+1) + z ];
					bwttowrite[j] = bwtusedcnts [ j * actgparts + z ];
					
					bwtdecoders[j] = UNIQUE_PTR_MOVE(
						rl_decoder::unique_ptr_type(
							new rl_decoder(bwtfilenames[j],bwtoffset)
						)
					);
					
					#if 0
					std::cerr << "block=" << j << " offset=" << bwtoffset << " bwttowrite=" << bwttowrite[j] << std::endl;
					#endif
				}
				
				uint64_t const totalbwt = std::accumulate(bwttowrite.begin(),bwttowrite.end(),0ull);

				#if defined(HUFRL)
				::libmaus::huffman::RLEncoderStd bwtenc(gpartfrag,0                    ,totalbwt,rlencoderblocksize);
				#else
				::libmaus::gamma::GammaRLEncoder bwtenc(gpartfrag,albits /* alphabet */,totalbwt,rlencoderblocksize);
				#endif
				
				// start writing loop
				uint64_t written = 0;
				while ( written < totalbwt )
				{
					// determine file we next read/decode from
					// this is the leftmost one with gap value 0 and still data to write
					uint64_t writeindex = bwtdecoders.size()-1;
					for ( 
						uint64_t i = 0; 
						i < gapcur.size(); 
						++i
					)
						if ( 
							(! gapcur[i]) 
							&& 
							(bwttowrite[i])
						)
						{
							writeindex = i;
							break;
						}
				
					// sanity check
					if ( ! bwttowrite[writeindex] )
					{
						assert ( bwttowrite[writeindex] );
					}
					
					// adjust counters
					written++;
					bwttowrite[writeindex]--;
					// get next gap value if block is not the last one
					if ( bwttowrite[writeindex] && writeindex < gapcur.size() )
						gapcur[writeindex] = gapdecoders[writeindex]->decode();
					
					// copy symbol
					uint64_t const sym = bwtdecoders[writeindex]->decode();
					bwtenc.encode(sym);

					// adjust gap values of blocks to the left
					for ( uint64_t i = 0; i < writeindex; ++i )
						if ( bwttowrite[i] )
						{
							assert ( gapcur[i] > 0 );
							gapcur[i]--;
						}							
				}
				// std::cerr << "(1)";
				
				// all data should have been written now
				for ( uint64_t i = 0; i < bwttowrite.size(); ++i )
					assert ( !bwttowrite[i] );
					
				bwtenc.flush();
			}

			return gpartfrags;
		}
	}
	
	template<typename gap_iterator>
	static uint64_t mergeIsa(
		std::string const & oldmergedisaname, // mergedisaname
		std::string const & newmergedisaname, // blockresults.files.sampledisa
		std::string const & mergedmergedisaname, // newmergedisaname
		uint64_t const blockstart, // blockstart
		gap_iterator & Gc,
		uint64_t const gn
	)
	{
		::libmaus::timing::RealTimeClock rtc;
		
		// merge sampled inverse suffix arrays		
		std::cerr << "[V] merging sampled inverse suffix arrays...";
		rtc.start();
		::libmaus::aio::SynchronousGenericInput<uint64_t> SGIISAold(oldmergedisaname,16*1024);
		::libmaus::aio::SynchronousGenericInput<uint64_t> SGIISAnew(newmergedisaname,16*1024);
		::libmaus::aio::SynchronousGenericOutput<uint64_t> SGOISA(mergedmergedisaname,16*1024);
		
		uint64_t s = 0;
		uint64_t blockp0rank = 0;
		// scan the gap array
		for ( uint64_t i = 0; i < gn; ++i )
		{
			s += Gc.get(); // *(Gc++);

			// while old suffixes are smaller than next new suffix
			int64_t re;
			while ( (re=SGIISAold.peek()) >= 0 && re < static_cast<int64_t>(s) )
			{
				// rank (add number of new suffixes processed before)
				uint64_t const r = SGIISAold.get() + i;
				// position (copy as is)
				uint64_t const p = SGIISAold.get();
			
				SGOISA . put ( r );
				SGOISA . put ( p );
				
				// std::cerr << "copying r=" << r << " p=" << p << std::endl;
			}

			if ( SGIISAnew.peek() == static_cast<int64_t>(i) )
			{
				// add number of old suffixes to rank
				uint64_t const r = SGIISAnew.get() + s;
				// keep absolute position as is
				uint64_t const p = SGIISAnew.get();
				
				SGOISA . put ( r );
				SGOISA . put ( p );
				
				if ( p == blockstart )
				{
					blockp0rank = r;
					//std::cerr << "p=" << p << " r=" << r << std::endl;
				}			

				// std::cerr << "copying r=" << r << " p=" << p << std::endl;
			}
		}

		SGOISA.flush();
		std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
	
		return blockp0rank;
	}

	static void saveGapFile(
		::libmaus::autoarray::AutoArray<uint32_t> const & G,
		std::string const & gapfile
	)
	{
		::libmaus::timing::RealTimeClock rtc;
	
		std::cerr << "[V] saving gap file...";
		rtc.start();
		#if defined(HUFGAP)
		::libmaus::util::Histogram gaphist;
		for ( uint64_t i = 0; i < G.size(); ++i )
			gaphist ( G[i] );
		::libmaus::huffman::GapEncoder GE(gapfile,gaphist,G.size());
		GE.encode(G.begin(),G.end());
		GE.flush();
		#else
		::libmaus::gamma::GammaGapEncoder GE(gapfile);
		GE.encode(G.begin(),G.end());
		#endif
		std::cerr << "done in time " << rtc.getElapsedSeconds() << std::endl;	
	}

	static void concatenateGT(
		std::vector<std::string> const & gtpartnames,
		std::string const & newgtpart, // blockresults.files.gt
		std::string const & newmergedgtname
	)
	{
		::libmaus::timing::RealTimeClock rtc;
		rtc.start();
		
		std::vector< std::string > allfiles(gtpartnames.begin(),gtpartnames.end());
		allfiles.push_back(newgtpart);
		
		#if 0
		// encoder for new gt stream
		libmaus::bitio::BitVectorOutput GTHEFref(newmergedgtname);
		libmaus::bitio::BitVectorInput BVI(allfiles);
		uint64_t const tn = libmaus::bitio::BitVectorInput::getLength(allfiles);
		for ( uint64_t i = 0; i < tn; ++i )
			GTHEFref.writeBit(BVI.readBit());
		GTHEFref.flush();
		#endif

		libmaus::bitio::BitVectorOutput GTHEF(newmergedgtname);
		
		unsigned int prevbits = 0;
		uint64_t prev = 0;

		// append part streams
		for ( uint64_t z = 0; z < allfiles.size(); ++z )
		{
			uint64_t const n = libmaus::bitio::BitVectorInput::getLength(allfiles[z]);
			uint64_t const fullwords = n / 64;
			uint64_t const restbits = n-fullwords*64;
			libmaus::aio::SynchronousGenericInput<uint64_t> SGI(allfiles[z],8192);
			uint64_t v = 0;
			
			if ( !prevbits )
			{
				
				// copy full words
				for ( uint64_t i = 0; i < fullwords; ++i )
				{
					bool const ok = SGI.getNext(v);
					assert ( ok );
					GTHEF.SGO.put(v);
				}
				
				
				// get next word
				if ( restbits )
				{
					SGI.getNext(prev);
					// move bits to top of word
					prev <<= (64-restbits);
					// number of bits left
					prevbits = restbits;
				}
				else
				{
					prevbits = 0;
					prev = 0;
				}
			}
			else
			{
				// process full words
				for ( uint64_t i = 0; i < fullwords; ++i )
				{
					bool const ok = SGI.getNext(v);
					assert ( ok );
					
					GTHEF.SGO.put(prev | (v >> prevbits));
					prev = v << (64-prevbits);
				}
				
				// if there are bits left in this stream
				if ( restbits )
				{
					// get next word
					SGI.getNext(v);
					// move to top
					v <<= (64 - restbits);

					// rest with previous rest fill more than a word
					if ( restbits + prevbits >= 64 )
					{
						GTHEF.SGO.put(prev | (v >> prevbits));
						prev = v << (64-prevbits);
						prevbits = restbits + prevbits - 64;
					}
					else
					{
						prev = prev | (v >> prevbits);
						prevbits = prevbits + restbits;
					}
				}
				// leave as is
				else
				{
				
				}
			}
			
			remove(gtpartnames[z].c_str());
		}
		
		for ( uint64_t i = 0; i < prevbits; ++i )
			GTHEF.writeBit((prev >> (63-i)) & 1);
		
		// flush gt stream
		GTHEF.flush();
		
		std::cerr << "[V] concatenated bit vectors in time " << rtc.getElapsedSeconds() << std::endl;
	}
	
	struct GapArrayComputationResult
	{
		::libmaus::autoarray::AutoArray<uint32_t> G;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		GapArrayComputationResult()
		: zactive(0)
		{
		
		}
		
		GapArrayComputationResult(
			::libmaus::autoarray::AutoArray<uint32_t> & rG, 
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: G(rG), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};

	struct GapArrayByteComputationResult
	{
		libmaus::suffixsort::GapArrayByte::shared_ptr_type G;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		GapArrayByteComputationResult()
		: zactive(0)
		{
		
		}
		
		GapArrayByteComputationResult(
			libmaus::suffixsort::GapArrayByte::shared_ptr_type rG, 
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: G(rG), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};

	struct SparseGapArrayComputationResult
	{
		// name of gap file
		std::vector<std::string> fn;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		SparseGapArrayComputationResult()
		: zactive(0)
		{
		
		}
		
		SparseGapArrayComputationResult(
			std::vector < std::string > const & rfn,
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: fn(rfn), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};
	
	static libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ensureWaveletTreeGenerated(::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults)
	{
		// generate wavelet tree from request if necessary
		if ( ! libmaus::util::GetFileSize::fileExists(blockresults.getFiles().getHWT()) )
		{
			libmaus::timing::RealTimeClock rtc; rtc.start();
			std::cerr << "[V] Generating HWT for gap file computation...";
			assert ( libmaus::util::GetFileSize::fileExists(blockresults.getFiles().getHWTReq() ) );	
			RlToHwtTermRequest::unique_ptr_type ureq(RlToHwtTermRequest::load(blockresults.getFiles().getHWTReq()));
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type	ptr(ureq->dispatch());
			remove ( blockresults.getFiles().getHWTReq().c_str() );
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			return UNIQUE_PTR_MOVE(ptr);			
		}
		else
		{
			libmaus::aio::CheckedInputStream CIS(blockresults.getFiles().getHWT());
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(new libmaus::wavelet::ImpCompactHuffmanWaveletTree(CIS));
			return UNIQUE_PTR_MOVE(ptr);
		}
	}

	static GapArrayComputationResult computeGapArray(
		std::string const & fn, // name of text file
		uint64_t const fs, // length of text file in symbols
		uint64_t const blockstart, // start offset
		uint64_t const cblocksize, // block size
		uint64_t const nextblockstart, // start of next block (mod fs)
		uint64_t const mergeprocrightend, // right end of merged area
		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults, // information on block
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		::libmaus::lf::DArray * const accD, // accumulated symbol freqs for block
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks // lf starting points
	)
	{
		// gap array
		::libmaus::autoarray::AutoArray<uint32_t> G(cblocksize+1);
		
		// set up lf mapping
		::libmaus::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
		IHWL.D = D.D;
		assert ( cblocksize == IHWL.n );
		
		// rank of position 0 in this block (for computing new gt array/stream)
		uint64_t const lp0 = blockresults.getBlockP0Rank();
		
		// last symbol in this block
		int64_t const firstblocklast = input_types_type::linear_wrapper::getSymbolAtPosition(fn,(nextblockstart+fs-1)%fs);

		/**
		 * array of absolute positions
		 **/
		uint64_t const zactive = zblocks.size();
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus::timing::RealTimeClock rtc; 
		rtc.start();
		#if defined(_OPENMP) && defined(LIBMAUS_HAVE_SYNC_OPS)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t z = 0; z < static_cast<int64_t>(zactive); ++z )
		{
			::libmaus::timing::RealTimeClock subsubrtc; subsubrtc.start();

			::libmaus::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
			
			std::string const gtpartname = newmergedgtname + "_" + ::libmaus::util::NumberSerialisation::formatNumber(z,4) + ".gt";
			::libmaus::util::TempFileRemovalContainer::addTempFile(gtpartname);
			gtpartnames[z] = gtpartname;
			#if 0
			::libmaus::huffman::HuffmanEncoderFileStd GTHEF(gtpartname);
			#endif
			libmaus::bitio::BitVectorOutput GTHEF(gtpartname);

			#if 0
			::libmaus::bitio::BitStreamFileDecoder gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			#endif
			libmaus::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			
			typename input_types_type::circular_reverse_wrapper CRWR(fn,zblock.getZAbsPos() % fs);
			uint64_t r = zblock.getZRank();

			uint64_t const zlen = zabsblockpos [ z ] - zabsblockpos [z+1];
			
			for ( uint64_t i = 0; i < zlen; ++i )
			{
				GTHEF.writeBit(r > lp0);
				
				int64_t const sym = CRWR.get();			
				bool const gtf = gtfile.readBit();

				r = IHWL.step(sym,r) + ((sym == firstblocklast)?gtf:0);
			
				#if defined(LIBMAUS_HAVE_SYNC_OPS)
				__sync_fetch_and_add(G.begin()+r,1);
				#else
				G[r]++;
				#endif
			}
			
			GTHEF.flush();
		}
		std::cerr << "[V] computed gap array in time " << rtc.getElapsedSeconds() << std::endl;

		return GapArrayComputationResult(G,gtpartnames,zactive,zabsblockpos);
	}

	static GapArrayComputationResult computeGapArray(
		std::string const & fn, // name of text file
		uint64_t const fs, // length of text file in symbols
		MergeStrategyMergeGapRequest const & msmgr, // merge request
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		::libmaus::lf::DArray * const accD // accumulated symbol freqs for block
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		// use gap object's zblocks vector
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
		return computeGapArray(fn,fs,blockstart,cblocksize,nextblockstart,mergeprocrightend,
			blockresults,mergedgtname,newmergedgtname,accD,zblocks);
	}

	static GapArrayByteComputationResult computeGapArrayByte(
		std::string const & fn, // name of text file
		uint64_t const fs, // length of text file in symbols
		uint64_t const blockstart, // start offset
		uint64_t const cblocksize, // block size
		uint64_t const nextblockstart, // start of next block (mod fs)
		uint64_t const mergeprocrightend, // right end of merged area
		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults, // information on block
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		std::string const & gapoverflowtmpfilename, // gap overflow tmp file
		::libmaus::lf::DArray * const accD, // accumulated symbol freqs for block
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks // lf starting points
	)
	{
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif

		// gap array
		::libmaus::suffixsort::GapArrayByte::shared_ptr_type pG(
			new ::libmaus::suffixsort::GapArrayByte(
				cblocksize+1,
				512, /* number of overflow words per thread */
				numthreads,
				gapoverflowtmpfilename
			)
		);
		::libmaus::suffixsort::GapArrayByte & G = *pG;
		
		// set up lf mapping
		::libmaus::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
		IHWL.D = D.D;
		assert ( cblocksize == IHWL.n );
		
		// rank of position 0 in this block (for computing new gt array/stream)
		uint64_t const lp0 = blockresults.getBlockP0Rank();
		
		// last symbol in this block
		int64_t const firstblocklast = input_types_type::linear_wrapper::getSymbolAtPosition(fn,(nextblockstart+fs-1)%fs);

		/**
		 * array of absolute positions
		 **/
		uint64_t const zactive = zblocks.size();
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus::timing::RealTimeClock rtc; 
		rtc.start();
		#if defined(_OPENMP) && defined(LIBMAUS_HAVE_SYNC_OPS)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t z = 0; z < static_cast<int64_t>(zactive); ++z )
		{
			::libmaus::timing::RealTimeClock subsubrtc; subsubrtc.start();

			::libmaus::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
			
			std::string const gtpartname = newmergedgtname + "_" + ::libmaus::util::NumberSerialisation::formatNumber(z,4) + ".gt";
			::libmaus::util::TempFileRemovalContainer::addTempFile(gtpartname);
			gtpartnames[z] = gtpartname;
			#if 0
			::libmaus::huffman::HuffmanEncoderFileStd GTHEF(gtpartname);
			#endif
			libmaus::bitio::BitVectorOutput GTHEF(gtpartname);

			#if 0
			::libmaus::bitio::BitStreamFileDecoder gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			#endif
			libmaus::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			
			typename input_types_type::circular_reverse_wrapper CRWR(fn,zblock.getZAbsPos() % fs);
			uint64_t r = zblock.getZRank();

			uint64_t const zlen = zabsblockpos [ z ] - zabsblockpos [z+1];
			
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif
			
			for ( uint64_t i = 0; i < zlen; ++i )
			{
				GTHEF.writeBit(r > lp0);
				
				int64_t const sym = CRWR.get();			
				bool const gtf = gtfile.readBit();

				r = IHWL.step(sym,r) + ((sym == firstblocklast)?gtf:0);
			
				if ( G(r) )
					G(r,tid);
			}
			
			GTHEF.flush();
		}
		std::cerr << "[V] computed gap array in time " << rtc.getElapsedSeconds() << std::endl;
		
		G.flush();

		return GapArrayByteComputationResult(pG,gtpartnames,zactive,zabsblockpos);
	}

	static GapArrayByteComputationResult computeGapArrayByte(
		std::string const & fn, // name of text file
		uint64_t const fs, // length of text file in symbols
		MergeStrategyMergeGapRequest const & msmgr, // merge request
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		std::string const & gapoverflowtmpfilename, // gap overflow tmp file
		::libmaus::lf::DArray * const accD // accumulated symbol freqs for block
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		// use gap object's zblocks vector
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
		return computeGapArrayByte(fn,fs,blockstart,cblocksize,nextblockstart,mergeprocrightend,
			blockresults,mergedgtname,newmergedgtname,gapoverflowtmpfilename,accD,zblocks);
	}
	
	struct ZNext
	{
		uint64_t znext;
		uint64_t znextcount;
		libmaus::parallel::OMPLock lock;
		
		ZNext(uint64_t const rznextcount) : znext(0), znextcount(rznextcount) {}
		
		bool getNext(uint64_t & next)
		{
			libmaus::parallel::ScopeLock slock(lock);
			
			if ( znext == znextcount )
				return false;
			else
			{
				next = znext++;
				return true;
			}
		}
	};

	static SparseGapArrayComputationResult computeSparseGapArray(
		std::string const & fn,
		uint64_t const fs,
		uint64_t const blockstart,
		uint64_t const cblocksize,
		uint64_t const nextblockstart,
		uint64_t const mergeprocrightend,
		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults,
		std::vector<std::string> const & mergedgtname,
		std::string const & newmergedgtname,
		::libmaus::lf::DArray * const accD,
		//
		std::string const & outputgapfilename,
		std::string const & tmpfileprefix,
		uint64_t const maxmem,
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks
	)
	{
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		uint64_t const memperthread = (maxmem + numthreads-1)/numthreads;
		uint64_t const wordsperthread = ( memperthread + sizeof(uint64_t) - 1 ) / sizeof(uint64_t);
		// uint64_t const wordsperthread = 600;
		uint64_t const parcheck = 64*1024;
		
		libmaus::autoarray::AutoArray< libmaus::autoarray::AutoArray<uint64_t> > GG(numthreads,false);
		for ( uint64_t i = 0; i < numthreads; ++i )
			GG[i] = libmaus::autoarray::AutoArray<uint64_t>(wordsperthread,false);
		
		libmaus::util::TempFileNameGenerator tmpgen(tmpfileprefix,3);
		libmaus::gamma::SparseGammaGapMultiFileLevelSet SGGFS(tmpgen,numthreads);
	
		// set up lf mapping
		::libmaus::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = 
			ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
		IHWL.D = D.D;
		assert ( cblocksize == IHWL.n );
		
		// rank of position 0 in this block (for computing new gt array/stream)
		uint64_t const lp0 = blockresults.getBlockP0Rank();
		
		// last symbol in this block
		int64_t const firstblocklast = input_types_type::linear_wrapper::getSymbolAtPosition(fn,(nextblockstart+fs-1)%fs);
		
		/**
		 * array of absolute positions
		 **/
		uint64_t const zactive = zblocks.size();
		::libmaus::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus::timing::RealTimeClock rtc; 
		rtc.start();
		
		ZNext znext(zactive);
		
		uint64_t termcnt = 0;
		libmaus::parallel::OMPLock termcntlock;
		
		libmaus::parallel::PosixSemaphore qsem; // queue semaphore
		libmaus::parallel::PosixSemaphore tsem; // term semaphore
		libmaus::parallel::PosixSemaphore globsem; // meta semaphore for both above
		libmaus::parallel::LockedBool termflag(false);
		libmaus::parallel::LockedBool qterm(false);
		
		SGGFS.registerMergePackSemaphore(&qsem);
		SGGFS.registerMergePackSemaphore(&globsem);
		SGGFS.registerTermSemaphore(&tsem);
		SGGFS.registerTermSemaphore(&globsem);
		SGGFS.setTermSemCnt(numthreads);
				
		#pragma omp parallel
		{
			uint64_t z = 0;
			
			while ( znext.getNext(z) )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif
				
				uint64_t * const Ga = GG[tid].begin();
				uint64_t *       Gc = Ga;
				uint64_t * const Ge = GG[tid].end();
			
				::libmaus::timing::RealTimeClock subsubrtc; subsubrtc.start();
				
				::libmaus::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
				
				std::string const gtpartname = newmergedgtname + "_" + ::libmaus::util::NumberSerialisation::formatNumber(z,4) + ".gt";
				::libmaus::util::TempFileRemovalContainer::addTempFile(gtpartname);
				gtpartnames[z] = gtpartname;
				libmaus::bitio::BitVectorOutput GTHEF(gtpartname);

				libmaus::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
				
				typename input_types_type::circular_reverse_wrapper CRWR(fn,zblock.getZAbsPos() % fs);
				uint64_t r = zblock.getZRank();

				uint64_t const zlen = zabsblockpos [ z ] - zabsblockpos [z+1];
				uint64_t const fullblocks = zlen/(Ge-Ga);
				uint64_t const rest = zlen - fullblocks * (Ge-Ga);
				
				for ( uint64_t b = 0; b < fullblocks; ++b )
				{
					Gc = Ga;
				
					while ( Gc != Ge )
					{
						uint64_t * const Te = Gc + std::min(parcheck,static_cast<uint64_t>(Ge-Gc));
					
						for ( ; Gc != Te; ++Gc )
						{
							GTHEF.writeBit(r > lp0);			
							int64_t const sym = CRWR.get();			
							bool const gtf = gtfile.readBit();
							r = IHWL.step(sym,r) + ((sym == firstblocklast)?gtf:0);
							*Gc = r;
						}

						while ( globsem.trywait() )
						{
							qsem.wait();
							SGGFS.checkMergeSingle();
						}
					}

					std::string const tfn = tmpgen.getFileName();
					libmaus::gamma::SparseGammaGapBlockEncoder::encodeArray(Ga,Gc,tfn);
					SGGFS.putFile(std::vector<std::string>(1,tfn));		

					while ( globsem.trywait() )
					{
						qsem.wait();
						SGGFS.checkMergeSingle();
					}
				}
				
				if ( rest )
				{
					Gc = Ga;
				
					while ( Gc != Ga + rest )
					{
						uint64_t * const Te = Gc + std::min(parcheck,static_cast<uint64_t>((Ga+rest)-Gc));
					
						for ( ; Gc != Te; ++Gc )
						{
							GTHEF.writeBit(r > lp0);			
							int64_t const sym = CRWR.get();			
							bool const gtf = gtfile.readBit();
							r = IHWL.step(sym,r) + ((sym == firstblocklast)?gtf:0);
							*Gc = r;				
						}

						while ( globsem.trywait() )
						{
							qsem.wait();
							SGGFS.checkMergeSingle();
						}
					}

					std::string const tfn = tmpgen.getFileName();
					libmaus::gamma::SparseGammaGapBlockEncoder::encodeArray(Ga,Gc,tfn);
					SGGFS.putFile(std::vector<std::string>(1,tfn));		

					while ( globsem.trywait() )
					{
						qsem.wait();
						SGGFS.checkMergeSingle();
					}
				}
				
				GTHEF.flush();

				while ( globsem.trywait() )
				{
					qsem.wait();
					SGGFS.checkMergeSingle();
				}
			}	
			
			{
				termcntlock.lock();
				if ( ++termcnt == numthreads )
					termflag.set(true);					
				termcntlock.unlock();				
			}
			
			bool running = true;
			while ( running )
			{
				if ( termflag.get() && (!qterm.get()) && SGGFS.isMergingQueueEmpty() )
				{
					for ( uint64_t i = 0; i < numthreads; ++i )
					{
						tsem.post();
						globsem.post();
					}
					
					qterm.set(true);
				}
				
				// std::cerr << "waiting for glob...";
				globsem.wait();
				// std::cerr << "done." << std::endl;

				if ( qsem.trywait() )
					SGGFS.checkMergeSingle();
				else
				{
					bool const tsemok = tsem.trywait();
					assert ( tsemok );
					running = false;
				}
			}
		}
		
		std::vector<std::string> const outputgapfilenames = SGGFS.mergeToDense(outputgapfilename,cblocksize+1);
		
		std::cerr << "[V] computed gap array in time " << rtc.getElapsedSeconds() << std::endl;

		return SparseGapArrayComputationResult(outputgapfilenames,gtpartnames,zactive,zabsblockpos);
	}

	static SparseGapArrayComputationResult computeSparseGapArray(
		std::string const & fn, // name of text file
		uint64_t const fs, // length of text file in symbols
		MergeStrategyMergeGapRequest const & msmgr, // merge request
		uint64_t const ihwtspace,
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		::libmaus::lf::DArray * const accD, // accumulated symbol freqs for block
		std::string const & outputgapfilename,
		std::string const & tmpfileprefix,
		uint64_t const maxmem
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		std::vector < ::libmaus::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
		return computeSparseGapArray(fn,fs,blockstart,cblocksize,nextblockstart,mergeprocrightend,
			blockresults,mergedgtname,newmergedgtname,accD,
			outputgapfilename,
			tmpfileprefix,
			maxmem-ihwtspace,
			zblocks
		);
	}
	
	static std::vector<std::string> stringVectorAppend(std::vector<std::string> V, std::vector<std::string> const & W)
	{
		for ( uint64_t i = 0; i < W.size(); ++i )
			V.push_back(W[i]);
		return V;
	}

	static void mergeBlocks(
		MergeStrategyMergeInternalBlock & mergereq,
		std::string const fn,
		uint64_t const fs,
		std::string const tmpfilenamebase,
		uint64_t const rlencoderblocksize,
		uint64_t const lfblockmult,
		uint64_t const numthreads,
		::std::map<int64_t,uint64_t> const & /* chist */,
		uint64_t const bwtterm,
		std::string const & huftreefilename
	)
	{
		std::cerr << "[V] Merging BWT blocks MergeStrategyMergeInternalBlock." << std::endl;
		
		assert ( mergereq.children.size() > 1 );
		assert ( mergereq.children.size() == mergereq.gaprequests.size()+1 );

		// std::cerr << "[V] Merging BWT blocks with gapmembound=" << gapmembound << std::endl;

		/*
		 * remove unused file
		 */
		remove ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
		// fill result structure
		result.setBlockStart( mergereq.children[0]->sortresult.getBlockStart() );
		result.setCBlockSize( 0 );
		for ( uint64_t i = 0; i < mergereq.children.size(); ++i )
			result.setCBlockSize( result.getCBlockSize() + mergereq.children[i]->sortresult.getCBlockSize() );
		// set up temp file names
		// of output bwt,
		// sampled inverse suffix array filename,
		// gt bit array,
		// huffman shaped wavelet tree and
		// histogram
		result.setTempPrefixAndRegisterAsTemp(tmpfilenamebase + "_out",0 /* no preset bwt file names */, 0 /* no preset gt file names */);

		// if we merge only two blocks together, then we do not need to write the gap array to disk
		if ( mergereq.children.size() == 2 )
		{
			// std::cerr << "** WHITEBOX INTERNAL 1 **" << std::endl;
		
			std::string const & sblockhist = mergereq.children[1]->sortresult.getFiles().getHist();
			// load char histogram for last/second block (for merging)
			::libmaus::lf::DArray::unique_ptr_type accD(new ::libmaus::lf::DArray(
				sblockhist
				)
			);
			// first block
			::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = 
				mergereq.children[0]->sortresult;

			// start of first block
			uint64_t const blockstart = blockresults.getBlockStart();
			// size of first block
			uint64_t const cblocksize = blockresults.getCBlockSize();

			// compute gap array
			GapArrayComputationResult const GACR = computeGapArray(
				fn,fs,*(mergereq.gaprequests[0]),
				mergereq.children[1]->sortresult.getFiles().getGT(), // previous gt files
				tmpfilenamebase + "_gparts", // new gt files
				accD.get()
			);

			#if 0
			// concatenate gt vectors
			concatenateGT(
				GACR.gtpartnames,
				blockresults.getFiles().getGT(),
				result.getFiles().getGT()
			);
			#endif
			
			std::vector<std::string> oldgtnames;
			for ( uint64_t i = 0; i < blockresults.getFiles().getGT().size(); ++i )
			{
				std::ostringstream ostr;
				ostr << tmpfilenamebase << "_renamed_" << std::setw(6) << std::setfill('0') << i << std::setw(0) << ".gt";
				std::string const renamed = ostr.str();
				oldgtnames.push_back(ostr.str());
				::libmaus::util::TempFileRemovalContainer::addTempFile(renamed);
				rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
			}

			result.setGT(stringVectorAppend(GACR.gtpartnames,oldgtnames));

			// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
			libmaus::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
			result.setBlockP0Rank( mergeIsa(
				mergereq.children[1]->sortresult.getFiles().getSampledISA(), // old sampled isa
				blockresults.getFiles().getSampledISA(), // new sampled isa
				result.getFiles().getSampledISA(),blockstart,mergeGO/*GACR.G.begin()*/,cblocksize+1 /* GACR.G.size() */
			) );
			
			::libmaus::timing::RealTimeClock rtc; rtc.start();
			std::cerr << "[V] merging BWTs...";
			
			//
			uint64_t const totalsuf = result.getCBlockSize();
			// number of packets
			uint64_t const numpack = numthreads;
			// suffixes per thread
			uint64_t const tpacksize = (totalsuf + numpack-1)/numpack;
			//
			uint64_t ilow = 0;
			// intervals on G
			std::vector < std::pair<uint64_t,uint64_t> > wpacks;
			std::vector < std::string > encfilenames;
			// prefix sums over G
			std::vector < uint64_t > P;
			P.push_back(0);
			
			// std::cerr << "(computing work packets...";
			::libmaus::timing::RealTimeClock wprtc; wprtc.start();
			while ( ilow != GACR.G.size() )
			{
				uint64_t s = 0;
				uint64_t ihigh = ilow;
				
				while ( ihigh != GACR.G.size() && s < tpacksize )
					s += (GACR.G[ihigh++]+1);
				
				uint64_t const p = s-(ihigh-ilow);
				
				if ( ihigh+1 == GACR.G.size() && GACR.G[ihigh] == 0 )
					ihigh++;

				// std::cerr << "[" << ilow << "," << ihigh << ")" << std::endl;
				
				assert ( p == std::accumulate(GACR.G.begin()+ilow,GACR.G.begin()+ihigh,0ull) );

				P.push_back(P.back() + p);
				wpacks.push_back(std::pair<uint64_t,uint64_t>(ilow,ihigh));
				encfilenames.push_back(
					tmpfilenamebase 
					// result.getFiles().getBWT() 
					+ "_"
					+ ::libmaus::util::NumberSerialisation::formatNumber(encfilenames.size(),6)
					+ ".bwt"
				);
				::libmaus::util::TempFileRemovalContainer::addTempFile(encfilenames.back());
				ilow = ihigh;
			}
			assert ( wpacks.size() <= numthreads );
			// std::cerr << "done,time=" << wprtc.getElapsedSeconds() << ")";
			
			// std::cerr << "(setting up IDDs...";
			wprtc.start();
			#if !defined(HUFRL)
			unsigned int const albits = rl_decoder::getAlBits(mergereq.children[0]->sortresult.getFiles().getBWT());
			#endif
			::libmaus::huffman::IndexDecoderDataArray IDD0(
				mergereq.children[0]->sortresult.getFiles().getBWT());
			::libmaus::huffman::IndexDecoderDataArray IDD1(
				mergereq.children[1]->sortresult.getFiles().getBWT());
			
			::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV0 = ::libmaus::huffman::IndexLoader::loadAccIndex(
				mergereq.children[0]->sortresult.getFiles().getBWT()
			);
			::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV1 = ::libmaus::huffman::IndexLoader::loadAccIndex(
				mergereq.children[1]->sortresult.getFiles().getBWT()
			);
						
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for ( int64_t b = 0; b < static_cast<int64_t>(wpacks.size()); ++b )
			{
				uint64_t const ilow = wpacks[b].first;
				uint64_t const ihigh = wpacks[b].second;
				
				if ( ilow != ihigh )
				{
					bool const islast = (ihigh == GACR.G.size());
					std::string const encfilename = encfilenames[b];

					#if defined(HUFRL)
					rl_decoder leftrlin(IDD0,ilow);
					rl_decoder rightrlin(IDD1,P[b]);
					#else
					rl_decoder leftrlin(IDD0,IECV0.get(),ilow);
					rl_decoder rightrlin(IDD1,IECV1.get(),P[b]);
					#endif
					
					uint64_t const outsuf = (ihigh-ilow)-(islast?1:0) + (P[b+1]-P[b]);

					#if defined(HUFRL)
					::libmaus::huffman::RLEncoderStd bwtenc(encfilename,0     ,outsuf,rlencoderblocksize);
					#else
					::libmaus::gamma::GammaRLEncoder bwtenc(encfilename,albits,outsuf,rlencoderblocksize);
					#endif
				
					if ( islast )
					{
						for ( uint64_t j = ilow; j < ihigh-1; ++j )
						{
							for ( uint64_t i = 0; i < GACR.G[j]; ++i )
								bwtenc.encode(rightrlin.decode());
								
							bwtenc.encode(leftrlin.decode());
						}

						for ( uint64_t i = 0; i < GACR.G[cblocksize]; ++i )
							bwtenc.encode(rightrlin.decode());						
					}
					else
					{
						for ( uint64_t j = ilow; j < ihigh; ++j )
						{
							for ( uint64_t i = 0; i < GACR.G[j]; ++i )
								bwtenc.encode(rightrlin.decode());
								
							bwtenc.encode(leftrlin.decode());
						}
					}
					
					bwtenc.flush();
				}
			}
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;

			#if 0
			std::cerr << "[V] concatenating bwt parts...";			
			rtc.start();
			#if defined(HUFRL)
			::libmaus::huffman::RLEncoderStd::concatenate(encfilenames,result.getFiles().getBWT());
			#else
			::libmaus::gamma::GammaRLEncoder::concatenate(encfilenames,result.getFiles().getBWT());
			#endif
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			result.setBWT(encfilenames);
			
			#if 0
			std::cerr << "[V] removing tmp files...";			
			rtc.start();
			for ( uint64_t i = 0; i < encfilenames.size(); ++i )
				remove ( encfilenames[i].c_str() );
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			// save histogram
			// std::cerr << "[V] saving histogram...";
			rtc.start();
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));
			// std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
		}
		else
		{
			// std::cerr << "** WHITEBOX INTERNAL 2 **" << std::endl;
			
			std::vector < std::string > gapfilenames;
			std::vector < std::vector<std::string> > bwtfilenames;
			for ( uint64_t bb = 0; bb < mergereq.children.size(); ++bb )
			{
				// gap file name
				if ( bb+1 < mergereq.children.size() )
				{
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					::libmaus::util::TempFileRemovalContainer::addTempFile(newgapname);
					gapfilenames.push_back(newgapname);
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}

			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus::lf::DArray::unique_ptr_type accD(new ::libmaus::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
				// gap file
				std::string const gapfile = gapfilenames[bx];

				// start of this block
				uint64_t const blockstart = blockresults.getBlockStart();
				// size of this block
				uint64_t const cblocksize = blockresults.getCBlockSize();

				// compute gap array
				GapArrayComputationResult const GACR = computeGapArray(
					fn,fs,*(mergereq.gaprequests[bx]),
					mergedgtname,
					newmergedgtname,
					accD.get()
				);

				// save the gap file			
				saveGapFile(GACR.G,gapfile);

				// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
				libmaus::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
				result.setBlockP0Rank( mergeIsa(mergedisaname,blockresults.getFiles().getSampledISA(),newmergedisaname,blockstart,mergeGO /*GACR.G.begin()*/,cblocksize+1 /*GACR.G.size()*/) );
				
				#if 0
				// concatenate gt vectors
				concatenateGT(GACR.gtpartnames,blockresults.getFiles().getGT(),newmergedgtname);
				#endif

				// rename files
				std::vector<std::string> oldgtnames;
				for ( uint64_t i = 0; i < blockresults.getFiles().getGT().size(); ++i )
				{
					std::ostringstream ostr;
					ostr << tmpfilenamebase 
						<< "_renamed_" 
						<< std::setw(6) << std::setfill('0') << bx << std::setw(0) 
						<< "_"
						<< std::setw(6) << std::setfill('0') << i << std::setw(0) 
						<< ".gt";
					std::string const renamed = ostr.str();
					oldgtnames.push_back(ostr.str());
					::libmaus::util::TempFileRemovalContainer::addTempFile(renamed);
					rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				// result.setGT(stringVectorAppend(GACR.gtpartnames,blockresults.getFiles().getGT()));
				
				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					remove ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			// rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus::timing::RealTimeClock mprtc;
			mprtc.start();
			result.setBWT(parallelGapFragMerge(
				bwtfilenames,
				stringVectorPack(gapfilenames),
				// result.getFiles().getBWT(),
				tmpfilenamebase+"_gpart",
				numthreads,
				lfblockmult,rlencoderblocksize));
			std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;

			for ( uint64_t i = 0; i < gapfilenames.size(); ++i )
				remove ( gapfilenames[i].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					remove ( bwtfilenames[i][j].c_str() );
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			RlToHwtBase<true>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			RlToHwtBase<false>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus::aio::CheckedOutputStream hwtreqCOS(result.getFiles().getHWTReq());
		RlToHwtTermRequest::serialise(
			hwtreqCOS,
			result.getFiles().getBWT(),
			result.getFiles().getHWT(),
			tmpfilenamebase + "_wt",
			huftreefilename,
			bwtterm,
			result.getBlockP0Rank(),
			input_types_type::utf8Wavelet()
		);
		hwtreqCOS.flush();
		hwtreqCOS.close();

		// remove obsolete files
		for ( uint64_t b = 0; b < mergereq.children.size(); ++b )
			mergereq.children[b]->sortresult.removeFiles();

		mergereq.releaseChildren();
	}

	static void mergeBlocks(
		MergeStrategyMergeInternalSmallBlock & mergereq,
		std::string const fn,
		uint64_t const fs,
		std::string const tmpfilenamebase,
		uint64_t const rlencoderblocksize,
		uint64_t const lfblockmult,
		uint64_t const numthreads,
		::std::map<int64_t,uint64_t> const & /* chist */,
		uint64_t const bwtterm,
		std::string const & huftreefilename
	)
	{
		assert ( mergereq.children.size() > 1 );
		assert ( mergereq.children.size() == mergereq.gaprequests.size()+1 );

		std::cerr << "[V] Merging BWT blocks MergeStrategyMergeInternalSmallBlock." << std::endl;

		/*
		 * remove unused file
		 */
		remove ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
		// fill result structure
		result.setBlockStart( mergereq.children[0]->sortresult.getBlockStart() );
		result.setCBlockSize( 0 );
		for ( uint64_t i = 0; i < mergereq.children.size(); ++i )
			result.setCBlockSize( result.getCBlockSize() + mergereq.children[i]->sortresult.getCBlockSize() );
		// set up temp file names
		// of output bwt,
		// sampled inverse suffix array filename,
		// gt bit array,
		// huffman shaped wavelet tree and
		// histogram
		result.setTempPrefixAndRegisterAsTemp(tmpfilenamebase + "_out",0 /* no preset bwt file names */, 0 /* no preset gt file names */);

		// if we merge only two blocks together, then we do not need to write the gap array to disk
		if ( mergereq.children.size() == 2 )
		{
			// std::cerr << "** WHITEBOX INTERNAL SMALL 1 **" << std::endl;
		
			std::string const & sblockhist = mergereq.children[1]->sortresult.getFiles().getHist();
			// load char histogram for last/second block (for merging)
			::libmaus::lf::DArray::unique_ptr_type accD(new ::libmaus::lf::DArray(
				sblockhist
				)
			);
			// first block
			::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = 
				mergereq.children[0]->sortresult;

			// start of first block
			uint64_t const blockstart = blockresults.getBlockStart();
			// size of first block
			uint64_t const cblocksize = blockresults.getCBlockSize();

			// compute gap array
			GapArrayByteComputationResult const GACR = computeGapArrayByte(
				fn,fs,*(mergereq.gaprequests[0]),
				mergereq.children[1]->sortresult.getFiles().getGT(), // previous gt files
				tmpfilenamebase + "_gparts", // new gt files
				tmpfilenamebase + "_gapoverflow",
				accD.get()
			);

			#if 0
			// concatenate gt vectors
			concatenateGT(
				GACR.gtpartnames,
				blockresults.getFiles().getGT(),
				result.getFiles().getGT()
			);
			#endif
			
			std::vector<std::string> oldgtnames;
			for ( uint64_t i = 0; i < blockresults.getFiles().getGT().size(); ++i )
			{
				std::ostringstream ostr;
				ostr << tmpfilenamebase << "_renamed_" << std::setw(6) << std::setfill('0') << i << std::setw(0) << ".gt";
				std::string const renamed = ostr.str();
				oldgtnames.push_back(ostr.str());
				::libmaus::util::TempFileRemovalContainer::addTempFile(renamed);
				rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
			}

			result.setGT(stringVectorAppend(GACR.gtpartnames,oldgtnames));

			// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
			libmaus::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap0dec(GACR.G->getDecoder());
			libmaus::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap0decbuf(new libmaus::suffixsort::GapArrayByteDecoderBuffer(*pgap0dec,8192));
			// libmaus::suffixsort::GapArrayByteDecoderBuffer::iterator gap0decbufit = pgap0decbuf->begin();
			result.setBlockP0Rank( mergeIsa(
				mergereq.children[1]->sortresult.getFiles().getSampledISA(), // old sampled isa
				blockresults.getFiles().getSampledISA(), // new sampled isa
				result.getFiles().getSampledISA(),blockstart,*pgap0decbuf/*gap0decbufit*//*GACR.G.begin()*/,cblocksize+1 /* GACR.G.size() */
			) );
			pgap0decbuf.reset();
			pgap0dec.reset();
			
			::libmaus::timing::RealTimeClock rtc; rtc.start();
			std::cerr << "[V] merging BWTs...";
			
			//
			uint64_t const totalsuf = result.getCBlockSize();
			// number of packets
			uint64_t const numpack = numthreads;
			// suffixes per thread
			uint64_t const tpacksize = (totalsuf + numpack-1)/numpack;
			//
			uint64_t ilow = 0;
			// intervals on G
			std::vector < std::pair<uint64_t,uint64_t> > wpacks;
			std::vector < std::string > encfilenames;
			// prefix sums over G
			std::vector < uint64_t > P;
			P.push_back(0);
			
			// std::cerr << "(computing work packets...";
			libmaus::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap1dec(GACR.G->getDecoder());
			libmaus::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap1decbuf(new libmaus::suffixsort::GapArrayByteDecoderBuffer(*pgap1dec,8192));
			libmaus::suffixsort::GapArrayByteDecoderBuffer::iterator gap1decbufit = pgap1decbuf->begin();
			::libmaus::timing::RealTimeClock wprtc; wprtc.start();
			while ( ilow != (cblocksize+1) )
			{
				uint64_t s = 0;
				uint64_t ihigh = ilow;
				
				while ( ihigh != (cblocksize+1) && s < tpacksize )
				{
					s += *(gap1decbufit++) + 1; // (GACR.G[ihigh++]+1);
					ihigh += 1;
				}
				
				uint64_t const p = s-(ihigh-ilow);
				
				if ( ihigh+1 == (cblocksize+1) && (*gap1decbufit == 0) /* GACR.G[ihigh] == 0 */ )
				{
					ihigh++;
					gap1decbufit++;
				}

				// std::cerr << "[" << ilow << "," << ihigh << ")" << std::endl;
				
				{
					libmaus::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap2dec(GACR.G->getDecoder(ilow));
					libmaus::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap2decbuf(new libmaus::suffixsort::GapArrayByteDecoderBuffer(*pgap2dec,8192));
					libmaus::suffixsort::GapArrayByteDecoderBuffer::iterator gap2decbufit = pgap2decbuf->begin();
					
					uint64_t a = 0;
					for ( uint64_t ia = ilow; ia < ihigh; ++ia )
						a += *(gap2decbufit++);
					
					assert ( p == a );
				}

				P.push_back(P.back() + p);
				wpacks.push_back(std::pair<uint64_t,uint64_t>(ilow,ihigh));
				encfilenames.push_back(
					tmpfilenamebase 
					// result.getFiles().getBWT() 
					+ "_"
					+ ::libmaus::util::NumberSerialisation::formatNumber(encfilenames.size(),6)
					+ ".bwt"
				);
				::libmaus::util::TempFileRemovalContainer::addTempFile(encfilenames.back());
				ilow = ihigh;
			}
			assert ( wpacks.size() <= numthreads );
			// std::cerr << "done,time=" << wprtc.getElapsedSeconds() << ")";
			pgap1decbuf.reset();
			pgap1dec.reset();
			
			// std::cerr << "(setting up IDDs...";
			wprtc.start();
			#if !defined(HUFRL)
			unsigned int const albits = rl_decoder::getAlBits(mergereq.children[0]->sortresult.getFiles().getBWT());
			#endif
			::libmaus::huffman::IndexDecoderDataArray IDD0(
				mergereq.children[0]->sortresult.getFiles().getBWT());
			::libmaus::huffman::IndexDecoderDataArray IDD1(
				mergereq.children[1]->sortresult.getFiles().getBWT());
			
			::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV0 = ::libmaus::huffman::IndexLoader::loadAccIndex(
				mergereq.children[0]->sortresult.getFiles().getBWT()
			);
			::libmaus::huffman::IndexEntryContainerVector::unique_ptr_type IECV1 = ::libmaus::huffman::IndexLoader::loadAccIndex(
				mergereq.children[1]->sortresult.getFiles().getBWT()
			);
						
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for ( int64_t b = 0; b < static_cast<int64_t>(wpacks.size()); ++b )
			{
				uint64_t const ilow = wpacks[b].first;
				uint64_t const ihigh = wpacks[b].second;
				
				if ( ilow != ihigh )
				{
					bool const islast = (ihigh == (cblocksize+1));
					std::string const encfilename = encfilenames[b];

					#if defined(HUFRL)
					rl_decoder leftrlin(IDD0,ilow);
					rl_decoder rightrlin(IDD1,P[b]);
					#else
					rl_decoder leftrlin(IDD0,IECV0.get(),ilow);
					rl_decoder rightrlin(IDD1,IECV1.get(),P[b]);
					#endif
					
					uint64_t const outsuf = (ihigh-ilow)-(islast?1:0) + (P[b+1]-P[b]);

					#if defined(HUFRL)
					::libmaus::huffman::RLEncoderStd bwtenc(encfilename,0     ,outsuf,rlencoderblocksize);
					#else
					::libmaus::gamma::GammaRLEncoder bwtenc(encfilename,albits,outsuf,rlencoderblocksize);
					#endif

					libmaus::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap3dec(GACR.G->getDecoder(ilow));
					libmaus::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap3decbuf(new libmaus::suffixsort::GapArrayByteDecoderBuffer(*pgap3dec,8192));
					libmaus::suffixsort::GapArrayByteDecoderBuffer::iterator gap3decbufit = pgap3decbuf->begin();
				
					if ( islast )
					{
						for ( uint64_t j = ilow; j < ihigh-1; ++j )
						{
							uint64_t const GACRGj = *(gap3decbufit++);
							
							for ( uint64_t i = 0; i < GACRGj; ++i )
								bwtenc.encode(rightrlin.decode());
								
							bwtenc.encode(leftrlin.decode());
						}
						
						assert ( ihigh == cblocksize+1 );

						uint64_t const GACRGcblocksize = *(gap3decbufit++);
						for ( uint64_t i = 0; i < GACRGcblocksize; ++i )
							bwtenc.encode(rightrlin.decode());						
					}
					else
					{
						for ( uint64_t j = ilow; j < ihigh; ++j )
						{
							uint64_t const GACRGj = *(gap3decbufit++);

							for ( uint64_t i = 0; i < GACRGj; ++i )
								bwtenc.encode(rightrlin.decode());
								
							bwtenc.encode(leftrlin.decode());
						}
					}
					
					bwtenc.flush();
				}
			}
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;

			#if 0
			std::cerr << "[V] concatenating bwt parts...";			
			rtc.start();
			#if defined(HUFRL)
			::libmaus::huffman::RLEncoderStd::concatenate(encfilenames,result.getFiles().getBWT());
			#else
			::libmaus::gamma::GammaRLEncoder::concatenate(encfilenames,result.getFiles().getBWT());
			#endif
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			result.setBWT(encfilenames);
			
			#if 0
			std::cerr << "[V] removing tmp files...";			
			rtc.start();
			for ( uint64_t i = 0; i < encfilenames.size(); ++i )
				remove ( encfilenames[i].c_str() );
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			// save histogram
			std::cerr << "[V] saving histogram...";
			rtc.start();
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
		}
		else
		{
			// std::cerr << "** WHITEBOX INTERNAL 2 **" << std::endl;
			
			std::vector < std::string > gapfilenames;
			std::vector < std::vector<std::string> > bwtfilenames;
			for ( uint64_t bb = 0; bb < mergereq.children.size(); ++bb )
			{
				// gap file name
				if ( bb+1 < mergereq.children.size() )
				{
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					::libmaus::util::TempFileRemovalContainer::addTempFile(newgapname);
					gapfilenames.push_back(newgapname);
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}

			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus::lf::DArray::unique_ptr_type accD(new ::libmaus::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
				std::string const newmergedgapoverflow = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".gapoverflow";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedgapoverflow);
				// gap file
				std::string const gapfile = gapfilenames[bx];

				// start of this block
				uint64_t const blockstart = blockresults.getBlockStart();
				// size of this block
				uint64_t const cblocksize = blockresults.getCBlockSize();

				// compute gap array
				GapArrayByteComputationResult const GACR = computeGapArrayByte(
					fn,fs,*(mergereq.gaprequests[bx]),
					mergedgtname,
					newmergedgtname,
					newmergedgapoverflow,
					accD.get()
				);

				// save the gap file			
				#if defined(HUFGAP)
				GACR.G->saveHufGapArray(gapfile);
				#else
				GACR.G->saveGammaGapArray(gapfile);
				#endif

				// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
				// libmaus::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
				libmaus::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap0dec(GACR.G->getDecoder());
				libmaus::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap0decbuf(new libmaus::suffixsort::GapArrayByteDecoderBuffer(*pgap0dec,8192));
				result.setBlockP0Rank( mergeIsa(mergedisaname,blockresults.getFiles().getSampledISA(),newmergedisaname,blockstart,*pgap0decbuf/*mergeGO*/ /*GACR.G.begin()*/,cblocksize+1 /*GACR.G.size()*/) );
				
				#if 0
				// concatenate gt vectors
				concatenateGT(GACR.gtpartnames,blockresults.getFiles().getGT(),newmergedgtname);
				#endif

				// rename files
				std::vector<std::string> oldgtnames;
				for ( uint64_t i = 0; i < blockresults.getFiles().getGT().size(); ++i )
				{
					std::ostringstream ostr;
					ostr << tmpfilenamebase 
						<< "_renamed_" 
						<< std::setw(6) << std::setfill('0') << bx << std::setw(0) 
						<< "_"
						<< std::setw(6) << std::setfill('0') << i << std::setw(0) 
						<< ".gt";
					std::string const renamed = ostr.str();
					oldgtnames.push_back(ostr.str());
					::libmaus::util::TempFileRemovalContainer::addTempFile(renamed);
					rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				// result.setGT(stringVectorAppend(GACR.gtpartnames,blockresults.getFiles().getGT()));
				
				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					remove ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			// rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus::timing::RealTimeClock mprtc;
			mprtc.start();
			result.setBWT(parallelGapFragMerge(
				bwtfilenames,
				stringVectorPack(gapfilenames),
				// result.getFiles().getBWT(),
				tmpfilenamebase+"_gpart",
				numthreads,
				lfblockmult,rlencoderblocksize));
			std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;

			for ( uint64_t i = 0; i < gapfilenames.size(); ++i )
				remove ( gapfilenames[i].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					remove ( bwtfilenames[i][j].c_str() );
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			RlToHwtBase<true>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			RlToHwtBase<false>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus::aio::CheckedOutputStream hwtreqCOS(result.getFiles().getHWTReq());
		RlToHwtTermRequest::serialise(
			hwtreqCOS,
			result.getFiles().getBWT(),
			result.getFiles().getHWT(),
			tmpfilenamebase + "_wt",
			huftreefilename,
			bwtterm,
			result.getBlockP0Rank(),
			input_types_type::utf8Wavelet()
		);
		hwtreqCOS.flush();
		hwtreqCOS.close();

		// remove obsolete files
		for ( uint64_t b = 0; b < mergereq.children.size(); ++b )
			mergereq.children[b]->sortresult.removeFiles();

		mergereq.releaseChildren();
	}

	static void mergeBlocks(
		MergeStrategyMergeExternalBlock & mergereq,
		std::string const fn,
		uint64_t const fs,
		std::string const tmpfilenamebase,
		uint64_t const rlencoderblocksize,
		uint64_t const lfblockmult,
		uint64_t const numthreads,
		::std::map<int64_t,uint64_t> const & /* chist */,
		uint64_t const bwtterm,
		uint64_t const mem,
		std::string const & huftreefilename
	)
	{
		std::cerr << "[V] Merging BWT blocks MergeStrategyMergeExternalBlock." << std::endl;
		
		assert ( mergereq.children.size() > 1 );
		assert ( mergereq.children.size() == mergereq.gaprequests.size()+1 );

		/*
		 * remove unused file
		 */
		remove ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
		// fill result structure
		result.setBlockStart( mergereq.children[0]->sortresult.getBlockStart() );
		result.setCBlockSize ( 0 );
		for ( uint64_t i = 0; i < mergereq.children.size(); ++i )
			result.setCBlockSize( result.getCBlockSize() + mergereq.children[i]->sortresult.getCBlockSize() );
		// set up
		// filenames of output bwt,
		// sampled inverse suffix array filename,
		// gt bit array,
		// huffman shaped wavelet tree and
		// histogram		
		result.setTempPrefixAndRegisterAsTemp(tmpfilenamebase + "_out",0,0);

		{
			std::vector < std::string > gapfilenameprefixes;
			std::vector < std::vector < std::string > > gapfilenames;
			std::vector < std::vector < std::string > > bwtfilenames;
			for ( uint64_t bb = 0; bb < mergereq.children.size(); ++bb )
			{
				// gap file name
				if ( bb+1 < mergereq.children.size() )
				{
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					gapfilenameprefixes.push_back(newgapname);	
					gapfilenames.push_back(std::vector<std::string>());
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}


			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus::lf::DArray::unique_ptr_type accD(new ::libmaus::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// std::cerr << "** WHITEBOX EXTERNAL **" << std::endl;
				
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
				// gap file
				std::string const gapfilenameprefix = gapfilenameprefixes[bx];

				// start of this block
				uint64_t const blockstart = blockresults.getBlockStart();
				// size of this block
				uint64_t const cblocksize = blockresults.getCBlockSize();

				// std::cerr << "*** compute sparse ***" << std::endl;

				SparseGapArrayComputationResult const GACR = computeSparseGapArray(
					fn,fs,*(mergereq.gaprequests[bx]),
					mergereq.children[mergereq.gaprequests[bx]->into]->getIHWTSpaceBytes(),
					mergedgtname,newmergedgtname,accD.get(),
					gapfilenameprefix,tmpfilenamebase+"_sparsegap",
					mem
				);
				gapfilenames[bx] = GACR.fn;

				// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
				libmaus::gamma::GammaGapDecoder GGD(gapfilenames[bx]);
				result.setBlockP0Rank( mergeIsa(mergedisaname,blockresults.getFiles().getSampledISA(),newmergedisaname,blockstart,GGD,cblocksize+1 /*GACR.G.size()*/) );
				
				// concatenate gt vectors
				//concatenateGT(GACR.gtpartnames,blockresults.getFiles().getGT(),newmergedgtname);

				// rename files
				std::vector<std::string> oldgtnames;
				for ( uint64_t i = 0; i < blockresults.getFiles().getGT().size(); ++i )
				{
					std::ostringstream ostr;
					ostr << tmpfilenamebase 
						<< "_renamed_" 
						<< std::setw(6) << std::setfill('0') << bx << std::setw(0) 
						<< "_"
						<< std::setw(6) << std::setfill('0') << i << std::setw(0) 
						<< ".gt";
					std::string const renamed = ostr.str();
					oldgtnames.push_back(ostr.str());
					::libmaus::util::TempFileRemovalContainer::addTempFile(renamed);
					rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					remove ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			//rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus::timing::RealTimeClock mprtc;
			mprtc.start();
			result.setBWT(
				parallelGapFragMerge(bwtfilenames,gapfilenames,tmpfilenamebase+"_gpart",numthreads,
					lfblockmult,rlencoderblocksize)
			);
			std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;

			for ( uint64_t i = 0; i < gapfilenames.size(); ++i )
				for ( uint64_t j = 0; j < gapfilenames[i].size(); ++j )
					remove ( gapfilenames[i][j].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					remove ( bwtfilenames[i][j].c_str() );		
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			RlToHwtBase<true>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			RlToHwtBase<false>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus::aio::CheckedOutputStream hwtreqCOS(result.getFiles().getHWTReq());
		RlToHwtTermRequest::serialise(hwtreqCOS,
			result.getFiles().getBWT(),
			result.getFiles().getHWT(),
			tmpfilenamebase + "_wt",
			huftreefilename,
			bwtterm,
			result.getBlockP0Rank(),
			input_types_type::utf8Wavelet()
		);
		hwtreqCOS.flush();
		hwtreqCOS.close();

		// remove obsolete files
		for ( uint64_t b = 0; b < mergereq.children.size(); ++b )
			mergereq.children[b]->sortresult.removeFiles();

		mergereq.releaseChildren();
	}
	
	static void sortIsaFile(std::string const & mergedisaname, uint64_t const blockmem)
	{
		// sort sampled inverse suffix array file
		std::string const mergeisatmp = mergedisaname+".tmp";
		::libmaus::util::TempFileRemovalContainer::addTempFile(mergeisatmp);
		std::string const mergeisatmpout = mergedisaname+".tmp.out";
		::libmaus::util::TempFileRemovalContainer::addTempFile(mergeisatmpout);
		// uint64_t const blockmem = 5*blocksize;
		// uint64_t const blockels = (blockmem + 2*sizeof(uint64_t)-1)/(2*sizeof(uint64_t));
		::libmaus::sorting::PairFileSorting::sortPairFile(
			std::vector<std::string>(1,mergedisaname),mergeisatmp,true /* second comp */,
			true,true,mergeisatmpout,blockmem/2/*par*/,true /* parallel */);
		remove ( (mergeisatmp).c_str() );
		rename ( mergeisatmpout.c_str(), mergedisaname.c_str() );
	
	}

	static uint64_t readBlockRanksSize(std::string const & mergedisaname)
	{
		return ::libmaus::util::GetFileSize::getFileSize(mergedisaname)/(2*sizeof(uint64_t));		
	}

	static ::libmaus::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > readBlockRanks(std::string const & mergedisaname)
	{
		// read sampled isa
		uint64_t const nsisa = readBlockRanksSize(mergedisaname);
		::libmaus::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > blockranks(nsisa,false);
		::libmaus::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type SGIsisa(new ::libmaus::aio::SynchronousGenericInput<uint64_t>(mergedisaname,16*1024));
		for ( uint64_t i = 0; i < nsisa; ++i )
		{
			int64_t const r = SGIsisa->get();
			assert ( r >= 0 );
			int64_t const p = SGIsisa->get();
			assert ( p >= 0 );
			blockranks [ i ] = std::pair<uint64_t,uint64_t>(r,p);
		}
		SGIsisa.reset();
		
		return blockranks;
	}
	
	static void checkSampledSA(
		std::string const & fn,
		uint64_t const fs,
		std::string const mergedsaname,
		uint64_t const sasamplingrate,
		uint64_t const numthreads,
		uint64_t const lfblockmult
	)
	{
		::libmaus::parallel::OMPLock cerrlock;
		// number of sampled suffix array elements
		uint64_t const nsa = (fs + sasamplingrate - 1) / sasamplingrate;
		
		// check that this matches what we have in the file
		assert ( ::libmaus::util::GetFileSize::getFileSize(mergedsaname) / (sizeof(uint64_t)) ==  nsa + 2 );
		
		if ( nsa && nsa-1 )
		{
			uint64_t const checkpos = nsa-1;
			uint64_t const satcheckpacks = numthreads * lfblockmult;
			uint64_t const sacheckblocksize = (checkpos + satcheckpacks-1) / satcheckpacks;
			uint64_t const sacheckpacks = ( checkpos + sacheckblocksize - 1 ) / sacheckblocksize;
			
			std::cerr << "[V] checking suffix array on text...";
			::libmaus::parallel::SynchronousCounter<uint64_t> SC;
			int64_t lastperc = -1;
			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for ( int64_t t = 0; t < static_cast<int64_t>(sacheckpacks); ++t )
			{
				uint64_t const low = t * sacheckblocksize;
				uint64_t const high = std::min(low+sacheckblocksize,checkpos);
				uint64_t const cnt = high-low;
				
				// std::cerr << "low=" << low << " high=" << high << " nsa=" << nsa << " cnt=" << cnt << std::endl;

				::libmaus::aio::SynchronousGenericInput<uint64_t> SGIsa(mergedsaname,16*1024,low+2,cnt+1);
				
				typename input_types_type::circular_suffix_comparator CSC(fn,fs);
				
				int64_t const fp = SGIsa.get();
				assert ( fp >= 0 );
				
				uint64_t prev = fp;
				
				while ( SGIsa.peek() >= 0 )
				{
					int64_t const p = SGIsa.get();
					assert ( p >= 0 );

					// std::cerr << "pre " << prev.first << " SA[" << r << "]=" << p << std::endl;
					
					bool const ok = CSC (prev,p);
					assert ( ok );
					
					prev = p;
				}
				
				uint64_t const sc = ++SC;
				int64_t const newperc = (sc*100) / sacheckpacks;
				
				cerrlock.lock();
				if ( newperc != lastperc )
				{
					std::cerr << "(" << newperc << ")";
					lastperc = newperc;
				}
				cerrlock.unlock();
			}
			std::cerr << "done." << std::endl;
		}
	
	}
	
	static void computeSampledSA(
		std::string const & fn,
		uint64_t const fs,
		::libmaus::lf::ImpCompactHuffmanWaveletLF const & IHWT,
		std::string const & mergedisaname,
		std::string const & outfn,
		std::string const & tmpfilenamebase,
		uint64_t const numthreads,
		uint64_t const lfblockmult,
		uint64_t const sasamplingrate,
		uint64_t const isasamplingrate,
		uint64_t const blockmem
	)
	{
		// read size of sampled isa
		uint64_t const nsisa = readBlockRanksSize(mergedisaname);

		::libmaus::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type SGIsisasa(new ::libmaus::aio::SynchronousGenericInput<uint64_t>(mergedisaname,16*1024));
		int64_t const fr = SGIsisasa->get(); assert ( fr != -1 );
		int64_t const fp = SGIsisasa->get(); assert ( fp != -1 );
		
		std::cerr << "[V] computing sampled suffix array parts...";
	
		std::pair<uint64_t,uint64_t> const isa0(fr,fp);
		std::pair<uint64_t,uint64_t> isapre(isa0);
		
		::std::vector< std::string > satempfilenames(numthreads);
		::std::vector< std::string > isatempfilenames(numthreads);
		::libmaus::autoarray::AutoArray < ::libmaus::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type > SAF(numthreads);
		::libmaus::autoarray::AutoArray < ::libmaus::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type > ISAF(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			satempfilenames[i] = ( tmpfilenamebase + ".sampledsa_" + ::libmaus::util::NumberSerialisation::formatNumber(i,6) );
			::libmaus::util::TempFileRemovalContainer::addTempFile(satempfilenames[i]);
			SAF[i] = UNIQUE_PTR_MOVE(::libmaus::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type(
				new ::libmaus::aio::SynchronousGenericOutput<uint64_t>(satempfilenames[i],8*1024)));

			isatempfilenames[i] = ( tmpfilenamebase + ".sampledisa_" + ::libmaus::util::NumberSerialisation::formatNumber(i,6) );
			::libmaus::util::TempFileRemovalContainer::addTempFile(isatempfilenames[i]);
			ISAF[i] = UNIQUE_PTR_MOVE(::libmaus::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type(
				new ::libmaus::aio::SynchronousGenericOutput<uint64_t>(isatempfilenames[i],8*1024)));
		}
		
		std::vector < std::pair< std::pair<uint64_t,uint64_t>, std::pair<uint64_t,uint64_t> > > WV;
		int64_t lastperc = -1;

		if ( nsisa > 1 )
		{
			for ( int64_t i = 1; i <= static_cast<int64_t>(nsisa); ++i )
			{
				int64_t const nr = (i == static_cast<int64_t>(nsisa)) ? isa0.first : SGIsisasa->get();
				int64_t const np = (i == static_cast<int64_t>(nsisa)) ? isa0.second : ((nr != -1) ? SGIsisasa->get() : -1);
				assert ( np >= 0 );
				
				std::pair<uint64_t,uint64_t> isai(nr,np);
				
				WV.push_back(std::pair< std::pair<uint64_t,uint64_t>, std::pair<uint64_t,uint64_t> >(isai,isapre));
				
				isapre.first = nr;
				isapre.second = np;

				if ( ((WV.size() % (lfblockmult*numthreads)) == 0) || i == static_cast<int64_t>(nsisa) )
				{
					#if defined(_OPENMP)
					#pragma omp parallel for
					#endif
					for ( int64_t j = 0; j < static_cast<int64_t>(WV.size()); ++j )
					{
						#if defined(_OPENMP)
						uint64_t const tid = omp_get_thread_num();
						#else
						uint64_t const tid = 0;
						#endif
						checkBwtBlockDecode(WV[j].first,WV[j].second,fn,fs,IHWT,*SAF[tid],*ISAF[tid],sasamplingrate,isasamplingrate);
					}
					
					WV.resize(0);
				}
				
				int64_t const newperc = ((i)*100) / (nsisa);
				if ( newperc != lastperc )
				{
					std::cerr << "(" << newperc << ")";
					lastperc = newperc;
				}
			}
		}
		else
		{
			assert ( fp == 0 );
			checkBwtBlockDecode(
				std::pair<uint64_t,uint64_t>(fr,0),
				std::pair<uint64_t,uint64_t>(fr,0),
				fn,fs,IHWT,*SAF[0],*ISAF[0],sasamplingrate,isasamplingrate,
				fs);
		}
		
		SGIsisasa.reset();
		
		for ( uint64_t i = 0; i < SAF.size(); ++i )
		{
			SAF[i]->flush();
			SAF[i].reset();
		}
		for ( uint64_t i = 0; i < ISAF.size(); ++i )
		{
			ISAF[i]->flush();
			ISAF[i].reset();
		}
			
		std::cerr << "done." << std::endl;
		
		std::cerr << "[V] sorting and merging sampled suffix array parts...";
		std::string const mergedsaname = ::libmaus::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".sa";
		{
		::libmaus::aio::CheckedOutputStream::unique_ptr_type pmergedsa(new ::libmaus::aio::CheckedOutputStream(mergedsaname));
		// write sampling rate
		::libmaus::serialize::Serialize<uint64_t>::serialize(*pmergedsa,sasamplingrate);
		::libmaus::serialize::Serialize<uint64_t>::serialize(*pmergedsa,(fs + sasamplingrate-1)/sasamplingrate);
		std::string const mergesatmp = mergedsaname + ".tmp";
		::libmaus::util::TempFileRemovalContainer::addTempFile(mergesatmp);
		::libmaus::sorting::PairFileSorting::sortPairFile(
			satempfilenames,mergesatmp,
			false /* second comp */,
			false /* keep first */,
			true /* keep second */,
			*pmergedsa /* output stream */,
			blockmem/2/*par*/,
			true /* parallel */,
			true /* delete input */
		);
		pmergedsa->flush();
		pmergedsa.reset();
		remove(mergesatmp.c_str());
		}
		std::cerr << "done." << std::endl;		

		std::cerr << "[V] sorting and merging sampled inverse suffix array parts...";
		std::string const mergedisaoutname = ::libmaus::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".isa";
		::libmaus::aio::CheckedOutputStream::unique_ptr_type pmergedisa(new ::libmaus::aio::CheckedOutputStream(mergedisaoutname));
		// write sampling rate
		::libmaus::serialize::Serialize<uint64_t>::serialize(*pmergedisa,isasamplingrate);
		::libmaus::serialize::Serialize<uint64_t>::serialize(*pmergedisa,(fs+isasamplingrate-1)/isasamplingrate);
		std::string const mergeisatmp = mergedisaoutname + ".tmp";
		::libmaus::util::TempFileRemovalContainer::addTempFile(mergeisatmp);
		::libmaus::sorting::PairFileSorting::sortPairFile(
			isatempfilenames,mergeisatmp,
			false /* second comp */,
			false /* keep first */,
			true /* keep second */,
			*pmergedisa /* output stream */,
			blockmem/2/*par*/,
			true /* parallel */,
			true /* delete input */
		);
		std::cerr << "done." << std::endl;		
		remove(mergeisatmp.c_str());

		#if 0
		// check sampled suffix array by pairwise comparison on text
		// this can take quadratic time in fs
		checkSampledSA(fn,fs,mergedsaname,sasamplingrate,numthreads,lfblockmult);
		#endif
	}

	static uint64_t getDefaultSaSamplingRate()
	{
		return 32;
	}
	
	static uint64_t getDefaultIsaSamplingRate()
	{
		return 256*1024;
	}

	static uint64_t getDefaultBlockSize(uint64_t const mem, uint64_t const threads)
	{
		return std::max(static_cast<uint64_t>(0.95 * mem / ( 5 * threads )),static_cast<uint64_t>(1));
	}

	static uint64_t getDefaultBlockAlign()
	{
		return 16*1024ull;
	}

	static uint64_t getDefaultWordsPerThread()
	{
		return 64ull*1024ull;
	}

	static uint64_t getDefaultNumThreads()
	{
		#if defined(_OPENMP)
		return omp_get_max_threads();
		#else
		return 1;
		#endif
	}

	static uint64_t getDefaultMem()
	{
		return 2ull * 1024ull * 1024ull * 1024ull;
	}
	
	static std::map<int64_t,uint64_t> mergeMaps(
		std::map<int64_t,uint64_t> const & A,
		std::map<int64_t,uint64_t> const & B)
	{
		std::map<int64_t,uint64_t>::const_iterator aita = A.begin(), aite = A.end();
		std::map<int64_t,uint64_t>::const_iterator bita = B.begin(), bite = B.end();
		std::map<int64_t,uint64_t> C;
		
		while ( aita != aite && bita != bite )
		{
			if ( aita->first < bita->first )
			{
				C[aita->first] = aita->second;
				++aita;
			}
			else if ( bita->first < aita->first )
			{
				C[bita->first] = bita->second;
				++bita;
			}
			else
			{
				C[aita->first] = aita->second + bita->second;
				++aita; ++bita;
			}
		}
		
		while ( aita != aite )
		{
			C[aita->first] = aita->second;
			++aita;			
		}
		
		while ( bita != bite )
		{
			C[bita->first] = bita->second;
			++bita;		
		}
		
		return C;
	}
	
	struct HashHistogram
	{
		libmaus::autoarray::AutoArray<uint64_t> L;
		libmaus::util::ExtendingSimpleCountingHash<uint64_t,uint64_t> H;
		libmaus::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > P;
		uint64_t p;
		
		HashHistogram(uint64_t const lowsize = 256, uint64_t const bigsize = (1ull<<16) ) : L(lowsize), H(libmaus::math::nextTwoPow(bigsize)), P(64ull*1024ull), p(0) {}
		
		void clear()
		{
			std::fill(L.begin(),L.end(),0ull);
			H.clear();
			p = 0;
		}
		
		void checkSize()
		{
			H.checkExtend();
		}
		
		void increment(uint64_t const k)
		{
			if ( k < L.size() )
				__sync_fetch_and_add(L.begin()+k,1);
			else
				H.insert(k);
		}
		
		void ensureSize(uint64_t const s)
		{
			while ( H.getTableSize() < s )
				H.extendInternal();
		}
		
		void extract()
		{
			uint64_t nonzero = 0;
			for ( uint64_t i = 0; i < L.size(); ++i )
				if ( L[i] )
					nonzero++;
			for ( libmaus::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::key_type const * k = H.begin(); k != H.end(); ++k )
				if ( *k != libmaus::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::unused() )
					nonzero++;
					
			if ( P.size() < nonzero )
				P.resize(nonzero);
				
			p = 0;

			for ( uint64_t i = 0; i < L.size(); ++i )
				if ( L[i] )
					P[p++] = std::pair<uint64_t,uint64_t>(i,L[i]);
					
			for ( libmaus::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::key_type const * k = H.begin(); k != H.end(); ++k )
				if ( *k != libmaus::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::unused() )
					P[p++] = std::pair<uint64_t,uint64_t>(*k,H.cntbegin() [ (k-H.begin()) ] );
					
			std::sort(P.begin(),P.begin()+p);
		}
	};

	static void getBlockSymFreqsHash(std::string const fn, uint64_t const glow, uint64_t ghigh, HashHistogram & H)
	{
		typedef typename input_types_type::linear_wrapper stream_type;
		typedef typename input_types_type::base_input_stream::char_type char_type;
		typedef typename ::libmaus::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
		uint64_t const fs = ghigh-glow;
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif

		uint64_t const loopsize = 16ull*1024ull;
		uint64_t const elperthread = 8*loopsize;
		H.ensureSize(elperthread * numthreads);
		H.clear();
		
		uint64_t const symsperfrag = (fs + numthreads - 1)/numthreads;
		// uint64_t const numfrags = (fs + symsperfrag - 1)/symsperfrag;
		uint64_t const loops = (symsperfrag + loopsize -1)/loopsize;

		::libmaus::autoarray::AutoArray<unsigned_char_type> GB(loopsize*numthreads,false);

		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif		
		for ( int64_t t = 0; t < static_cast<int64_t>(numthreads); ++t )
		{
			uint64_t const low  = std::min(glow + t * symsperfrag,ghigh);
			uint64_t const high = std::min(low+symsperfrag       ,ghigh);
			uint64_t const size = high-low;
			uint64_t todo = size;
			
			unsigned_char_type * const B = GB.begin() + t * loopsize;

			stream_type CIS(fn);
			CIS.seekg(low);
			
			for ( uint64_t i = 0; i < loops; ++i )
			{
				uint64_t const toread = std::min(todo,loopsize);
				
				if ( toread )
				{
					CIS.read(reinterpret_cast<char_type *>(B),toread);
					assert ( CIS.gcount() == static_cast<int64_t>(toread) );
				}

				for ( uint64_t i = 0; i < toread; ++i )
					H.increment(B[i]);

				todo -= toread;

				H.checkSize();
			}			
		}
		
		H.extract();
	}
	
	static std::map<int64_t,uint64_t> getBlockSymFreqs(std::string const fn, uint64_t const glow, uint64_t ghigh)
	{
		typedef typename input_types_type::linear_wrapper stream_type;
		typedef typename input_types_type::base_input_stream::char_type char_type;
		typedef typename ::libmaus::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
		uint64_t const fs = ghigh-glow;
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		uint64_t const symsperfrag = (fs + numthreads - 1)/numthreads;
		uint64_t const numfrags = (fs + symsperfrag - 1)/symsperfrag;

		libmaus::util::HistogramSet HS(numfrags,256);

		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif		
		for ( int64_t t = 0; t < static_cast<int64_t>(numfrags); ++t )
		{
			uint64_t const low = glow + t * symsperfrag;
			uint64_t const high = std::min(low+symsperfrag,ghigh);
			uint64_t const size = high-low;
			uint64_t todo = size;

			stream_type CIS(fn);
			CIS.seekg(low);

			::libmaus::autoarray::AutoArray<unsigned_char_type> B(16*1024,false);
			libmaus::util::Histogram & H = HS[t];
			
			while ( todo )
			{
				uint64_t const toread = std::min(todo,B.size());
				CIS.read(reinterpret_cast<char_type *>(B.begin()),toread);
				assert ( CIS.gcount() == static_cast<int64_t>(toread) );
				for ( uint64_t i = 0; i < toread; ++i )
					H (B[i]);
				todo -= toread;
			}
		}

		::libmaus::util::Histogram::unique_ptr_type PH(UNIQUE_PTR_MOVE(HS.merge()));
		::libmaus::util::Histogram & H(*PH);

		return H.getByType<int64_t>();
	}
				
	static int computeBwt(::libmaus::util::ArgInfo const & arginfo)
	{
		uint64_t mcnt = 0;
		::libmaus::util::TempFileRemovalContainer::setup();
		uint64_t const rlencoderblocksize = 16*1024;
		::libmaus::parallel::OMPLock cerrlock;

		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else	
		uint64_t const numthreads = 1;
		#endif

		// total memory available
		uint64_t const mem = std::max(static_cast<uint64_t>(1),arginfo.getValueUnsignedNumeric<uint64_t>("mem",getDefaultMem()));
		// base for tmp file names
		std::string const tmpfilenamebase = arginfo.getDefaultTmpFileName();
		// file name of serialised character histogram
		std::string const chistfilename = tmpfilenamebase + ".chist";
		// file name of serialised huffman tree
		std::string const huftreefilename = tmpfilenamebase + ".huftree";
		::libmaus::util::TempFileRemovalContainer::addTempFile(chistfilename);
		::libmaus::util::TempFileRemovalContainer::addTempFile(huftreefilename);
		// output file name
		std::string const outfn = arginfo.getValue<std::string>("outputfilename",tmpfilenamebase+".bwt");
		// final inverse suffix array sampling rate
		uint64_t const isasamplingrate = ::libmaus::math::nextTwoPow(arginfo.getValue<uint64_t>("isasamplingrate",getDefaultIsaSamplingRate()));
		// final suffix array sampling rate
		uint64_t const sasamplingrate = ::libmaus::math::nextTwoPow(arginfo.getValue<uint64_t>("sasamplingrate",getDefaultSaSamplingRate()));
		// ISA sampling rate during block merging
		uint64_t const preisasamplingrate = ::libmaus::math::nextTwoPow(arginfo.getValue<uint64_t>("preisasamplingrate",256*1024));
		// target block alignment
		uint64_t const tblockalign = std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("blockalign",getDefaultBlockAlign()));
		// block alignment, a multiple of preisasamplingrate
		uint64_t const blockalign = ((tblockalign + preisasamplingrate - 1)/preisasamplingrate)*preisasamplingrate;
		// target block size
		uint64_t const tblocksize = std::max(static_cast<uint64_t>(1),arginfo.getValueUnsignedNumeric<uint64_t>("blocksize",getDefaultBlockSize(mem,numthreads)));
		std::string const fn = arginfo.getRestArg<std::string>(0);
		/* check if input file exists */
		if ( ! ::libmaus::util::GetFileSize::fileExists(fn) )
		{
			::libmaus::exception::LibMausException se;
			se.getStream() << "File " << fn << " does not exit." << std::endl;
			se.finish();
			throw se;
		}
		/* get file size */
		uint64_t const fs = input_types_type::linear_wrapper::getFileSize(fn);
		
		/* check that file is not empty */
		if ( ! fs )
		{
			::libmaus::exception::LibMausException se;
			se.getStream() << "File " << fn << " is empty." << std::endl;
			se.finish();
			throw se;		
		}
		
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		
		// target number of blocks
		uint64_t const tnumblocks = ( fs + tblocksize - 1 ) / tblocksize;
		// actual block size (multiple of blockalign)
		uint64_t const blocksize = (( (( fs + tnumblocks - 1 ) / tnumblocks) + blockalign - 1 ) / blockalign ) * blockalign;
		// actual number of blocks
		uint64_t const numblocks = ( fs + blocksize - 1 ) / blocksize;
		// maximum block size
		// uint64_t const maxblocksize = (numblocks > 1) ? blocksize : fs;

		std::cerr << "[V] sorting file " << fn << " of size " << fs << " with block size " << blocksize << std::endl;
		
		// there should be at least one block
		assert ( numblocks );

		::libmaus::suffixsort::BwtMergeTempFileNameSetVector blocktmpnames(tmpfilenamebase, numblocks, numthreads /* bwt */, numthreads /* gt */);

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		std::cerr << "[V] processing file " << fn << " of size " << fs << " block size " << blocksize 
			<< " number of blocks " << numblocks << std::endl;
		
		std::cerr << "[V] computing symbol frequences" << std::endl;
		std::map<int64_t,uint64_t> chistnoterm;	
		// std::vector< std::map<int64_t,uint64_t> > blockfreqvec(numblocks);
		libmaus::timing::RealTimeClock rtc;
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			// block id
			uint64_t const b = numblocks-bb-1;
			// start of block in file
			uint64_t const blockstart = b*blocksize;
			// size of this block
			uint64_t const cblocksize = std::min(blocksize,fs-blockstart);

			std::map<int64_t,uint64_t> const blockfreqs =
				getBlockSymFreqs(fn,blockstart,blockstart+cblocksize);
			
			std::string const freqstmpfilename = blocktmpnames[b].getHist() + ".freqs";
			libmaus::util::TempFileRemovalContainer::addTempFile(freqstmpfilename);	
			libmaus::aio::CheckedOutputStream freqCOS(freqstmpfilename);
			libmaus::util::NumberMapSerialisation::serialiseMap(freqCOS,blockfreqs);
			freqCOS.flush();
			freqCOS.close();
			
			#if 0
			rtc.start();
			std::map<int64_t,uint64_t> const blockfreqs =
				getBlockSymFreqs(fn,blockstart,blockstart+cblocksize);
			double const t0 = rtc.getElapsedSeconds();

			HashHistogram HH;
			rtc.start();
			getBlockSymFreqsHash(fn,blockstart,blockstart+cblocksize,HH);
			double const t1 = rtc.getElapsedSeconds();
			
			std::cerr << "(" << blockfreqs.size() << "," << HH.p << ")";
			std::cerr << "[" << t0 << "," << t1 << "]";
			#endif

			chistnoterm = mergeMaps(chistnoterm,blockfreqs);

			// blockfreqvec[b] = blockfreqs;
		}
		
		#if 0
		exit(0);
		#endif

		/* get symbol count histogram */
		int64_t const bwtterm = chistnoterm.size() ? (chistnoterm.rbegin()->first+1) : 0;
		::std::map<int64_t,uint64_t> chist = chistnoterm;
		chist[bwtterm] = 1;
		uint64_t const maxsym = chist.size() ? chist.rbegin()->first : 0;

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		libmaus::aio::CheckedOutputStream::unique_ptr_type chistCOS(new libmaus::aio::CheckedOutputStream(chistfilename));
		(*chistCOS) << ::libmaus::util::NumberMapSerialisation::serialiseMap(chist);
		chistCOS->flush();
		chistCOS->close();
		chistCOS.reset();
		
		std::cerr << "[V] computed symbol frequences, input alphabet size is " << chistnoterm.size() << std::endl;

		std::cerr << "[V] bwtterm=" << bwtterm << std::endl;

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		libmaus::huffman::HuffmanTree::unique_ptr_type uhnode(new libmaus::huffman::HuffmanTree(chist.begin(),chist.size(),false,true,true));
		
		libmaus::aio::CheckedOutputStream::unique_ptr_type huftreeCOS(new libmaus::aio::CheckedOutputStream(huftreefilename));
		uhnode->serialise(*huftreeCOS);
		huftreeCOS->flush();
		huftreeCOS->close();
		huftreeCOS.reset();

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		::libmaus::huffman::HuffmanTree::EncodeTable::unique_ptr_type EC(
			new ::libmaus::huffman::HuffmanTree::EncodeTable(*uhnode));

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		std::vector < MergeStrategyBlock::shared_ptr_type > stratleafs(numblocks);

		bool const computeTermSymbolHwt = arginfo.getValue<int>("computeTermSymbolHwt",false);;
	
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			uint64_t const b = numblocks-bb-1;
			// start of block in file
			uint64_t const blockstart = b*blocksize;
			// size of this block
			uint64_t const cblocksize = std::min(blocksize,fs-blockstart);
			// symbol frequency map
			// std::map<int64_t,uint64_t> const & blockfreqs = blockfreqvec[b];
			std::string const freqstmpfilename = blocktmpnames[b].getHist() + ".freqs";
			libmaus::aio::CheckedInputStream freqCIS(freqstmpfilename);
			std::map<int64_t,uint64_t> const blockfreqs = 
				libmaus::util::NumberMapSerialisation::deserialiseMap<libmaus::aio::CheckedInputStream,int64_t,uint64_t>(freqCIS);
			freqCIS.close();
			remove(freqstmpfilename.c_str());
			// 
			uint64_t const sourcelengthbits = input_types_type::getSourceLengthBits(fn,blockstart,blockstart+cblocksize,blockfreqs);
			//
			uint64_t const sourcelengthbytes = input_types_type::getSourceLengthBytes(fn,blockstart,blockstart+cblocksize,blockfreqs);
			//
			uint64_t const sourcetextindexbits = input_types_type::getSourceTextIndexBits(fn,blockstart,blockstart+cblocksize,blockfreqs);

			MergeStrategyBlock::shared_ptr_type PMSB = MergeStrategyBaseBlock::construct(
				blockstart,blockstart+cblocksize,
				*EC,
				blockfreqs,
				sourcelengthbits,sourcelengthbytes,
				sourcetextindexbits
			);
			
			/* set up and register sort request */
			::libmaus::suffixsort::BwtMergeZBlockRequestVector zreqvec;
			dynamic_cast<MergeStrategyBaseBlock *>(PMSB.get())->sortreq = 
				BwtMergeBlockSortRequest(
					input_types_type::getType(),
					fn,fs,
					chistfilename,
					huftreefilename,
					bwtterm,
					maxsym,
					blocktmpnames[b].serialise(),
					tmpfilenamebase,
					rlencoderblocksize,
					preisasamplingrate,
					blockstart,cblocksize,
					zreqvec,
					computeTermSymbolHwt
				);

			// std::cerr << *PMSB;
			
			stratleafs[b] = PMSB;

			#if defined(FERAMANZGEN_MEMORY_DEBUG)
			std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
			#endif
		}

		uhnode.reset();
		EC.reset();

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif


		// word available per computation thread
		uint64_t const wordsperthread = std::max(static_cast<uint64_t>(1),arginfo.getValueUnsignedNumeric<uint64_t>("wordsperthread",getDefaultWordsPerThread()));

		std::cerr << "[V] constructing merge tree" << std::endl;

		// construct merge tree and register z requests
		MergeStrategyBlock::shared_ptr_type mergetree = constructMergeTree(stratleafs,mem,numthreads,wordsperthread);
		
		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		// inner node queue
		std::deque<MergeStrategyBlock *> itodo;

		std::cerr << "[V] sorting blocks" << std::endl;

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		// sort single blocks
		BaseBlockSorting::unique_ptr_type BBS(new BaseBlockSorting(stratleafs,mem,numthreads,itodo));
		BBS->start();
		BBS->join();
		BBS.reset();

		std::cerr << "[V] sorted blocks" << std::endl;
		
		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		#if 0
		uint64_t maxhwtsize = 0;
		for ( uint64_t i = 0; i < stratleafs.size(); ++i )
			maxhwtsize = std::max(maxhwtsize,::libmaus::util::GetFileSize::getFileSize(stratleafs[i]->sortresult.getFiles().getHWT()));
		#endif
		
		std::cerr << "[V] filling gap request objects" << std::endl;
		mergetree->fillGapRequestObjects(numthreads);

		#if defined(FERAMANZGEN_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		std::cerr << "[V]" << std::string(80,'-') << std::endl;
		std::cerr << *mergetree;
		std::cerr << "[V]" << std::string(80,'-') << std::endl;
		
		uint64_t mtmpid = 0;
		uint64_t const lfblockmult = 1; // std::max(static_cast<uint64_t>(1),arginfo.getValue<uint64_t>("lfblockmult", (maxblocksize < (1024*1024)) ? 1 : 4));
		
		/**
		 * process merge tree
		 **/
		while ( itodo.size() )
		{
			MergeStrategyBlock * p = itodo.front();
			itodo.pop_front();
			
			std::cerr << "[V] processing ";
			p->printBase(std::cerr);
			std::cerr << std::endl;

			std::ostringstream tmpstr;
			tmpstr << tmpfilenamebase << "_" << std::setfill('0') << std::setw(6) << (mtmpid++);
			
			if ( dynamic_cast<MergeStrategyMergeInternalBlock *>(p) )
			{
				mergeBlocks(
					*(dynamic_cast<MergeStrategyMergeInternalBlock *>(p)),
					fn,fs,
					tmpstr.str(),
					rlencoderblocksize,
					lfblockmult,
					numthreads,
					chist,bwtterm,
					huftreefilename
				);
			}
			else if ( dynamic_cast<MergeStrategyMergeInternalSmallBlock *>(p) )
			{
				mergeBlocks(
					*(dynamic_cast<MergeStrategyMergeInternalSmallBlock *>(p)),
					fn,fs,
					tmpstr.str(),
					rlencoderblocksize,
					lfblockmult,
					numthreads,
					chist,bwtterm,
					huftreefilename
				);
			}
			else if ( dynamic_cast<MergeStrategyMergeExternalBlock *>(p) )
			{
				mergeBlocks(
					*(dynamic_cast<MergeStrategyMergeExternalBlock *>(p)),
					fn,fs,
					tmpstr.str(),
					rlencoderblocksize,
					lfblockmult,
					numthreads,
					chist,bwtterm,
					mem,
					huftreefilename
				);
			}

			bool const pfinished = p->parent && p->parent->childFinished();
					
			if ( pfinished )
				itodo.push_back(p->parent);
				
			#if defined(FERAMANZGEN_MEMORY_DEBUG)
			std::cerr << "[M"<< (mcnt++) << "] " << libmaus::util::MemUsage() << " " << libmaus::autoarray::AutoArrayMemUsage() << std::endl;
			#endif
		}

		#if 0
		uint64_t const memperthread = maxblocksize*sizeof(uint32_t) + maxhwtsize;
		#endif

		uint64_t const memperthread = (mem + numthreads-1)/numthreads;

		::libmaus::suffixsort::BwtMergeBlockSortResult const mergeresult = mergetree->sortresult;
		
		#if defined(HUFRL)
		::libmaus::huffman::RLEncoderStd::concatenate(mergeresult.getFiles().getBWT(),outfn,true /* removeinput */);
		#else
		::libmaus::gamma::GammaRLEncoder::concatenate(mergeresult.getFiles().getBWT(),outfn,true /* removeinput */);
		#endif
		//rename ( mergeresult.getFiles().getBWT().c_str(), outfn.c_str() );
		
		// remove hwt request for term symbol hwt
		remove ( mergeresult.getFiles().getHWTReq().c_str() );
		
		// remove hwt (null op)
		std::string const debhwt = ::libmaus::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hwt" + ".deb";
		rename ( mergeresult.getFiles().getHWT().c_str(), debhwt.c_str() );
		remove ( debhwt.c_str() );
		
		// remove gt files
		for ( uint64_t i = 0; i < mergeresult.getFiles().getGT().size(); ++i )
			remove ( mergeresult.getFiles().getGT()[i].c_str() );
		
		std::string const mergedisaname = mergeresult.getFiles().getSampledISA();
			
		std::cerr << "[V] computing Huffman shaped wavelet tree of final BWT...";	
		std::string const outhwt = ::libmaus::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hwt";
		libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type pICHWT;
		if ( input_types_type::utf8Wavelet() )
		{
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(
				RlToHwtBase<true>::rlToHwt(outfn, outhwt, tmpfilenamebase+"_finalhwttmp")
			);
			pICHWT = UNIQUE_PTR_MOVE(tICHWT);
		}
		else
		{
			libmaus::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(
				RlToHwtBase<false>::rlToHwt(outfn, outhwt, tmpfilenamebase+"_finalhwttmp")
			);		
			pICHWT = UNIQUE_PTR_MOVE(tICHWT);
		}
		std::cerr << "done, " << std::endl;
		
		std::cerr << "[V] loading Huffman shaped wavelet tree of final BWT...";	
		::libmaus::lf::ImpCompactHuffmanWaveletLF IHWT(pICHWT);
		std::cerr << "done." << std::endl;

		// sort the sampled isa file	
		uint64_t const blockmem = memperthread; // memory per thread
		sortIsaFile(mergedisaname,blockmem);

		// compute sampled suffix array and sampled inverse suffix array
		computeSampledSA(
			fn,fs,IHWT,mergedisaname,outfn,tmpfilenamebase,
			numthreads,lfblockmult,sasamplingrate,isasamplingrate,blockmem
		);

		// serialise character histogram
		std::string const outhist = ::libmaus::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hist";
		::libmaus::aio::CheckedOutputStream::unique_ptr_type Phistout(new ::libmaus::aio::CheckedOutputStream(outhist));
		::libmaus::util::NumberMapSerialisation::serialiseMap(*Phistout,chistnoterm);
		Phistout->flush();
		Phistout->close();
		Phistout.reset();
		
		// std::cerr << "[V] mergeresult.blockp0rank=" << mergeresult.blockp0rank << std::endl;
		
		return EXIT_SUCCESS;
	}
};

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		libmaus::timing::RealTimeClock rtc; rtc.start();

		#if defined(_OPENMP)
		unsigned int const maxthreads = omp_get_max_threads();
		unsigned int const numthreads = arginfo.getValue<unsigned int>("numthreads", BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultNumThreads());
		omp_set_num_threads(numthreads);
		#endif
		
		if ( arginfo.helpRequested() || ! arginfo.restargs.size() )
		{
			::libmaus::exception::LibMausException se;
			
			std::ostream & str = se.getStream();
			
			str << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			str << std::endl;
			str << "usage: " << arginfo.progname << " [options] <inputfile>" << std::endl;
			str << std::endl;
			str << "options:" << std::endl;
			str << "inputtype=[<bytestream>] (bytestream,compactstream,pac,pacterm,lz4,utf-8)" << std::endl;
			str << "outputfilename=[<"<< arginfo.getDefaultTmpFileName()+".bwt" << ">] (name of output .bwt file)" << std::endl;
			str << "sasamplingrate=[" << BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultSaSamplingRate() << "] sampling rate for sampled suffix array"<< std::endl;
			str << "isasamplingrate=[" << BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultIsaSamplingRate() << "] sampling rate for sampled inverse suffix array"<< std::endl;
			// str << "blocksize=[" << BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultBlockSize(BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultMem(),numthreads) << "] block size" << std::endl;
			str << "mem=[" << BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultMem() << "] memory target" << std::endl;
			#if defined(_OPENMP)
			str << "numthreads=[" << BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::getDefaultNumThreads() << "] number of threads" << std::endl;
			#endif
			// blocksize
			
			se.finish();
			throw se;
		}
		
		
		std::string const inputtype = arginfo.getValue<std::string>("inputtype","bytestream");
		
		if ( inputtype == "compactstream" )
			BwtMergeSort<libmaus::suffixsort::CompactInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "pac" )
			BwtMergeSort<libmaus::suffixsort::PacInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "pacterm" )
			BwtMergeSort<libmaus::suffixsort::PacTermInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "lz4" )
			BwtMergeSort<libmaus::suffixsort::Lz4InputTypes>::computeBwt(arginfo);
		else if ( inputtype == "utf-8" )
		{
			// compute index of file for random access, if it does not already exist
			std::string const fn = arginfo.getRestArg<std::string>(0);
			std::string const idxfn = fn + ".idx";
			if ( ! ::libmaus::util::GetFileSize::fileExists(idxfn) )
			{
				::libmaus::util::Utf8BlockIndex::unique_ptr_type index =
					UNIQUE_PTR_MOVE(::libmaus::util::Utf8BlockIndex::constructFromUtf8File(fn));
				::libmaus::aio::CheckedOutputStream COS(idxfn);
				index->serialise(COS);
				COS.flush();
				COS.close();
			}
			BwtMergeSort<libmaus::suffixsort::Utf8InputTypes>::computeBwt(arginfo);
		}
		else if ( inputtype == "bytestream" )
			BwtMergeSort<libmaus::suffixsort::ByteInputTypes>::computeBwt(arginfo);			
		else
		{
			libmaus::exception::LibMausException se;
			se.getStream() << "Unknown input type " << inputtype << std::endl;
			se.finish();
			throw se;
		}
			
		#if defined(_OPENMP)
		omp_set_num_threads(maxthreads);		
		#endif
		
		std::cerr << "[M] " << libmaus::util::MemUsage() << " runtime " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
