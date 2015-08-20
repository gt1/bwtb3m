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

#include <libmaus2/suffixsort/BwtMergeBlockSortRequest.hpp>
#include <libmaus2/wavelet/RlToHwtTermRequest.hpp>
#include <libmaus2/suffixsort/BwtMergeBlockSortRequestBase.hpp>
#include <libmaus2/wavelet/RlToHwtBase.hpp>

#include <libmaus2/wavelet/ImpCompactHuffmanWaveletTree.hpp>
#include <libmaus2/wavelet/ImpExternalWaveletGeneratorCompactHuffman.hpp>
#include <libmaus2/wavelet/ImpExternalWaveletGeneratorCompactHuffmanParallel.hpp>

#include <iostream>
#include <iomanip>

#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/aio/CircularWrapper.hpp>
#include <libmaus2/aio/FileFragment.hpp>
#include <libmaus2/aio/ReorderConcatGenericInput.hpp>

#include <libmaus2/autoarray/AutoArray.hpp>

#include <libmaus2/bitio/BitStreamFileDecoder.hpp>
#include <libmaus2/bitio/CompactDecoderBuffer.hpp>
#include <libmaus2/bitio/BitVectorInput.hpp>
#include <libmaus2/bitio/BitVectorOutput.hpp>

#include <libmaus2/gamma/GammaGapEncoder.hpp>
#include <libmaus2/gamma/GammaGapDecoder.hpp>
#include <libmaus2/gamma/SparseGammaGapFileSet.hpp>
#include <libmaus2/gamma/SparseGammaGapFileLevelSet.hpp>
#include <libmaus2/gamma/SparseGammaGapMultiFileLevelSet.hpp>
#include <libmaus2/gamma/SparseGammaGapEncoder.hpp>
#include <libmaus2/gamma/SparseGammaGapDecoder.hpp>

#include <libmaus2/huffman/GapDecoder.hpp>
#include <libmaus2/huffman/GapEncoder.hpp>
#include <libmaus2/huffman/HuffmanTree.hpp>

#include <libmaus2/gamma/GammaGapEncoder.hpp>
#include <libmaus2/gamma/GammaGapDecoder.hpp>
#include <libmaus2/gamma/GammaRLEncoder.hpp>
#include <libmaus2/gamma/GammaRLDecoder.hpp>

#include <libmaus2/huffman/huffman.hpp>
#include <libmaus2/huffman/HuffmanEncoderFile.hpp>
#include <libmaus2/huffman/RLEncoder.hpp>
#include <libmaus2/huffman/RLDecoder.hpp>

#include <libmaus2/lf/DArray.hpp>
#include <libmaus2/lf/LF.hpp>
#include <libmaus2/lf/ImpCompactHuffmanWaveletLF.hpp>

#include <libmaus2/math/numbits.hpp>

#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/parallel/LockedBool.hpp>

#include <libmaus2/sorting/PairFileSorting.hpp>

#include <libmaus2/suffixsort/BwtMergeBlockSortResult.hpp>
#include <libmaus2/suffixsort/BwtMergeTempFileNameSet.hpp>
#include <libmaus2/suffixsort/BwtMergeTempFileNameSetVector.hpp>
#include <libmaus2/suffixsort/BwtMergeZBlock.hpp>
#include <libmaus2/suffixsort/BwtMergeZBlockRequest.hpp>
#include <libmaus2/suffixsort/BwtMergeZBlockRequestVector.hpp>
#include <libmaus2/suffixsort/ByteInputTypes.hpp>
#include <libmaus2/suffixsort/CircularBwt.hpp>
#include <libmaus2/suffixsort/CircularSuffixComparator.hpp>
#include <libmaus2/suffixsort/CompactInputTypes.hpp>
#include <libmaus2/suffixsort/divsufsort.hpp>
#include <libmaus2/suffixsort/GapArrayByte.hpp>
#include <libmaus2/suffixsort/GapMergePacket.hpp>
#include <libmaus2/suffixsort/PacInputTypes.hpp>
#include <libmaus2/suffixsort/PacTermInputTypes.hpp>
#include <libmaus2/suffixsort/Utf8InputTypes.hpp>
#include <libmaus2/suffixsort/Lz4InputTypes.hpp>

#include <libmaus2/wavelet/ImpExternalWaveletGeneratorHuffman.hpp>
#include <libmaus2/wavelet/ImpExternalWaveletGeneratorHuffmanParallel.hpp>
#include <libmaus2/wavelet/ImpHuffmanWaveletTree.hpp>
#include <libmaus2/wavelet/Utf8ToImpCompactHuffmanWaveletTree.hpp>
#include <libmaus2/wavelet/Utf8ToImpHuffmanWaveletTree.hpp>

#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/BorderArray.hpp>
#include <libmaus2/util/FileTempFileContainer.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/util/HistogramSet.hpp>
#include <libmaus2/util/KMP.hpp>
#include <libmaus2/util/MemUsage.hpp>
#include <libmaus2/util/NumberMapSerialisation.hpp>
#include <libmaus2/util/OctetString.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/util/StringSerialisation.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/util/Utf8String.hpp>
#include <libmaus2/util/SimpleCountingHash.hpp>
#include <libmaus2/util/SuccinctBorderArray.hpp>

libmaus2::parallel::OMPLock gcerrlock;
libmaus2::parallel::OMPLock glock;

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
	typedef libmaus2::util::shared_ptr<MergeStrategyBlock>::type shared_ptr_type;

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
	::libmaus2::suffixsort::BwtMergeBlockSortResult sortresult;
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
	virtual uint64_t directSortSpace() const
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
	virtual void fillQueryObjects(libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV) = 0;
	// fill zblock ranks from finished base blocks for t threads
	virtual void fillGapRequestObjects(uint64_t const t) = 0;
	
	virtual bool isLeaf() const = 0;
	virtual void setParent(MergeStrategyBlock * rparent) = 0;
	virtual bool childFinished() = 0;
};

struct MergeStrategyBaseBlock : public MergeStrategyBlock
{
	libmaus2::suffixsort::BwtMergeBlockSortRequest sortreq;
	std::vector<uint64_t> querypos;

	MergeStrategyBaseBlock() : MergeStrategyBlock() {}

	MergeStrategyBaseBlock(
		uint64_t const rlow, uint64_t const rhigh, 
		::libmaus2::huffman::HuffmanTree::EncodeTable & EC,
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
		::libmaus2::huffman::HuffmanTree::EncodeTable & EC,
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
		::libmaus2::suffixsort::BwtMergeZBlockRequestVector zreqvec;
		zreqvec.resize(querypos.size());
		for ( uint64_t i = 0; i < querypos.size(); ++i )
			zreqvec[i] = ::libmaus2::suffixsort::BwtMergeZBlockRequest(querypos[i]);
		sortreq.zreqvec = zreqvec;
	}

	void fillQueryObjects(libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV)
	{
		libmaus2::autoarray::AutoArray < ::libmaus2::suffixsort::BwtMergeZBlock > const & Z = sortresult.getZBlocks();
		
		for ( uint64_t i = 0; i < Z.size(); ++i )
		{
			uint64_t const p = Z[i].getZAbsPos();
			uint64_t const r = Z[i].getZRank();
			
			// search objects with matching position
			typedef libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject>::iterator it;
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
		libmaus2::exception::LibMausException se;
		se.getStream() << "childFinished called on an object of struct MergeStrategyBaseBlock" << std::endl;
		se.finish();
		throw se;
	}

	/**
	 * @return space required in bytes for sorting block in internal memory using divsufsort
	 **/
	virtual uint64_t directSortSpace() const
	{
		if ( ((high-low)+sortreq.lcpnext) < (1ull << 31) )
		{
			return 
				((high-low)+sortreq.lcpnext)*(sizeof(uint32_t)) + // suffix array
				(high-low+7)/8 + // GT bit vector
				(sourcetextindexbits+7)/8 +
				sourcelengthbytes;
		}
		else
		{
			return ((high-low)+sortreq.lcpnext)*(sizeof(uint64_t)) + // suffix array
				(high-low+7)/8 + // GT bit vector
				(sourcetextindexbits+7)/8 +
				sourcelengthbytes;
		}
	}
};

// #define HUFGAP
#define HUFRL

#if defined(HUFRL)
typedef ::libmaus2::huffman::RLEncoderStd rl_encoder;
typedef ::libmaus2::huffman::RLDecoder rl_decoder;
#else
typedef ::libmaus2::gamma::GammaRLEncoder rl_encoder;
typedef ::libmaus2::gamma::GammaRLDecoder rl_decoder;
#endif

/**
 * sorting thread for base blocks
 **/
struct BaseBlockSortThread : public libmaus2::parallel::PosixThread
{
	typedef BaseBlockSortThread this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	//! thread id
	uint64_t tid;

	/**
	 * semaphore. delivers a message whenever there is sufficient
	 * free space to process the next element
	 **/
	libmaus2::parallel::PosixSemaphore & P;
	//! vector of blocks to be processed
	std::vector < MergeStrategyBlock::shared_ptr_type > & V;

	//! next package to be processed
	volatile uint64_t & next;
	//! amount of free memory
	volatile uint64_t & freemem;
	//! number of finished threads
	volatile uint64_t & finished;
	//! lock for the above
	libmaus2::parallel::PosixMutex & freememlock;
	//! inner node queue
	std::deque<MergeStrategyBlock *> & itodo;
	//! pending
	std::deque<uint64_t> & pending;
		
	BaseBlockSortThread(
		uint64_t rtid,
		libmaus2::parallel::PosixSemaphore & rP,
		std::vector < MergeStrategyBlock::shared_ptr_type > & rV,
		uint64_t & rnext,
		uint64_t & rfreemem,
		uint64_t & rfinished,
		libmaus2::parallel::PosixMutex & rfreememlock,
		std::deque<MergeStrategyBlock *> & ritodo,
		std::deque<uint64_t> & rpending
	) : tid(rtid), P(rP), V(rV), next(rnext), freemem(rfreemem), finished(rfinished), freememlock(rfreememlock),
	    itodo(ritodo), pending(rpending)
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
				libmaus2::parallel::ScopePosixMutex scopelock(freememlock);

				if ( pending.size() )
				{
					pack = pending.front();	
					pending.pop_front();					
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
					block->sortresult = ::libmaus2::suffixsort::BwtMergeBlockSortResult::load(block->sortreq.dispatch<rl_encoder>());
				}
				catch(std::exception const & ex)
				{
					libmaus2::parallel::ScopePosixMutex scopelock(freememlock);
					std::cerr << tid << " failed " << pack << " " << ex.what() << std::endl;
				}

				{
					// get lock
					libmaus2::parallel::ScopePosixMutex scopelock(freememlock);
					
					std::cerr << "[V] [" << tid << "] sorted block " << pack << std::endl;

					if ( V[pack]->parent )
					{
						bool const pfinished = V[pack]->parent->childFinished();
					
						if ( pfinished )
						{
							itodo.push_back(V[pack]->parent);
						}
					}
					
					// "free" memory
					freemem += V[pack]->directSortSpace();
					
					// post if there is room for another active sorting thread
					while ( next < V.size() && freemem >= V[next]->directSortSpace() )
					{
						freemem -= V[next]->directSortSpace();
						pending.push_back(static_cast<uint64_t>(next));
						next += 1;
						P.post();
					}
					
					if ( next == V.size() )
						P.post();
				}
			}
		}
		
		{
		libmaus2::parallel::ScopePosixMutex scopelock(freememlock);
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
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	std::vector < MergeStrategyBlock::shared_ptr_type > & V;
	
	libmaus2::parallel::PosixSemaphore P;
	uint64_t next;
	uint64_t freemem;
	uint64_t finished;
	libmaus2::parallel::PosixMutex freememlock;
	//! inner node queue
	std::deque<MergeStrategyBlock *> & itodo;
	std::deque<uint64_t> pending;

	libmaus2::autoarray::AutoArray<BaseBlockSortThread::unique_ptr_type> threads;

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
				libmaus2::exception::LibMausException se;
				se.getStream() << "Memory provided is " << freemem << " but " 
					<< V[i]->directSortSpace() << " are required for sorting block " << i << std::endl;
				se.finish();
				throw se;
			}
	
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			BaseBlockSortThread::unique_ptr_type tthreadsi(
				new BaseBlockSortThread(i,P,V,next,freemem,finished,freememlock,itodo,pending)
			);
			threads[i] = UNIQUE_PTR_MOVE(tthreadsi);

		}
	}
	
	void setup()
	{
		while ( next < V.size() && freemem >= V[next]->directSortSpace() )
		{
			freemem -= V[next]->directSortSpace();
			pending.push_back(next);
			next += 1;
		}
		uint64_t const p = pending.size();
		for ( uint64_t i = 0; i < p; ++i )
			P.post();
	}
	
	void start(uint64_t const stacksize)
	{
		for ( uint64_t i = 0; i < threads.size(); ++i )
			threads[i]->startStack(stacksize);
			
		setup();
	}

	void start()
	{
		for ( uint64_t i = 0; i < threads.size(); ++i )
			threads[i]->start();
			
		setup();
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
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::vector<MergeStrategyBlock::shared_ptr_type> const * pchildren;
	uint64_t into;
	std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > zblocks;
	
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

	libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> getQueryPositionObjects(uint64_t const t)
	{
		std::vector<uint64_t> const Q = getQueryPositions(t);
		libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> VV(Q.size());
		
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
	void fillQueryObjects(libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> & VV)
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
			libmaus2::autoarray::AutoArray<MergeStrategyMergeGapRequestQueryObject> VV = gaprequests[i]->getQueryPositionObjects(/* VV, */t);
			std::sort(VV.begin(),VV.end());
			// fill the objects (add up ranks for each query position)
			children[gaprequests[i]->into]->fillQueryObjects(VV);
			// size of vector
			uint64_t const vn = VV.size();

			// push z blocks back to front
			for ( uint64_t i = 0; i < VV.size(); ++i )
			{
				uint64_t const ii = vn-i-1;
				VV[i].o->zblocks.push_back(::libmaus2::suffixsort::BwtMergeZBlock(VV[ii].p,VV[ii].r));
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
		libmaus2::exception::LibMausException se;
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
		typedef typename ::libmaus2::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
		uint64_t const fs = stream_type::getFileSize(fn);
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		uint64_t const symsperfrag = (fs + numthreads - 1)/numthreads;
		uint64_t const numfrags = (fs + symsperfrag - 1)/symsperfrag;

		libmaus2::util::HistogramSet HS(numfrags,256);
		
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

			::libmaus2::autoarray::AutoArray<unsigned_char_type> B(16*1024,false);
			libmaus2::util::Histogram & H = HS[t];
			
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

		::libmaus2::util::Histogram::unique_ptr_type PH(HS.merge());
		::libmaus2::util::Histogram & H(*PH);

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
		typedef ::libmaus2::huffman::BitInputBuffer4 sbis_type;			
		::libmaus2::aio::InputStreamInstance::unique_ptr_type istr;
		sbis_type::raw_input_ptr_type ript;
		sbis_type::unique_ptr_type SBIS;
		
		::libmaus2::aio::InputStreamInstance::unique_ptr_type tistr(new ::libmaus2::aio::InputStreamInstance(fn));
		istr = UNIQUE_PTR_MOVE(tistr);
		sbis_type::raw_input_ptr_type tript(new sbis_type::raw_input_type(*istr));
		ript = UNIQUE_PTR_MOVE(tript);
		sbis_type::unique_ptr_type tSBIS(new sbis_type(ript,static_cast<uint64_t>(64*1024)));
		SBIS = UNIQUE_PTR_MOVE(tSBIS);
		
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
		::libmaus2::lf::ImpCompactHuffmanWaveletLF const & IHWT,
		::libmaus2::aio::SynchronousGenericOutput<uint64_t> & SGO,
		::libmaus2::aio::SynchronousGenericOutput<uint64_t> & ISGO,
		uint64_t const sasamplingrate,
		uint64_t const isasamplingrate,
		int64_t const ibs = -1
	)
	{
		assert ( ::libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(sasamplingrate) == 1 );
		assert ( ::libmaus2::rank::PopCnt8<sizeof(unsigned long)>::popcnt8(isasamplingrate) == 1 );
		
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

			rl_encoder rlenc(outputfilename,0 /* alphabet */,0,rlencoderblocksize);
		
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
				
				// ::libmaus2::util::GetFileSize::copy(bwtfilenames[0][i],outputfilename);
				libmaus2::aio::OutputStreamFactoryContainer::rename ( bwtfilenames[0][i].c_str(), outputfilename.c_str() );
				
				outputfilenames.push_back(outputfilename);
			}
			
			return outputfilenames;
		}
		// at least one gap file, merge
		else
		{
			#if defined(HUFGAP)
			typedef ::libmaus2::huffman::GapDecoder gapfile_decoder_type;
			#else
			typedef ::libmaus2::gamma::GammaGapDecoder gapfile_decoder_type;
			#endif

			unsigned int const albits = rl_decoder::haveAlphabetBits() ? rl_decoder::getAlBits(bwtfilenames.front()) : 0;
			
			uint64_t const firstblockgapfilesize = gapfilenames.size() ? gapfile_decoder_type::getLength(gapfilenames[0]) : 0;
			assert ( firstblockgapfilesize );

			// uint64_t const fs = rl_decoder::getLength(bwtfilenames);
			uint64_t const fs = getRlLength(bwtfilenames);

			// first gap file meta information
			::libmaus2::huffman::IndexDecoderDataArray::unique_ptr_type Pfgapidda(new ::libmaus2::huffman::IndexDecoderDataArray(gapfilenames[0]));
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
			::libmaus2::autoarray::AutoArray< ::libmaus2::suffixsort::GapMergePacket> gmergepackets(maxgparts);
			// actual merge parts
			uint64_t actgparts = 0;
			uint64_t gs = 0;
			
			// std::cerr << "fs=" << fs << " tgparts=" << tgparts << " gpartsize=" << gpartsize << std::endl;
						
			// compute number of suffixes per gblock
			while ( hlow != firstblockgapfilesize )
			{
				uint64_t const gsrest = (fs-gs);
				uint64_t const gskip = std::min(gsrest,gpartsize);
				libmaus2::huffman::KvInitResult kvinit;
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
					::libmaus2::exception::LibMausException se;
					se.getStream() << "acgtgparts=" << actgparts << " >= maxgparts=" << maxgparts << std::endl;
					se.finish();
					throw se;
					// assert ( actgparts < maxgparts );
				}
				#endif
				
				// save interval on first gap array and number of suffixes on this interval
				gmergepackets[actgparts++] = ::libmaus2::suffixsort::GapMergePacket(hlow,hhigh,s);
				
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
			::libmaus2::autoarray::AutoArray<uint64_t> spref(actgparts+1,false);
			for ( uint64_t i = 0; i < actgparts; ++i )
				spref[i] = gmergepackets[i].s;
			spref.prefixSums();
			
			// compute prefix sums over number of suffixes per block used for each gpart
			::libmaus2::autoarray::AutoArray < uint64_t > bwtusedcntsacc( (actgparts+1)*(gapfilenames.size()+1), false );
			
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
				::libmaus2::autoarray::AutoArray < gapfile_decoder_type::unique_ptr_type > gapdecs(gapfilenames.size());

				for ( uint64_t j = 0; j < gapfilenames.size(); ++j )
				{
					::libmaus2::huffman::KvInitResult kvinitresult;
					gapfile_decoder_type::unique_ptr_type tgapdecsj(
						new gapfile_decoder_type(
							gapfilenames[j],
							lspref,kvinitresult
						)						
					);
					gapdecs[j] = UNIQUE_PTR_MOVE(tgapdecsj);
					
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
			::libmaus2::autoarray::AutoArray < uint64_t > bwtusedcnts( (gapfilenames.size()+1) * actgparts, false );
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
				::libmaus2::util::TempFileRemovalContainer::addTempFile(gpartfrag);
				gpartfrags[z] = gpartfrag;
			
				#if 0
				std::cerr << gmergepackets[z] << ",spref=" << spref[z] << std::endl;
				#endif
				uint64_t lspref = spref[z];
				::libmaus2::autoarray::AutoArray < gapfile_decoder_type::unique_ptr_type > gapdecoders(gapfilenames.size());
				::libmaus2::autoarray::AutoArray< uint64_t > gapcur(gapfilenames.size());

				// set up gap file decoders at the proper offsets
				for ( uint64_t j = 0; j < gapfilenames.size(); ++j )
				{
					// sum up number of suffixes in later blocks for this gpart
					uint64_t suflat = 0;
					for ( uint64_t k = j+1; k < bwtfilenames.size(); ++k )
						suflat += bwtusedcnts [ k*actgparts + z ];
				
					::libmaus2::huffman::KvInitResult kvinitresult;
					gapfile_decoder_type::unique_ptr_type tgapdecodersj(
						new gapfile_decoder_type(
							gapfilenames[j],
							lspref,kvinitresult
						)					
					);
					gapdecoders[j] = UNIQUE_PTR_MOVE(tgapdecodersj);
					if ( suflat )
						gapcur[j] = gapdecoders[j]->decode();
					else
						gapcur[j] = 0;
						
					lspref = kvinitresult.voffset + kvinitresult.kvtarget;

					if ( j == 0 )
						assert ( kvinitresult.kvtarget == 0 );
				}				

				::libmaus2::autoarray::AutoArray < uint64_t > bwttowrite(bwtfilenames.size(),false);
				::libmaus2::autoarray::AutoArray < rl_decoder::unique_ptr_type > bwtdecoders(bwtfilenames.size());
				
				for ( uint64_t j = 0; j < bwtfilenames.size(); ++j )
				{
					uint64_t const bwtoffset = bwtusedcntsacc [ j * (actgparts+1) + z ];
					bwttowrite[j] = bwtusedcnts [ j * actgparts + z ];
					
					rl_decoder::unique_ptr_type tbwtdecodersj(
						new rl_decoder(bwtfilenames[j],bwtoffset)					
					);
					bwtdecoders[j] = UNIQUE_PTR_MOVE(tbwtdecodersj);
					
					#if 0
					std::cerr << "block=" << j << " offset=" << bwtoffset << " bwttowrite=" << bwttowrite[j] << std::endl;
					#endif
				}
				
				uint64_t const totalbwt = std::accumulate(bwttowrite.begin(),bwttowrite.end(),0ull);

				rl_encoder bwtenc(gpartfrag,albits /* alphabet */,totalbwt,rlencoderblocksize);
				
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
		::libmaus2::timing::RealTimeClock rtc;
		
		// merge sampled inverse suffix arrays		
		std::cerr << "[V] merging sampled inverse suffix arrays...";
		rtc.start();
		::libmaus2::aio::SynchronousGenericInput<uint64_t> SGIISAold(oldmergedisaname,16*1024);
		::libmaus2::aio::SynchronousGenericInput<uint64_t> SGIISAnew(newmergedisaname,16*1024);
		::libmaus2::aio::SynchronousGenericOutput<uint64_t> SGOISA(mergedmergedisaname,16*1024);
		
		// sum over old (RHS block) suffixes
		uint64_t s = 0;
		// rank of position zero in merged block
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
			}

			// is next sample in next (LHS) block for rank i?
			if ( SGIISAnew.peek() == static_cast<int64_t>(i) )
			{
				// add number of old suffixes to rank
				uint64_t const r = SGIISAnew.get() + s;
				// keep absolute position as is
				uint64_t const p = SGIISAnew.get();
				
				SGOISA . put ( r );
				SGOISA . put ( p );
			
				// check whether this rank is for the leftmost position in the merged block	
				if ( p == blockstart )
					blockp0rank = r;
			}
		}
		
		assert ( SGIISAnew.peek() < 0 );
		assert ( SGIISAold.peek() < 0 );

		SGOISA.flush();
		std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
	
		return blockp0rank;
	}

	static void saveGapFile(
		::libmaus2::autoarray::AutoArray<uint32_t> const & G,
		std::string const & gapfile
	)
	{
		::libmaus2::timing::RealTimeClock rtc;
	
		std::cerr << "[V] saving gap file...";
		rtc.start();
		#if defined(HUFGAP)
		::libmaus2::util::Histogram gaphist;
		for ( uint64_t i = 0; i < G.size(); ++i )
			gaphist ( G[i] );
		::libmaus2::huffman::GapEncoder GE(gapfile,gaphist,G.size());
		GE.encode(G.begin(),G.end());
		GE.flush();
		#else
		::libmaus2::gamma::GammaGapEncoder GE(gapfile);
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
		::libmaus2::timing::RealTimeClock rtc;
		rtc.start();
		
		std::vector< std::string > allfiles(gtpartnames.begin(),gtpartnames.end());
		allfiles.push_back(newgtpart);
		
		#if 0
		// encoder for new gt stream
		libmaus2::bitio::BitVectorOutput GTHEFref(newmergedgtname);
		libmaus2::bitio::BitVectorInput BVI(allfiles);
		uint64_t const tn = libmaus2::bitio::BitVectorInput::getLength(allfiles);
		for ( uint64_t i = 0; i < tn; ++i )
			GTHEFref.writeBit(BVI.readBit());
		GTHEFref.flush();
		#endif

		libmaus2::bitio::BitVectorOutput GTHEF(newmergedgtname);
		
		unsigned int prevbits = 0;
		uint64_t prev = 0;

		// append part streams
		for ( uint64_t z = 0; z < allfiles.size(); ++z )
		{
			uint64_t const n = libmaus2::bitio::BitVectorInput::getLength(allfiles[z]);
			uint64_t const fullwords = n / 64;
			uint64_t const restbits = n-fullwords*64;
			libmaus2::aio::SynchronousGenericInput<uint64_t> SGI(allfiles[z],8192);
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
			
			libmaus2::aio::FileRemoval::removeFile(gtpartnames[z].c_str());
		}
		
		for ( uint64_t i = 0; i < prevbits; ++i )
			GTHEF.writeBit((prev >> (63-i)) & 1);
		
		// flush gt stream
		GTHEF.flush();
		
		std::cerr << "[V] concatenated bit vectors in time " << rtc.getElapsedSeconds() << std::endl;
	}
	
	struct GapArrayComputationResult
	{
		::libmaus2::autoarray::AutoArray<uint32_t> G;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		GapArrayComputationResult()
		: zactive(0)
		{
		
		}
		
		GapArrayComputationResult(
			::libmaus2::autoarray::AutoArray<uint32_t> & rG, 
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus2::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: G(rG), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};

	struct GapArrayByteComputationResult
	{
		libmaus2::suffixsort::GapArrayByte::shared_ptr_type G;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		GapArrayByteComputationResult()
		: zactive(0)
		{
		
		}
		
		GapArrayByteComputationResult(
			libmaus2::suffixsort::GapArrayByte::shared_ptr_type rG, 
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus2::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: G(rG), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};

	struct SparseGapArrayComputationResult
	{
		// name of gap file
		std::vector<std::string> fn;
		std::vector < std::string > gtpartnames;
		uint64_t zactive;
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos;
		
		SparseGapArrayComputationResult()
		: zactive(0)
		{
		
		}
		
		SparseGapArrayComputationResult(
			std::vector < std::string > const & rfn,
			std::vector < std::string > const & rgtpartnames, 
			uint64_t const rzactive,
			::libmaus2::autoarray::AutoArray<uint64_t> & rzabsblockpos
		)
		: fn(rfn), gtpartnames(rgtpartnames), zactive(rzactive), zabsblockpos(rzabsblockpos) {}
	};
	
	static libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ensureWaveletTreeGenerated(::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults)
	{
		// generate wavelet tree from request if necessary
		if ( ! libmaus2::util::GetFileSize::fileExists(blockresults.getFiles().getHWT()) )
		{
			libmaus2::timing::RealTimeClock rtc; rtc.start();
			std::cerr << "[V] Generating HWT for gap file computation...";
			assert ( libmaus2::util::GetFileSize::fileExists(blockresults.getFiles().getHWTReq() ) );	
			libmaus2::wavelet::RlToHwtTermRequest::unique_ptr_type ureq(libmaus2::wavelet::RlToHwtTermRequest::load(blockresults.getFiles().getHWTReq()));
			libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type	ptr(ureq->dispatch<rl_decoder>());
			libmaus2::aio::FileRemoval::removeFile ( blockresults.getFiles().getHWTReq().c_str() );
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			return UNIQUE_PTR_MOVE(ptr);			
		}
		else
		{
			libmaus2::aio::InputStreamInstance CIS(blockresults.getFiles().getHWT());
			libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ptr(new libmaus2::wavelet::ImpCompactHuffmanWaveletTree(CIS));
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
		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults, // information on block
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		::libmaus2::lf::DArray * const accD, // accumulated symbol freqs for block
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks // lf starting points
	)
	{
		// gap array
		::libmaus2::autoarray::AutoArray<uint32_t> G(cblocksize+1);
		
		// set up lf mapping
		::libmaus2::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
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
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus2::timing::RealTimeClock rtc; 
		rtc.start();
		#if defined(_OPENMP) && defined(LIBMAUS2_HAVE_SYNC_OPS)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t z = 0; z < static_cast<int64_t>(zactive); ++z )
		{
			::libmaus2::timing::RealTimeClock subsubrtc; subsubrtc.start();

			::libmaus2::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
			
			std::string const gtpartname = newmergedgtname + "_" + ::libmaus2::util::NumberSerialisation::formatNumber(z,4) + ".gt";
			::libmaus2::util::TempFileRemovalContainer::addTempFile(gtpartname);
			gtpartnames[z] = gtpartname;
			#if 0
			::libmaus2::huffman::HuffmanEncoderFileStd GTHEF(gtpartname);
			#endif
			libmaus2::bitio::BitVectorOutput GTHEF(gtpartname);

			#if 0
			::libmaus2::bitio::BitStreamFileDecoder gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			#endif
			libmaus2::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			
			typename input_types_type::circular_reverse_wrapper CRWR(fn,zblock.getZAbsPos() % fs);
			uint64_t r = zblock.getZRank();

			uint64_t const zlen = zabsblockpos [ z ] - zabsblockpos [z+1];
			
			for ( uint64_t i = 0; i < zlen; ++i )
			{
				GTHEF.writeBit(r > lp0);
				
				int64_t const sym = CRWR.get();			
				bool const gtf = gtfile.readBit();

				r = IHWL.step(sym,r) + ((sym == firstblocklast)?gtf:0);
			
				#if defined(LIBMAUS2_HAVE_SYNC_OPS)
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
		::libmaus2::lf::DArray * const accD // accumulated symbol freqs for block
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		// use gap object's zblocks vector
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
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
		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults, // information on block
		std::vector<std::string> const & mergedgtname, // previous gt file name
		std::string const & newmergedgtname, // new gt file name
		std::string const & gapoverflowtmpfilename, // gap overflow tmp file
		::libmaus2::lf::DArray * const accD, // accumulated symbol freqs for block
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks // lf starting points
	)
	{
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif

		// gap array
		::libmaus2::suffixsort::GapArrayByte::shared_ptr_type pG(
			new ::libmaus2::suffixsort::GapArrayByte(
				cblocksize+1,
				512, /* number of overflow words per thread */
				numthreads,
				gapoverflowtmpfilename
			)
		);
		::libmaus2::suffixsort::GapArrayByte & G = *pG;
		
		// set up lf mapping
		::libmaus2::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
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
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus2::timing::RealTimeClock rtc; 
		rtc.start();
		#if defined(_OPENMP) && defined(LIBMAUS2_HAVE_SYNC_OPS)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t z = 0; z < static_cast<int64_t>(zactive); ++z )
		{
			::libmaus2::timing::RealTimeClock subsubrtc; subsubrtc.start();

			::libmaus2::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
			
			std::string const gtpartname = newmergedgtname + "_" + ::libmaus2::util::NumberSerialisation::formatNumber(z,4) + ".gt";
			::libmaus2::util::TempFileRemovalContainer::addTempFile(gtpartname);
			gtpartnames[z] = gtpartname;
			#if 0
			::libmaus2::huffman::HuffmanEncoderFileStd GTHEF(gtpartname);
			#endif
			libmaus2::bitio::BitVectorOutput GTHEF(gtpartname);

			#if 0
			::libmaus2::bitio::BitStreamFileDecoder gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			#endif
			libmaus2::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
			
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
		::libmaus2::lf::DArray * const accD // accumulated symbol freqs for block
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		// use gap object's zblocks vector
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
		return computeGapArrayByte(fn,fs,blockstart,cblocksize,nextblockstart,mergeprocrightend,
			blockresults,mergedgtname,newmergedgtname,gapoverflowtmpfilename,accD,zblocks);
	}
	
	struct ZNext
	{
		uint64_t znext;
		uint64_t znextcount;
		libmaus2::parallel::OMPLock lock;
		
		ZNext(uint64_t const rznextcount) : znext(0), znextcount(rznextcount) {}
		
		bool getNext(uint64_t & next)
		{
			libmaus2::parallel::ScopeLock slock(lock);
			
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
		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults,
		std::vector<std::string> const & mergedgtname,
		std::string const & newmergedgtname,
		::libmaus2::lf::DArray * const accD,
		//
		std::string const & outputgapfilename,
		std::string const & tmpfileprefix,
		uint64_t const maxmem,
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks
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
		
		libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray<uint64_t> > GG(numthreads,false);
		for ( uint64_t i = 0; i < numthreads; ++i )
			GG[i] = libmaus2::autoarray::AutoArray<uint64_t>(wordsperthread,false);
		
		libmaus2::util::TempFileNameGenerator tmpgen(tmpfileprefix,3);
		libmaus2::gamma::SparseGammaGapMultiFileLevelSet SGGFS(tmpgen,numthreads);
	
		// set up lf mapping
		::libmaus2::lf::DArray D(static_cast<std::string const &>(blockresults.getFiles().getHist()));
		accD->merge(D);
		#if 0
		bool const hwtdelayed = 
			ensureWaveletTreeGenerated(blockresults);
		#endif
		::libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type ICHWL(ensureWaveletTreeGenerated(blockresults));
		::libmaus2::lf::ImpCompactHuffmanWaveletLF IHWL(ICHWL);
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
		::libmaus2::autoarray::AutoArray<uint64_t> zabsblockpos(zactive+1,false);
		for ( uint64_t z = 0; z < zactive; ++z )
			zabsblockpos[z] = zblocks[z].getZAbsPos();
		zabsblockpos [ zactive ] = blockstart + cblocksize;

		std::vector < std::string > gtpartnames(zactive);

		::libmaus2::timing::RealTimeClock rtc; 
		rtc.start();
		
		ZNext znext(zactive);
		
		uint64_t termcnt = 0;
		libmaus2::parallel::OMPLock termcntlock;
		
		libmaus2::parallel::PosixSemaphore qsem; // queue semaphore
		libmaus2::parallel::PosixSemaphore tsem; // term semaphore
		libmaus2::parallel::PosixSemaphore globsem; // meta semaphore for both above
		libmaus2::parallel::LockedBool termflag(false);
		libmaus2::parallel::LockedBool qterm(false);
		
		SGGFS.registerMergePackSemaphore(&qsem);
		SGGFS.registerMergePackSemaphore(&globsem);
		SGGFS.registerTermSemaphore(&tsem);
		SGGFS.registerTermSemaphore(&globsem);
		SGGFS.setTermSemCnt(numthreads);
		
		#if defined(_OPENMP)		
		#pragma omp parallel
		#endif
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
			
				::libmaus2::timing::RealTimeClock subsubrtc; subsubrtc.start();
				
				::libmaus2::suffixsort::BwtMergeZBlock const & zblock = zblocks[z];
				
				std::string const gtpartname = newmergedgtname + "_" + ::libmaus2::util::NumberSerialisation::formatNumber(z,4) + ".gt";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(gtpartname);
				gtpartnames[z] = gtpartname;
				libmaus2::bitio::BitVectorOutput GTHEF(gtpartname);

				libmaus2::bitio::BitVectorInput gtfile(mergedgtname, (mergeprocrightend - zblock.getZAbsPos()) );
				
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
					libmaus2::gamma::SparseGammaGapBlockEncoder::encodeArray(Ga,Gc,tfn);
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
					libmaus2::gamma::SparseGammaGapBlockEncoder::encodeArray(Ga,Gc,tfn);
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
		::libmaus2::lf::DArray * const accD, // accumulated symbol freqs for block
		std::string const & outputgapfilename,
		std::string const & tmpfileprefix,
		uint64_t const maxmem
	)
	{
		uint64_t const into = msmgr.into;

		std::vector<MergeStrategyBlock::shared_ptr_type> const & children =
			*(msmgr.pchildren);

		::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = children[into]->sortresult;

		uint64_t const blockstart = blockresults.getBlockStart();
		uint64_t const cblocksize = blockresults.getCBlockSize();
		uint64_t const nextblockstart = (blockstart+cblocksize)%fs;
		uint64_t const mergeprocrightend =
			children.at(children.size()-1)->sortresult.getBlockStart() +
			children.at(children.size()-1)->sortresult.getCBlockSize();
		std::vector < ::libmaus2::suffixsort::BwtMergeZBlock > const & zblocks = msmgr.zblocks;
		
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
		libmaus2::aio::FileRemoval::removeFile ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus2::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
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
			::libmaus2::lf::DArray::unique_ptr_type accD(new ::libmaus2::lf::DArray(
				sblockhist
				)
			);
			// first block
			::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = 
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
				::libmaus2::util::TempFileRemovalContainer::addTempFile(renamed);
				libmaus2::aio::OutputStreamFactoryContainer::rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
			}

			result.setGT(stringVectorAppend(GACR.gtpartnames,oldgtnames));

			// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
			libmaus2::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
			result.setBlockP0Rank( mergeIsa(
				mergereq.children[1]->sortresult.getFiles().getSampledISA(), // old sampled isa
				blockresults.getFiles().getSampledISA(), // new sampled isa
				result.getFiles().getSampledISA(),blockstart,mergeGO/*GACR.G.begin()*/,cblocksize+1 /* GACR.G.size() */
			) );
			
			::libmaus2::timing::RealTimeClock rtc; rtc.start();
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
			::libmaus2::timing::RealTimeClock wprtc; wprtc.start();
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
					+ ::libmaus2::util::NumberSerialisation::formatNumber(encfilenames.size(),6)
					+ ".bwt"
				);
				::libmaus2::util::TempFileRemovalContainer::addTempFile(encfilenames.back());
				ilow = ihigh;
			}
			assert ( wpacks.size() <= numthreads );
			// std::cerr << "done,time=" << wprtc.getElapsedSeconds() << ")";
			
			// std::cerr << "(setting up IDDs...";
			wprtc.start();
			
			unsigned int const albits = rl_decoder::haveAlphabetBits() ? rl_decoder::getAlBits(mergereq.children[0]->sortresult.getFiles().getBWT()) : 0;

			::libmaus2::huffman::IndexDecoderDataArray IDD0(
				mergereq.children[0]->sortresult.getFiles().getBWT());
			::libmaus2::huffman::IndexDecoderDataArray IDD1(
				mergereq.children[1]->sortresult.getFiles().getBWT());
			
			::libmaus2::huffman::IndexEntryContainerVector::unique_ptr_type IECV0 = ::libmaus2::huffman::IndexLoader::loadAccIndex(
				mergereq.children[0]->sortresult.getFiles().getBWT()
			);
			::libmaus2::huffman::IndexEntryContainerVector::unique_ptr_type IECV1 = ::libmaus2::huffman::IndexLoader::loadAccIndex(
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

					rl_decoder leftrlin(IDD0,IECV0.get(),ilow);
					rl_decoder rightrlin(IDD1,IECV1.get(),P[b]);
					
					uint64_t const outsuf = (ihigh-ilow)-(islast?1:0) + (P[b+1]-P[b]);

					rl_encoder bwtenc(encfilename,albits,outsuf,rlencoderblocksize);
				
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
			rl_encoder::concatenate(encfilenames,result.getFiles().getBWT());
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			result.setBWT(encfilenames);
			
			#if 0
			std::cerr << "[V] removing tmp files...";			
			rtc.start();
			for ( uint64_t i = 0; i < encfilenames.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile ( encfilenames[i].c_str() );
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
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					::libmaus2::util::TempFileRemovalContainer::addTempFile(newgapname);
					gapfilenames.push_back(newgapname);
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus2::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus2::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				libmaus2::aio::OutputStreamFactoryContainer::rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}

			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus2::lf::DArray::unique_ptr_type accD(new ::libmaus2::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
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
				libmaus2::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
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
					::libmaus2::util::TempFileRemovalContainer::addTempFile(renamed);
					libmaus2::aio::OutputStreamFactoryContainer::rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				// result.setGT(stringVectorAppend(GACR.gtpartnames,blockresults.getFiles().getGT()));
				
				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					libmaus2::aio::OutputStreamFactoryContainer::rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			// libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus2::timing::RealTimeClock mprtc;
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
				libmaus2::aio::FileRemoval::removeFile ( gapfilenames[i].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					libmaus2::aio::FileRemoval::removeFile ( bwtfilenames[i][j].c_str() );
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus2::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			libmaus2::wavelet::RlToHwtBase<true,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			libmaus2::wavelet::RlToHwtBase<false,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus2::util::TempFileRemovalContainer::addTempFile(result.getFiles().getHWTReq());
		{
		libmaus2::aio::OutputStreamInstance hwtreqCOS(result.getFiles().getHWTReq());
		libmaus2::wavelet::RlToHwtTermRequest::serialise(
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
		// hwtreqCOS.close();
		}

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
		libmaus2::aio::FileRemoval::removeFile ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus2::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
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
			::libmaus2::lf::DArray::unique_ptr_type accD(new ::libmaus2::lf::DArray(
				sblockhist
				)
			);
			// first block
			::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = 
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
				::libmaus2::util::TempFileRemovalContainer::addTempFile(renamed);
				libmaus2::aio::OutputStreamFactoryContainer::rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
			}

			result.setGT(stringVectorAppend(GACR.gtpartnames,oldgtnames));

			// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
			libmaus2::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap0dec(GACR.G->getDecoder());
			libmaus2::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap0decbuf(new libmaus2::suffixsort::GapArrayByteDecoderBuffer(*pgap0dec,8192));
			// libmaus2::suffixsort::GapArrayByteDecoderBuffer::iterator gap0decbufit = pgap0decbuf->begin();
			result.setBlockP0Rank( mergeIsa(
				mergereq.children[1]->sortresult.getFiles().getSampledISA(), // old sampled isa
				blockresults.getFiles().getSampledISA(), // new sampled isa
				result.getFiles().getSampledISA(),blockstart,*pgap0decbuf/*gap0decbufit*//*GACR.G.begin()*/,cblocksize+1 /* GACR.G.size() */
			) );
			pgap0decbuf.reset();
			pgap0dec.reset();
			
			::libmaus2::timing::RealTimeClock rtc; rtc.start();
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
			libmaus2::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap1dec(GACR.G->getDecoder());
			libmaus2::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap1decbuf(new libmaus2::suffixsort::GapArrayByteDecoderBuffer(*pgap1dec,8192));
			libmaus2::suffixsort::GapArrayByteDecoderBuffer::iterator gap1decbufit = pgap1decbuf->begin();
			::libmaus2::timing::RealTimeClock wprtc; wprtc.start();
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
				
				#if defined(GAP_ARRAY_BYTE_DEBUG)
				{
					// check obtained prefix sum
					libmaus2::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap2dec(GACR.G->getDecoder(ilow));
					libmaus2::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap2decbuf(new libmaus2::suffixsort::GapArrayByteDecoderBuffer(*pgap2dec,8192));
					libmaus2::suffixsort::GapArrayByteDecoderBuffer::iterator gap2decbufit = pgap2decbuf->begin();
					
					uint64_t a = 0;
					for ( uint64_t ia = ilow; ia < ihigh; ++ia )
						a += *(gap2decbufit++);
					
					assert ( p == a );
				}
				#endif

				P.push_back(P.back() + p);
				wpacks.push_back(std::pair<uint64_t,uint64_t>(ilow,ihigh));
				encfilenames.push_back(
					tmpfilenamebase 
					// result.getFiles().getBWT() 
					+ "_"
					+ ::libmaus2::util::NumberSerialisation::formatNumber(encfilenames.size(),6)
					+ ".bwt"
				);
				::libmaus2::util::TempFileRemovalContainer::addTempFile(encfilenames.back());
				ilow = ihigh;
			}
			assert ( wpacks.size() <= numthreads );
			// std::cerr << "done,time=" << wprtc.getElapsedSeconds() << ")";
			pgap1decbuf.reset();
			pgap1dec.reset();
			
			// std::cerr << "(setting up IDDs...";
			wprtc.start();
			unsigned int const albits = rl_decoder::haveAlphabetBits() ? rl_decoder::getAlBits(mergereq.children[0]->sortresult.getFiles().getBWT()) : 0;
			::libmaus2::huffman::IndexDecoderDataArray IDD0(
				mergereq.children[0]->sortresult.getFiles().getBWT());
			::libmaus2::huffman::IndexDecoderDataArray IDD1(
				mergereq.children[1]->sortresult.getFiles().getBWT());
			
			::libmaus2::huffman::IndexEntryContainerVector::unique_ptr_type IECV0 = ::libmaus2::huffman::IndexLoader::loadAccIndex(
				mergereq.children[0]->sortresult.getFiles().getBWT()
			);
			::libmaus2::huffman::IndexEntryContainerVector::unique_ptr_type IECV1 = ::libmaus2::huffman::IndexLoader::loadAccIndex(
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

					rl_decoder leftrlin(IDD0,IECV0.get(),ilow);
					rl_decoder rightrlin(IDD1,IECV1.get(),P[b]);
					
					uint64_t const outsuf = (ihigh-ilow)-(islast?1:0) + (P[b+1]-P[b]);

					rl_encoder bwtenc(encfilename,albits,outsuf,rlencoderblocksize);

					libmaus2::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap3dec(GACR.G->getDecoder(ilow));
					libmaus2::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap3decbuf(new libmaus2::suffixsort::GapArrayByteDecoderBuffer(*pgap3dec,8192));
					libmaus2::suffixsort::GapArrayByteDecoderBuffer::iterator gap3decbufit = pgap3decbuf->begin();
				
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
			rl_encoder::concatenate(encfilenames,result.getFiles().getBWT());
			std::cerr << "done, time " << rtc.getElapsedSeconds() << std::endl;
			#endif
			
			result.setBWT(encfilenames);
			
			#if 0
			std::cerr << "[V] removing tmp files...";			
			rtc.start();
			for ( uint64_t i = 0; i < encfilenames.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile ( encfilenames[i].c_str() );
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
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					::libmaus2::util::TempFileRemovalContainer::addTempFile(newgapname);
					gapfilenames.push_back(newgapname);
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus2::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus2::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				libmaus2::aio::OutputStreamFactoryContainer::rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}

			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus2::lf::DArray::unique_ptr_type accD(new ::libmaus2::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
				std::string const newmergedgapoverflow = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".gapoverflow";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedgapoverflow);
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
				// libmaus2::util::GetObject<uint32_t const *> mergeGO(GACR.G.begin());
				libmaus2::suffixsort::GapArrayByteDecoder::unique_ptr_type pgap0dec(GACR.G->getDecoder());
				libmaus2::suffixsort::GapArrayByteDecoderBuffer::unique_ptr_type pgap0decbuf(new libmaus2::suffixsort::GapArrayByteDecoderBuffer(*pgap0dec,8192));
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
					::libmaus2::util::TempFileRemovalContainer::addTempFile(renamed);
					libmaus2::aio::OutputStreamFactoryContainer::rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				// result.setGT(stringVectorAppend(GACR.gtpartnames,blockresults.getFiles().getGT()));
				
				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					libmaus2::aio::OutputStreamFactoryContainer::rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			// libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus2::timing::RealTimeClock mprtc;
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
				libmaus2::aio::FileRemoval::removeFile ( gapfilenames[i].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					libmaus2::aio::FileRemoval::removeFile ( bwtfilenames[i][j].c_str() );
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus2::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			libmaus2::wavelet::RlToHwtBase<true,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			libmaus2::wavelet::RlToHwtBase<false,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus2::util::TempFileRemovalContainer::addTempFile(result.getFiles().getHWTReq());
		{
		libmaus2::aio::OutputStreamInstance hwtreqCOS(result.getFiles().getHWTReq());
		libmaus2::wavelet::RlToHwtTermRequest::serialise(
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
		//hwtreqCOS.close();
		}

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
		std::string const sparsetmpfilenamebase,
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
		libmaus2::aio::FileRemoval::removeFile ( mergereq.children[mergereq.children.size()-1]->sortresult.getFiles().getHWT().c_str() );
		
		// get result object
		::libmaus2::suffixsort::BwtMergeBlockSortResult & result = mergereq.sortresult;
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
					std::string const newgapname = tmpfilenamebase + "_merging_" + ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) + ".gap";
					gapfilenameprefixes.push_back(newgapname);	
					gapfilenames.push_back(std::vector<std::string>());
				}

				// bwt name
				std::vector<std::string> newbwtnames;
				for ( uint64_t i = 0; i < mergereq.children[bb]->sortresult.getFiles().getBWT().size(); ++i )
				{
					std::string const newbwtname = tmpfilenamebase + "_merging_" 
						+ ::libmaus2::util::NumberSerialisation::formatNumber(bb,4) 
						+ "_"
						+ ::libmaus2::util::NumberSerialisation::formatNumber(i,4) 
						+ ".bwt";
					::libmaus2::util::TempFileRemovalContainer::addTempFile(newbwtname);
					newbwtnames.push_back(newbwtname);
				}
				bwtfilenames.push_back(newbwtnames);
			}

			// rename last bwt file set
			for ( uint64_t i = 0; i < mergereq.children.back()->sortresult.getFiles().getBWT().size(); ++i )
			{
				libmaus2::aio::OutputStreamFactoryContainer::rename ( 
					mergereq.children.back()->sortresult.getFiles().getBWT()[i].c_str(),
					bwtfilenames.back()[i].c_str() 
				);
			}


			std::vector<std::string> mergedgtname  = mergereq.children.back()->sortresult.getFiles().getGT();
			std::string mergedisaname = mergereq.children.back()->sortresult.getFiles().getSampledISA();

			// load char histogram for last block
			std::string const & lblockhist = mergereq.children.back()->sortresult.getFiles().getHist();
			::libmaus2::lf::DArray::unique_ptr_type accD(new ::libmaus2::lf::DArray(lblockhist));

			/**
			 * iteratively merge blocks together
			 **/
			for ( uint64_t bb = 0; bb+1 < mergereq.children.size(); ++bb )
			{
				// std::cerr << "** WHITEBOX EXTERNAL **" << std::endl;
				
				// block we merge into
				uint64_t const bx = mergereq.children.size()-bb-2;
				std::cerr << "[V] merging blocks " << bx+1 << " to end into " << bx << std::endl;
				::libmaus2::suffixsort::BwtMergeBlockSortResult const & blockresults = 
					mergereq.children[bx]->sortresult;

				// output files for this iteration
				std::string const newmergedgtname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".gt";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedgtname);
				std::string const newmergedisaname = tmpfilenamebase + "_merged_" + ::libmaus2::util::NumberSerialisation::formatNumber(bx,4) + ".sampledisa";
				::libmaus2::util::TempFileRemovalContainer::addTempFile(newmergedisaname);
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
					gapfilenameprefix,
					sparsetmpfilenamebase+"_sparsegap",
					mem
				);
				gapfilenames[bx] = GACR.fn;

				// merge sampled inverse suffix arrays, returns rank of position 0 (relative to block start)
				libmaus2::gamma::GammaGapDecoder GGD(gapfilenames[bx]);
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
					::libmaus2::util::TempFileRemovalContainer::addTempFile(renamed);
					libmaus2::aio::OutputStreamFactoryContainer::rename(blockresults.getFiles().getGT()[i].c_str(), renamed.c_str());
				}

				/*
				 * remove files we no longer need
				 */
				// files local to this block
				for ( uint64_t i = 0; i < blockresults.getFiles().getBWT().size(); ++i )
					libmaus2::aio::OutputStreamFactoryContainer::rename ( blockresults.getFiles().getBWT()[i].c_str(), bwtfilenames[bx][i].c_str() );
				blockresults.removeFilesButBwt();
				// previous stage gt bit vector
				for ( uint64_t i = 0; i < mergedgtname.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile ( mergedgtname[i].c_str() );
				
				// update current file names
				mergedgtname = stringVectorAppend(GACR.gtpartnames,oldgtnames);
				mergedisaname = newmergedisaname;
			}
			
			// renamed sampled inverse suffix array
			libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedisaname.c_str(), result.getFiles().getSampledISA().c_str() );
			// rename gt bit array filename
			// libmaus2::aio::OutputStreamFactoryContainer::rename ( mergedgtname.c_str(), result.getFiles().getGT().c_str() );
			result.setGT(mergedgtname);
			// save histogram
			accD->serialise(static_cast<std::string const & >(result.getFiles().getHist()));

			std::cerr << "[V] merging parts...";
			::libmaus2::timing::RealTimeClock mprtc;
			mprtc.start();
			result.setBWT(
				parallelGapFragMerge(bwtfilenames,gapfilenames,tmpfilenamebase+"_gpart",numthreads,
					lfblockmult,rlencoderblocksize)
			);
			std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;

			for ( uint64_t i = 0; i < gapfilenames.size(); ++i )
				for ( uint64_t j = 0; j < gapfilenames[i].size(); ++j )
					libmaus2::aio::FileRemoval::removeFile ( gapfilenames[i][j].c_str() );
			for ( uint64_t i = 0; i < bwtfilenames.size(); ++i )
				for ( uint64_t j = 0; j < bwtfilenames[i].size(); ++j )
					libmaus2::aio::FileRemoval::removeFile ( bwtfilenames[i][j].c_str() );		
		}

		#if 0
		std::cerr << "[V] computing term symbol hwt...";
		::libmaus2::timing::RealTimeClock mprtc;
		mprtc.start();
		if ( input_types_type::utf8Wavelet() )
			libmaus2::wavelet::RlToHwtBase<true,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		else
			libmaus2::wavelet::RlToHwtBase<false,rl_decoder>::rlToHwtTerm(result.getFiles().getBWT(),result.getFiles().getHWT(),tmpfilenamebase + "_wt",chist,bwtterm,result.getBlockP0Rank());
		std::cerr << "done, time " << mprtc.getElapsedSeconds() << std::endl;
		#endif

		libmaus2::util::TempFileRemovalContainer::addTempFile(result.getFiles().getHWTReq());
		{
		libmaus2::aio::OutputStreamInstance hwtreqCOS(result.getFiles().getHWTReq());
		libmaus2::wavelet::RlToHwtTermRequest::serialise(hwtreqCOS,
			result.getFiles().getBWT(),
			result.getFiles().getHWT(),
			tmpfilenamebase + "_wt",
			huftreefilename,
			bwtterm,
			result.getBlockP0Rank(),
			input_types_type::utf8Wavelet()
		);
		hwtreqCOS.flush();
		//hwtreqCOS.close();
		}

		// remove obsolete files
		for ( uint64_t b = 0; b < mergereq.children.size(); ++b )
			mergereq.children[b]->sortresult.removeFiles();

		mergereq.releaseChildren();
	}
	
	static void sortIsaFile(std::string const & mergedisaname, uint64_t const blockmem)
	{
		// sort sampled inverse suffix array file
		std::string const mergeisatmp = mergedisaname+".tmp";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(mergeisatmp);
		std::string const mergeisatmpout = mergedisaname+".tmp.out";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(mergeisatmpout);
		// uint64_t const blockmem = 5*blocksize;
		// uint64_t const blockels = (blockmem + 2*sizeof(uint64_t)-1)/(2*sizeof(uint64_t));
		::libmaus2::sorting::PairFileSorting::sortPairFile(
			std::vector<std::string>(1,mergedisaname),mergeisatmp,true /* second comp */,
			true,true,mergeisatmpout,blockmem/2/*par*/,true /* parallel */);
		libmaus2::aio::FileRemoval::removeFile ( (mergeisatmp).c_str() );
		libmaus2::aio::OutputStreamFactoryContainer::rename ( mergeisatmpout.c_str(), mergedisaname.c_str() );
	
	}

	static uint64_t readBlockRanksSize(std::string const & mergedisaname)
	{
		return ::libmaus2::util::GetFileSize::getFileSize(mergedisaname)/(2*sizeof(uint64_t));		
	}

	static ::libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > readBlockRanks(std::string const & mergedisaname)
	{
		// read sampled isa
		uint64_t const nsisa = readBlockRanksSize(mergedisaname);
		::libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > blockranks(nsisa,false);
		::libmaus2::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type SGIsisa(new ::libmaus2::aio::SynchronousGenericInput<uint64_t>(mergedisaname,16*1024));
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
		::libmaus2::parallel::OMPLock cerrlock;
		// number of sampled suffix array elements
		uint64_t const nsa = (fs + sasamplingrate - 1) / sasamplingrate;
		
		// check that this matches what we have in the file
		assert ( ::libmaus2::util::GetFileSize::getFileSize(mergedsaname) / (sizeof(uint64_t)) ==  nsa + 2 );
		
		if ( nsa && nsa-1 )
		{
			uint64_t const checkpos = nsa-1;
			uint64_t const satcheckpacks = numthreads * lfblockmult;
			uint64_t const sacheckblocksize = (checkpos + satcheckpacks-1) / satcheckpacks;
			uint64_t const sacheckpacks = ( checkpos + sacheckblocksize - 1 ) / sacheckblocksize;
			
			std::cerr << "[V] checking suffix array on text...";
			::libmaus2::parallel::SynchronousCounter<uint64_t> SC;
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

				::libmaus2::aio::SynchronousGenericInput<uint64_t> SGIsa(mergedsaname,16*1024,low+2,cnt+1);
				
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
		::libmaus2::lf::ImpCompactHuffmanWaveletLF const & IHWT,
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

		::libmaus2::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type SGIsisasa(new ::libmaus2::aio::SynchronousGenericInput<uint64_t>(mergedisaname,16*1024));
		int64_t const fr = SGIsisasa->get(); assert ( fr != -1 );
		int64_t const fp = SGIsisasa->get(); assert ( fp != -1 );
		
		std::cerr << "[V] computing sampled suffix array parts...";
	
		std::pair<uint64_t,uint64_t> const isa0(fr,fp);
		std::pair<uint64_t,uint64_t> isapre(isa0);
		
		::std::vector< std::string > satempfilenames(numthreads);
		::std::vector< std::string > isatempfilenames(numthreads);
		::libmaus2::autoarray::AutoArray < ::libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type > SAF(numthreads);
		::libmaus2::autoarray::AutoArray < ::libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type > ISAF(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			satempfilenames[i] = ( tmpfilenamebase + ".sampledsa_" + ::libmaus2::util::NumberSerialisation::formatNumber(i,6) );
			::libmaus2::util::TempFileRemovalContainer::addTempFile(satempfilenames[i]);
			::libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type tSAFi(
				new ::libmaus2::aio::SynchronousGenericOutput<uint64_t>(satempfilenames[i],8*1024)
			);
			SAF[i] = UNIQUE_PTR_MOVE(tSAFi);

			isatempfilenames[i] = ( tmpfilenamebase + ".sampledisa_" + ::libmaus2::util::NumberSerialisation::formatNumber(i,6) );
			::libmaus2::util::TempFileRemovalContainer::addTempFile(isatempfilenames[i]);
			::libmaus2::aio::SynchronousGenericOutput<uint64_t>::unique_ptr_type tISAFi(
				new ::libmaus2::aio::SynchronousGenericOutput<uint64_t>(isatempfilenames[i],8*1024)
			);
			ISAF[i] = UNIQUE_PTR_MOVE(tISAFi);
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
		std::string const mergedsaname = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".sa";
		{
		::libmaus2::aio::OutputStreamInstance::unique_ptr_type pmergedsa(new ::libmaus2::aio::OutputStreamInstance(mergedsaname));
		// write sampling rate
		::libmaus2::serialize::Serialize<uint64_t>::serialize(*pmergedsa,sasamplingrate);
		::libmaus2::serialize::Serialize<uint64_t>::serialize(*pmergedsa,(fs + sasamplingrate-1)/sasamplingrate);
		std::string const mergesatmp = mergedsaname + ".tmp";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(mergesatmp);
		::libmaus2::sorting::PairFileSorting::sortPairFile(
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
		libmaus2::aio::FileRemoval::removeFile(mergesatmp.c_str());
		}
		std::cerr << "done." << std::endl;		

		std::cerr << "[V] sorting and merging sampled inverse suffix array parts...";
		std::string const mergedisaoutname = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".isa";
		::libmaus2::aio::OutputStreamInstance::unique_ptr_type pmergedisa(new ::libmaus2::aio::OutputStreamInstance(mergedisaoutname));
		// write sampling rate
		::libmaus2::serialize::Serialize<uint64_t>::serialize(*pmergedisa,isasamplingrate);
		::libmaus2::serialize::Serialize<uint64_t>::serialize(*pmergedisa,(fs+isasamplingrate-1)/isasamplingrate);
		std::string const mergeisatmp = mergedisaoutname + ".tmp";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(mergeisatmp);
		::libmaus2::sorting::PairFileSorting::sortPairFile(
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
		libmaus2::aio::FileRemoval::removeFile(mergeisatmp.c_str());

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
	
	static bool getDefaultBWTOnly()
	{
		return false;
	}

	static uint64_t getDefaultBlockSize(uint64_t const mem, uint64_t const threads, uint64_t const fs)
	{
		uint64_t const memblocksize = std::max(static_cast<uint64_t>(0.95 * mem / ( 5 * threads )),static_cast<uint64_t>(1));
		uint64_t const fsblocksize = (fs + threads - 1)/threads;
		return std::min(memblocksize,fsblocksize);
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
		libmaus2::autoarray::AutoArray<uint64_t> L;
		libmaus2::util::ExtendingSimpleCountingHash<uint64_t,uint64_t> H;
		libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > P;
		uint64_t p;
		
		HashHistogram(uint64_t const lowsize = 256, uint64_t const bigsize = (1ull<<16) ) : L(lowsize), H(libmaus2::math::nextTwoPow(bigsize)), P(64ull*1024ull), p(0) {}
		
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
			for ( libmaus2::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::key_type const * k = H.begin(); k != H.end(); ++k )
				if ( *k != libmaus2::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::unused() )
					nonzero++;
					
			if ( P.size() < nonzero )
				P.resize(nonzero);
				
			p = 0;

			for ( uint64_t i = 0; i < L.size(); ++i )
				if ( L[i] )
					P[p++] = std::pair<uint64_t,uint64_t>(i,L[i]);
					
			for ( libmaus2::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::key_type const * k = H.begin(); k != H.end(); ++k )
				if ( *k != libmaus2::util::ExtendingSimpleCountingHash<uint64_t,uint64_t>::unused() )
					P[p++] = std::pair<uint64_t,uint64_t>(*k,H.cntbegin() [ (k-H.begin()) ] );
					
			std::sort(P.begin(),P.begin()+p);
		}
	};

	static void getBlockSymFreqsHash(std::string const fn, uint64_t const glow, uint64_t ghigh, HashHistogram & H)
	{
		typedef typename input_types_type::linear_wrapper stream_type;
		typedef typename input_types_type::base_input_stream::char_type char_type;
		typedef typename ::libmaus2::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
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

		::libmaus2::autoarray::AutoArray<unsigned_char_type> GB(loopsize*numthreads,false);

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
		typedef typename ::libmaus2::util::UnsignedCharVariant<char_type>::type unsigned_char_type;
		
		uint64_t const fs = ghigh-glow;
		
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		uint64_t const symsperfrag = (fs + numthreads - 1)/numthreads;
		uint64_t const numfrags = (fs + symsperfrag - 1)/symsperfrag;

		libmaus2::util::HistogramSet HS(numfrags,256);

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

			::libmaus2::autoarray::AutoArray<unsigned_char_type> B(16*1024,false);
			libmaus2::util::Histogram & H = HS[t];
			
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

		::libmaus2::util::Histogram::unique_ptr_type PH(HS.merge());
		::libmaus2::util::Histogram & H(*PH);

		return H.getByType<int64_t>();
	}
	
	static uint64_t getBlockStart(uint64_t const b, uint64_t const blocksize, uint64_t const fullblocks)
	{
		return (b < fullblocks) ? (b*blocksize) : (fullblocks * blocksize + (b-fullblocks) * (blocksize-1));
	}
	
	static uint64_t getBlockSize(uint64_t const b, uint64_t const blocksize, uint64_t const fullblocks)
	{
		return (b < fullblocks) ? blocksize : (blocksize-1);
	}
				
	static int computeBwt(::libmaus2::util::ArgInfo const & arginfo)
	{
		libmaus2::timing::RealTimeClock bwtclock;
		bwtclock.start();
	
		#if defined(BWTB3M_MEMORY_DEBUG)
		uint64_t mcnt = 0;
		#endif
		
		::libmaus2::util::TempFileRemovalContainer::setup();
		uint64_t const rlencoderblocksize = 16*1024;
		::libmaus2::parallel::OMPLock cerrlock;

		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else	
		uint64_t const numthreads = 1;
		#endif

		// compute BWT only? (no SA and ISA)
		bool const bwtonly = arginfo.getValue<unsigned int>("bwtonly",getDefaultBWTOnly());
		// total memory available
		uint64_t const mem = std::max(static_cast<uint64_t>(1),arginfo.getValueUnsignedNumeric<uint64_t>("mem",getDefaultMem()));
		// base for tmp file names
		std::string const tmpfilenamebase = arginfo.getUnparsedValue("tmpprefix",arginfo.getDefaultTmpFileName());
		// base for spare tmp file names
		std::string const sparsetmpfilenamebase = arginfo.getUnparsedValue("sparsetmpprefix",tmpfilenamebase);
		// file name of serialised character histogram
		std::string const chistfilename = tmpfilenamebase + ".chist";
		// file name of serialised huffman tree
		std::string const huftreefilename = tmpfilenamebase + ".huftree";
		::libmaus2::util::TempFileRemovalContainer::addTempFile(chistfilename);
		::libmaus2::util::TempFileRemovalContainer::addTempFile(huftreefilename);
		// output file name
		std::string const outfn = arginfo.getValue<std::string>("outputfilename",tmpfilenamebase+".bwt");
		// final inverse suffix array sampling rate
		uint64_t const isasamplingrate = ::libmaus2::math::nextTwoPow(arginfo.getValue<uint64_t>("isasamplingrate",getDefaultIsaSamplingRate()));
		// final suffix array sampling rate
		uint64_t const sasamplingrate = ::libmaus2::math::nextTwoPow(arginfo.getValue<uint64_t>("sasamplingrate",getDefaultSaSamplingRate()));

		// file name		
		std::string const fn = arginfo.getRestArg<std::string>(0);
		// check whether file exists
		if ( ! ::libmaus2::util::GetFileSize::fileExists(fn) )
		{
			::libmaus2::exception::LibMausException se;
			se.getStream() << "File " << fn << " does not exist or cannot be opened." << std::endl;
			se.finish();
			throw se;
		}

		/* get file size */
		uint64_t const fs = input_types_type::linear_wrapper::getFileSize(fn);

		/* check that file is not empty */
		if ( ! fs )
		{
			::libmaus2::exception::LibMausException se;
			se.getStream() << "File " << fn << " is empty." << std::endl;
			se.finish();
			throw se;		
		}

		// target block size
		uint64_t const tblocksize = 
			std::max(
				static_cast<uint64_t>(1),
				arginfo.getValueUnsignedNumeric<uint64_t>(
					"blocksize",getDefaultBlockSize(mem,numthreads,fs)
				)
		);
		// number of blocks
		uint64_t const numblocks = (fs + tblocksize - 1) / tblocksize;
		// final block size
		uint64_t const blocksize = (fs + numblocks - 1)/ numblocks;
		// full block product
		uint64_t const fullblockprod = numblocks * blocksize;
		// extraneous
		uint64_t const extrasyms = fullblockprod - fs;
		// check
		assert ( extrasyms < numblocks );
		// full blocks
		uint64_t const fullblocks = numblocks - extrasyms;
		// reduced blocks
		uint64_t const redblocks = numblocks - fullblocks;
		// check
		assert ( fullblocks * blocksize + redblocks * (blocksize-1) == fs );
		
		// next power of two
		uint64_t const blocksizenexttwo = ::libmaus2::math::nextTwoPow(blocksize);
		// prev power of two
		uint64_t const blocksizeprevtwo = (blocksize == blocksizenexttwo) ? blocksize : (blocksizenexttwo / 2);
		
		// ISA sampling rate during block merging
		uint64_t const preisasamplingrate = std::min(::libmaus2::math::nextTwoPow(arginfo.getValue<uint64_t>("preisasamplingrate",256*1024)),blocksizeprevtwo);

			
		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		std::cerr << "[V] sorting file " << fn << " of size " << fs << " with block size " << blocksize << " (" << numblocks << " blocks)" << " and " << numthreads << " threads" << std::endl;
		std::cerr << "[V] full blocks " << fullblocks << " reduced blocks " << redblocks << std::endl;
		
		// there should be at least one block
		assert ( numblocks );

		#if 0
		std::vector<int64_t> minblockperiods(numblocks);
		// compute periods
		libmaus2::parallel::OMPLock block;
		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			// block id
			uint64_t const b = numblocks-bb-1;
			// start of block in file
			uint64_t const blockstart = getBlockStart(b,blocksize,fullblocks);

			typedef typename input_types_type::base_input_stream base_input_stream;
			typedef typename base_input_stream::char_type char_type;
			typedef typename input_types_type::circular_wrapper circular_wrapper;

			uint64_t const readlen = 3 * blocksize;
			circular_wrapper textstr(fn,blockstart);
			libmaus2::autoarray::AutoArray<char_type> A(readlen,false);
			textstr.read(&A[0],A.size());

			block.lock();
			std::cerr << "\r" << std::string(80,' ') << "\r[Checking " << (bb+1) << "/" << numblocks << "]\r";
			block.unlock();

			libmaus2::util::BorderArray<uint32_t> SBA(A.begin(),A.size());
			
			uint64_t minper = std::numeric_limits<uint64_t>::max();
			for ( uint64_t i = (blocksize-1); i < A.size(); ++i )
			{
				uint64_t const period = (i+1)-SBA[i];

				// if period divides length of the prefix
				if ( 
					( ( (i+1) % period ) == 0 ) && ((i+1) / period > 1) 
					&&
					(period < minper)
				)
				{
					block.lock();
					std::cerr << "\nlen " << (i+1) << " block " << b << " period " << period << std::endl;
					block.unlock();
					minper = period;
				}
			}
			
			if ( minper == std::numeric_limits<uint64_t>::max() )
				minblockperiods[b] = -1;
			else
				minblockperiods[b] = minper;
				
		}
		std::cerr << std::endl;
		
		exit(0);
		#endif

		std::cerr << "[V] computing LCP between block suffixes and the following block start: ";
		std::vector<uint64_t> largelcpblocks;
		libmaus2::parallel::OMPLock largelcpblockslock;
		uint64_t const largelcpthres = 16*1024;
		std::vector<uint64_t> boundedlcpblockvalues(numblocks);
		libmaus2::parallel::SynchronousCounter<uint64_t> largelcpblockscomputed(0);
		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			// block id
			uint64_t const b = numblocks-bb-1;
			// start of block in file
			uint64_t const blockstart = getBlockStart(b,blocksize,fullblocks);
			// size of this block
			uint64_t const cblocksize = getBlockSize(b,blocksize,fullblocks);

			// start of next block
			uint64_t const nextblockstart = (blockstart + cblocksize) % fs;
		
			// find bounded lcp between this block and start of next
			uint64_t const blcp = libmaus2::suffixsort::BwtMergeBlockSortRequestBase::findSplitCommonBounded<input_types_type>(fn,blockstart,cblocksize,nextblockstart,fs,largelcpthres);

			if ( blcp >= largelcpthres )
			{
				libmaus2::parallel::ScopeLock slock(largelcpblockslock);
				largelcpblocks.push_back(b);
			}

			{
			libmaus2::parallel::ScopeLock slock(largelcpblockslock);
			uint64_t const finished = ++largelcpblockscomputed;
			std::cerr << "(" << static_cast<double>(finished)/numblocks << ")";
			}
			
			boundedlcpblockvalues[b] = blcp;
		}
		std::cerr << "done." << std::endl;
		
		std::sort(largelcpblocks.begin(),largelcpblocks.end());
		for ( uint64_t ib = 0; ib < largelcpblocks.size(); ++ib )
		{
			uint64_t const b = largelcpblocks[ib];
			
			std::cerr << "[V] Recomputing lcp value for block " << b << std::endl;

			// start of block in file
			uint64_t const blockstart = getBlockStart(b,blocksize,fullblocks);
			// size of this block
			uint64_t const cblocksize = getBlockSize(b,blocksize,fullblocks);

			// start of next block
			uint64_t const nextblockstart = (blockstart + cblocksize) % fs;
		
			// find bounded lcp between this block and start of next
			uint64_t const blcp = libmaus2::suffixsort::BwtMergeBlockSortRequestBase::findSplitCommon<input_types_type>(fn,blockstart,cblocksize,nextblockstart,fs);

			boundedlcpblockvalues[b] = blcp;
		}
		
		// exit(0);

		::libmaus2::suffixsort::BwtMergeTempFileNameSetVector blocktmpnames(tmpfilenamebase, numblocks, numthreads /* bwt */, numthreads /* gt */);

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		std::cerr << "[V] computing symbol frequences" << std::endl;
		std::map<int64_t,uint64_t> chistnoterm;	
		// std::vector< std::map<int64_t,uint64_t> > blockfreqvec(numblocks);
		libmaus2::timing::RealTimeClock rtc;
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			// block id
			uint64_t const b = numblocks-bb-1;
			// start of block in file
			uint64_t const blockstart = getBlockStart(b,blocksize,fullblocks);
			// size of this block
			uint64_t const cblocksize = getBlockSize(b,blocksize,fullblocks);

			std::map<int64_t,uint64_t> const blockfreqs =
				getBlockSymFreqs(fn,blockstart,blockstart+cblocksize);
			
			std::string const freqstmpfilename = blocktmpnames[b].getHist() + ".freqs";
			libmaus2::util::TempFileRemovalContainer::addTempFile(freqstmpfilename);	
			{
			libmaus2::aio::OutputStreamInstance freqCOS(freqstmpfilename);
			libmaus2::util::NumberMapSerialisation::serialiseMap(freqCOS,blockfreqs);
			freqCOS.flush();
			//freqCOS.close();
			}
			
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

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		libmaus2::aio::OutputStreamInstance::unique_ptr_type chistCOS(new libmaus2::aio::OutputStreamInstance(chistfilename));
		(*chistCOS) << ::libmaus2::util::NumberMapSerialisation::serialiseMap(chist);
		chistCOS->flush();
		// chistCOS->close();
		chistCOS.reset();
		
		std::cerr << "[V] computed symbol frequences, input alphabet size is " << chistnoterm.size() << std::endl;

		std::cerr << "[V] bwtterm=" << bwtterm << std::endl;

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		libmaus2::huffman::HuffmanTree::unique_ptr_type uhnode(new libmaus2::huffman::HuffmanTree(chist.begin(),chist.size(),false,true,true));
		
		libmaus2::aio::OutputStreamInstance::unique_ptr_type huftreeCOS(new libmaus2::aio::OutputStreamInstance(huftreefilename));
		uhnode->serialise(*huftreeCOS);
		huftreeCOS->flush();
		// huftreeCOS->close();
		huftreeCOS.reset();

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		::libmaus2::huffman::HuffmanTree::EncodeTable::unique_ptr_type EC(
			new ::libmaus2::huffman::HuffmanTree::EncodeTable(*uhnode));

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		std::vector < MergeStrategyBlock::shared_ptr_type > stratleafs(numblocks);

		bool const computeTermSymbolHwt = arginfo.getValue<int>("computeTermSymbolHwt",false);;
	
		for ( uint64_t bb = 0; bb < numblocks; ++bb )
		{
			uint64_t const b = numblocks-bb-1;

			// start of block in file
			uint64_t const blockstart = getBlockStart(b,blocksize,fullblocks);
			// size of this block
			uint64_t const cblocksize = getBlockSize(b,blocksize,fullblocks);

			// symbol frequency map
			// std::map<int64_t,uint64_t> const & blockfreqs = blockfreqvec[b];
			std::string const freqstmpfilename = blocktmpnames[b].getHist() + ".freqs";
			libmaus2::aio::InputStreamInstance freqCIS(freqstmpfilename);
			std::map<int64_t,uint64_t> const blockfreqs = 
				libmaus2::util::NumberMapSerialisation::deserialiseMap<libmaus2::aio::InputStreamInstance,int64_t,uint64_t>(freqCIS);
			// freqCIS.close();
			libmaus2::aio::FileRemoval::removeFile(freqstmpfilename.c_str());
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
			::libmaus2::suffixsort::BwtMergeZBlockRequestVector zreqvec;
			dynamic_cast<MergeStrategyBaseBlock *>(PMSB.get())->sortreq = 
				libmaus2::suffixsort::BwtMergeBlockSortRequest(
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
					computeTermSymbolHwt,
					boundedlcpblockvalues[b]
				);

			// std::cerr << *PMSB;
			
			stratleafs[b] = PMSB;

			#if defined(BWTB3M_MEMORY_DEBUG)
			std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
			#endif
		}

		uhnode.reset();
		EC.reset();

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif


		// word available per computation thread
		uint64_t const wordsperthread = std::max(static_cast<uint64_t>(1),arginfo.getValueUnsignedNumeric<uint64_t>("wordsperthread",getDefaultWordsPerThread()));

		std::cerr << "[V] constructing merge tree" << std::endl;

		// construct merge tree and register z requests
		MergeStrategyBlock::shared_ptr_type mergetree = constructMergeTree(stratleafs,mem,numthreads,wordsperthread);
		
		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		// inner node queue
		std::deque<MergeStrategyBlock *> itodo;

		std::cerr << "[V] sorting blocks" << std::endl;

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif
		
		// sort single blocks
		BaseBlockSorting::unique_ptr_type BBS(new BaseBlockSorting(stratleafs,mem,numthreads,itodo));
		BBS->start();
		BBS->join();
		BBS.reset();

		std::cerr << "[V] sorted blocks" << std::endl;
		
		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
		#endif

		#if 0
		uint64_t maxhwtsize = 0;
		for ( uint64_t i = 0; i < stratleafs.size(); ++i )
			maxhwtsize = std::max(maxhwtsize,::libmaus2::util::GetFileSize::getFileSize(stratleafs[i]->sortresult.getFiles().getHWT()));
		#endif
		
		std::cerr << "[V] filling gap request objects" << std::endl;
		mergetree->fillGapRequestObjects(numthreads);

		#if defined(BWTB3M_MEMORY_DEBUG)
		std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
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
			std::ostringstream sparsetmpstr;
			sparsetmpstr << sparsetmpfilenamebase << "_" << std::setfill('0') << std::setw(6) << (mtmpid++);
			
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
					sparsetmpstr.str(),
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
				
			#if defined(BWTB3M_MEMORY_DEBUG)
			std::cerr << "[M"<< (mcnt++) << "] " << libmaus2::util::MemUsage() << " " << libmaus2::autoarray::AutoArrayMemUsage() << std::endl;
			#endif
		}

		#if 0
		uint64_t const memperthread = maxblocksize*sizeof(uint32_t) + maxhwtsize;
		#endif

		uint64_t const memperthread = (mem + numthreads-1)/numthreads;

		::libmaus2::suffixsort::BwtMergeBlockSortResult const mergeresult = mergetree->sortresult;
		
		rl_encoder::concatenate(mergeresult.getFiles().getBWT(),outfn,true /* removeinput */);
		// libmaus2::aio::OutputStreamFactoryContainer::rename ( mergeresult.getFiles().getBWT().c_str(), outfn.c_str() );
		
		std::cerr << "[V] BWT computed in time " << bwtclock.formatTime(bwtclock.getElapsedSeconds()) << std::endl;

		// serialise character histogram
		std::string const outhist = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hist";
		::libmaus2::aio::OutputStreamInstance::unique_ptr_type Phistout(new ::libmaus2::aio::OutputStreamInstance(outhist));
		::libmaus2::util::NumberMapSerialisation::serialiseMap(*Phistout,chistnoterm);
		Phistout->flush();
		// Phistout->close();
		Phistout.reset();
		
		// remove hwt request for term symbol hwt
		libmaus2::aio::FileRemoval::removeFile ( mergeresult.getFiles().getHWTReq().c_str() );
		
		// remove hwt (null op)
		std::string const debhwt = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hwt" + ".deb";
		if ( libmaus2::aio::InputStreamFactoryContainer::tryOpen(mergeresult.getFiles().getHWT()) )
		{
			libmaus2::aio::OutputStreamFactoryContainer::rename ( mergeresult.getFiles().getHWT().c_str(), debhwt.c_str() );
			libmaus2::aio::FileRemoval::removeFile ( debhwt.c_str() );
		}
		
		// remove gt files
		for ( uint64_t i = 0; i < mergeresult.getFiles().getGT().size(); ++i )
			libmaus2::aio::FileRemoval::removeFile ( mergeresult.getFiles().getGT()[i].c_str() );
		
		if ( bwtonly )
		{
			std::string const mergedisaname = mergeresult.getFiles().getSampledISA();
			std::string const outisa = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".preisa";
			libmaus2::aio::OutputStreamFactoryContainer::rename(mergedisaname.c_str(),outisa.c_str());
		}
		else
		{	
			std::cerr << "[V] computing Huffman shaped wavelet tree of final BWT...";	
			std::string const outhwt = ::libmaus2::util::OutputFileNameTools::clipOff(outfn,".bwt") + ".hwt";
			libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type pICHWT;
			if ( input_types_type::utf8Wavelet() )
			{
				libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(
					libmaus2::wavelet::RlToHwtBase<true,rl_decoder>::rlToHwt(outfn, outhwt, tmpfilenamebase+"_finalhwttmp")
				);
				pICHWT = UNIQUE_PTR_MOVE(tICHWT);
			}
			else
			{
				libmaus2::wavelet::ImpCompactHuffmanWaveletTree::unique_ptr_type tICHWT(
					libmaus2::wavelet::RlToHwtBase<false,rl_decoder>::rlToHwt(outfn, outhwt, tmpfilenamebase+"_finalhwttmp")
				);		
				pICHWT = UNIQUE_PTR_MOVE(tICHWT);
			}
			std::cerr << "done, " << std::endl;
			
			std::cerr << "[V] loading Huffman shaped wavelet tree of final BWT...";	
			::libmaus2::lf::ImpCompactHuffmanWaveletLF IHWT(pICHWT);
			std::cerr << "done." << std::endl;

			// sort the sampled isa file	
			uint64_t const blockmem = memperthread; // memory per thread
			std::string const mergedisaname = mergeresult.getFiles().getSampledISA();
			sortIsaFile(mergedisaname,blockmem);

			// compute sampled suffix array and sampled inverse suffix array
			computeSampledSA(
				fn,fs,IHWT,mergedisaname,outfn,tmpfilenamebase,
				numthreads,lfblockmult,sasamplingrate,isasamplingrate,blockmem
			);

			// std::cerr << "[V] mergeresult.blockp0rank=" << mergeresult.blockp0rank << std::endl;
		}
		
		return EXIT_SUCCESS;
	}
};

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::timing::RealTimeClock rtc; rtc.start();

		#if defined(_OPENMP)
		unsigned int const maxthreads = omp_get_max_threads();
		unsigned int const numthreads = arginfo.getValue<unsigned int>("numthreads", BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultNumThreads());
		omp_set_num_threads(numthreads);
		#endif
		
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
			str << "sasamplingrate=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultSaSamplingRate() << "] sampling rate for sampled suffix array"<< std::endl;
			str << "isasamplingrate=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultIsaSamplingRate() << "] sampling rate for sampled inverse suffix array"<< std::endl;
			// str << "blocksize=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultBlockSize(BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultMem(),numthreads) << "] block size" << std::endl;
			str << "mem=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultMem() << "] memory target" << std::endl;
			#if defined(_OPENMP)
			str << "numthreads=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultNumThreads() << "] number of threads" << std::endl;
			#endif
			str << "bwtonly=[" << BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::getDefaultBWTOnly() << "] compute BWT only (no sampled suffix array and reverse)" << std::endl;
			str << std::string("tmpprefix=[") + arginfo.getDefaultTmpFileName() + " (prefix for tmp files)" << std::endl;
			str << "sparsetmpprefix=[tmpprefix] (prefix for sparse gap tmp files)" << std::endl;
			// blocksize
			
			se.finish();
			throw se;
		}
		
		
		std::string const inputtype = arginfo.getValue<std::string>("inputtype","bytestream");
		
		if ( inputtype == "compactstream" )
			BwtMergeSort<libmaus2::suffixsort::CompactInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "pac" )
			BwtMergeSort<libmaus2::suffixsort::PacInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "pacterm" )
			BwtMergeSort<libmaus2::suffixsort::PacTermInputTypes>::computeBwt(arginfo);
		else if ( inputtype == "lz4" )
			BwtMergeSort<libmaus2::suffixsort::Lz4InputTypes>::computeBwt(arginfo);
		else if ( inputtype == "utf-8" )
		{
			// compute index of file for random access, if it does not already exist
			std::string const fn = arginfo.getRestArg<std::string>(0);
			std::string const idxfn = fn + ".idx";
			if ( ! ::libmaus2::util::GetFileSize::fileExists(idxfn) )
			{
				::libmaus2::util::Utf8BlockIndex::unique_ptr_type index(::libmaus2::util::Utf8BlockIndex::constructFromUtf8File(fn));
				::libmaus2::aio::OutputStreamInstance COS(idxfn);
				index->serialise(COS);
				COS.flush();
			}
			BwtMergeSort<libmaus2::suffixsort::Utf8InputTypes>::computeBwt(arginfo);
		}
		else if ( inputtype == "bytestream" )
			BwtMergeSort<libmaus2::suffixsort::ByteInputTypes>::computeBwt(arginfo);			
		else
		{
			libmaus2::exception::LibMausException se;
			se.getStream() << "Unknown input type " << inputtype << std::endl;
			se.finish();
			throw se;
		}
			
		#if defined(_OPENMP)
		omp_set_num_threads(maxthreads);		
		#endif
		
		std::cerr << "[M] " << libmaus2::util::MemUsage() << " runtime " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
