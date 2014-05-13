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
#include <libmaus/huffman/RLDecoder.hpp>
#include <libmaus/util/ArgInfo.hpp>
#include <libmaus/util/utf8.hpp>
#include <libmaus/util/NumberMapSerialisation.hpp>
#include <libmaus/util/OutputFileNameTools.hpp>
#include <libmaus/util/TempFileRemovalContainer.hpp>
#include <libmaus/huffman/HuffmanTree.hpp>
#include <libmaus/rank/RunLengthBitVectorGenerator.hpp>
#include <libmaus/aio/CheckedInputOutputStream.hpp>
#include <libmaus/wavelet/ImpCompactHuffmanWaveletTree.hpp>
#include <iostream>

void hwtToRlHwt(::libmaus::util::ArgInfo const & arginfo)
{
	uint64_t const rlblocksize = arginfo.getValueUnsignedNumeric<uint64_t>("blocksize",64*1024);

	std::string const infn = arginfo.getRestArg<std::string>(0);	
        std::string const rlhwtname = libmaus::util::OutputFileNameTools::clipOff(infn,".bwt") + ".rlhwt";
        
        if ( ! libmaus::util::GetFileSize::fileExists(rlhwtname) )
        {
		std::string const histfn = libmaus::util::OutputFileNameTools::clipOff(infn,".bwt") + ".hist";

		libmaus::aio::CheckedInputStream histstr(histfn);
		std::map<int64_t,uint64_t> const hist = ::libmaus::util::NumberMapSerialisation::deserialiseMap<std::istream,int64_t,uint64_t>(histstr);
		std::vector<uint64_t> syms;
		uint64_t n = 0;

		for ( std::map<int64_t,uint64_t>::const_iterator ita = hist.begin();
			ita != hist.end(); ++ita )
		{
			syms.push_back(ita->first);
			n += ita->second;
		}

		libmaus::huffman::HuffmanTree H(syms);
		libmaus::huffman::HuffmanTree::EncodeTable E(H);
		
		std::vector<uint64_t> nodebitcnts(H.inner());
		
		for ( uint64_t i = 0; i < syms.size(); ++i )
		{
			uint64_t const s = syms[i];
			assert ( hist.find(s) != hist.end() );
			uint64_t const f = hist.find(s)->second;
			unsigned int const codelen = E.getCodeLength(s);
			uint64_t node = H.root();
			
			for ( unsigned int j = 0; j < codelen; ++j )
			{
				uint64_t const nodeid = node - H.leafs();
				nodebitcnts[nodeid] += f;
			
				bool const b = E.getBitFromTop(s,j);
				
				if ( b )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
		}

		std::string const tmpfileprefix = arginfo.getDefaultTmpFileName();
		std::string const rlfnprefix = tmpfileprefix + "_rl";
		libmaus::util::TempFileRemovalContainer::setup();
		std::vector<std::string> rlfilenames(H.inner());
		std::vector<std::string> rlidxfilenames(H.inner());
		libmaus::autoarray::AutoArray<libmaus::aio::CheckedOutputStream::unique_ptr_type> rloutfiles(H.inner());
		libmaus::autoarray::AutoArray<libmaus::aio::CheckedInputOutputStream::unique_ptr_type> rlidxoutfiles(H.inner());
		libmaus::autoarray::AutoArray<libmaus::rank::RunLengthBitVectorGenerator::unique_ptr_type> rlgens(H.inner());
		// libmaus::autoarray::AutoArray<>
		for ( uint64_t i = 0; i < H.inner(); ++i )
		{
			std::ostringstream fnostr;
			fnostr << rlfnprefix << "_" << std::setw(6) << std::setfill('0') << i << std::setw(0);
			rlfilenames[i] = fnostr.str() + ".rl";
			rlidxfilenames[i] = fnostr.str() + ".rlidx";
			
			libmaus::util::TempFileRemovalContainer::addTempFile(rlfilenames[i]);
			libmaus::util::TempFileRemovalContainer::addTempFile(rlidxfilenames[i]);
			
			libmaus::aio::CheckedOutputStream::unique_ptr_type trl(
				new libmaus::aio::CheckedOutputStream(rlfilenames[i])
			);
			rloutfiles[i] = UNIQUE_PTR_MOVE(trl);
			
			libmaus::aio::CheckedInputOutputStream::unique_ptr_type trlidx(
				new libmaus::aio::CheckedInputOutputStream(rlidxfilenames[i])
			);
			rlidxoutfiles[i] = UNIQUE_PTR_MOVE(trlidx);
			
			libmaus::rank::RunLengthBitVectorGenerator::unique_ptr_type tgen(
				new libmaus::rank::RunLengthBitVectorGenerator(
					*rloutfiles[i],
					*rlidxoutfiles[i],
					nodebitcnts[i],
					rlblocksize
				)
			);
			
			rlgens[i] = UNIQUE_PTR_MOVE(tgen);
		}
		
		libmaus::huffman::RLDecoder dec(std::vector<std::string>(1,infn));
		int64_t sym = -1;
		uint64_t fin = 0;
		while ( (sym=dec.get()) >= 0 )
		{
			unsigned int const codelen = E.getCodeLength(sym);
			uint64_t node = H.root();
			
			for ( unsigned int i = 0; i < codelen; ++i )
			{
				uint64_t const nodeid = node - H.leafs();
				bool const b = E.getBitFromTop(sym,i);
				
				rlgens[nodeid]->putbit(b);
				
				if ( b )
					node = H.rightChild(node);
				else
					node = H.leftChild(node);
			}
			
			if ( ++fin % (128*1024*1024) == 0 )
			{
				if ( isatty(STDERR_FILENO) )
					std::cerr << '\r' << std::string(60,' ') << '\r' <<
						static_cast<double>(fin)/n << "\t" << fin << "/" << n << std::flush;
				else
					std::cerr << static_cast<double>(fin)/n << "\t" << fin << "/" << n << '\n';
			}
		}

		if ( isatty(STDERR_FILENO) )
			std::cerr << '\r' << std::string(60,' ') << '\r' <<
				static_cast<double>(fin)/n << "\t" << fin << "/" << n << std::endl;
		else
			std::cerr << static_cast<double>(fin)/n << "\t" << fin << "/" << n << std::endl;


		for ( uint64_t i = 0; i < H.inner(); ++i )
		{
			rlgens[i]->flush();
			rlgens[i].reset();
			rlidxoutfiles[i]->flush();
			rlidxoutfiles[i].reset();
			rloutfiles[i]->flush();
			rloutfiles[i].reset();
			remove(rlidxfilenames[i].c_str());
		}
		
		libmaus::aio::CheckedOutputStream rlhwtCOS(rlhwtname);
		
		uint64_t p = 0;
		p += ::libmaus::util::NumberSerialisation::serialiseNumber(rlhwtCOS,n);
		p += H.serialise(rlhwtCOS);
		p += ::libmaus::util::NumberSerialisation::serialiseNumber(rlhwtCOS,H.inner());
		
		std::vector<uint64_t> nodepos(H.inner());
		for ( uint64_t i = 0; i < H.inner(); ++i )
		{
			nodepos[i] = p;

			uint64_t const fs = libmaus::util::GetFileSize::getFileSize(rlfilenames[i]);
			libmaus::aio::CheckedInputStream CIS(rlfilenames[i]);
			libmaus::util::GetFileSize::copy(CIS,rlhwtCOS,fs);
			
			p += fs;
			
			remove(rlfilenames[i].c_str());
		}

		uint64_t const ip = p;        
		p += ::libmaus::util::NumberSerialisation::serialiseNumberVector(rlhwtCOS,nodepos);
		p += ::libmaus::util::NumberSerialisation::serialiseNumber(rlhwtCOS,ip);
			
		rlhwtCOS.flush();
		rlhwtCOS.close();

		std::cerr << "[V] size of rlhwt is " << p << std::endl;
	}

	bool const verify = arginfo.getValue<unsigned int>("verify",false);
	
	if ( verify )
	{
		libmaus::aio::CheckedInputStream rlhwtCIS(rlhwtname);
		libmaus::wavelet::ImpCompactRLHuffmanWaveletTree ICRLHWT(rlhwtCIS);
		#if defined(_OPENMP)
		uint64_t const numthreads = omp_get_max_threads();
		#else
		uint64_t const numthreads = 1;
		#endif
		
		for ( uint64_t i = 0; i < ICRLHWT.dicts.size(); ++i )
		{
			double const avg = ICRLHWT.dicts[i]->getAvgBlockBitLength();
			std::cerr << "node " << i << " avg block bit length " 
				<< avg 
				<< " (" << avg / 8.0 << " bytes)"
				<< std::endl;

			libmaus::util::Histogram::unique_ptr_type hist(ICRLHWT.dicts[i]->getRunLengthHistogram());
			hist->printFrac(std::cerr);
		}
		
		uint64_t const n = ICRLHWT.size();
		uint64_t const packsize = (n + numthreads - 1)/numthreads;
		uint64_t const numpacks = (n + packsize-1)/packsize;
		libmaus::parallel::SynchronousCounter<uint64_t> fin(0);
		libmaus::timing::RealTimeClock rtc;
		rtc.start();
		
		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif
		for ( int64_t t = 0; t < static_cast<int64_t>(numpacks); ++t )
		{
			uint64_t const low = t * packsize;
			uint64_t const high = std::min(low+packsize,n);

			libmaus::huffman::RLDecoder debdec(std::vector<std::string>(1,infn),low);
			
			for ( uint64_t i = low; i < high; ++i )
			{
				int64_t const sym = ICRLHWT[i];
				int64_t const debsym = debdec.decode();
				assert ( sym == debsym );
				
				uint64_t lfin;
				if ( (lfin=++fin) % (1024*1024) == 0 )
				{
					std::cerr 
						<< '\r' << std::string(80,' ') << '\r' << lfin/(1024*1024) 
						<< " (" << static_cast<double>(lfin)/n << ")"
						<< " " << (lfin/rtc.getElapsedSeconds())/1.0e6
						<< std::flush;
				}
			}
		}
		
		std::cerr 
			<< '\r' << std::string(80,' ') << '\r' << fin.get()
			<< " (" << static_cast<double>(fin.get())/n << ")"
			<< " " << (fin.get()/rtc.getElapsedSeconds())/1.0e6
			<< std::endl;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus::util::ArgInfo const arginfo(argc,argv);
		
		if ( arginfo.helpRequested() || arginfo.restargs.size() < 1 )
		{
			::libmaus::exception::LibMausException se;
			std::ostream & str = se.getStream();
			str << "usage: " << argv[0] << " <in.bwt>" << std::endl;
			str << std::endl;
			
			se.finish();
			throw se;
		}
		
		hwtToRlHwt(arginfo);

		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
