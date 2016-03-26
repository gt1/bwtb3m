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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/wavelet/RlToHwtBase.hpp>
#include <libmaus2/huffman/RLDecoder.hpp>
#include <libmaus2/wavelet/Utf8ToImpCompactHuffmanWaveletTree.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/parallel/NumCpus.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const bwt = arg[0];
		std::string const hwt = ::libmaus2::util::OutputFileNameTools::clipOff(bwt,".bwt") + ".hwt";
		unsigned int const numthreads = libmaus2::parallel::NumCpus::getNumLogicalProcessors();

		::libmaus2::util::Histogram::unique_ptr_type mhist(libmaus2::wavelet::RlToHwtBase<true,libmaus2::huffman::RLDecoder>::computeRlSymHist(bwt,numthreads));
		::std::map<int64_t,uint64_t> const chist = mhist->getByType<int64_t>();
		::libmaus2::huffman::HuffmanTree H ( chist.begin(), chist.size(), false, true, true );
		libmaus2::wavelet::Utf8ToImpCompactHuffmanWaveletTree::constructWaveletTreeFromRl<libmaus2::huffman::RLDecoder,true>(
			bwt,hwt,hwt + "_tmp",H,libmaus2::wavelet::Utf8ToImpCompactHuffmanWaveletTree::getDefaultTPartSizeMax(),numthreads
		);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
