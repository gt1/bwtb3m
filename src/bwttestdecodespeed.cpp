/**
    bwtb3m
    Copyright (C) 2009-2015 German Tischler

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
#include <libmaus2/huffman/RLDecoder.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/fm/SampledISA.hpp>
#include <libmaus2/lf/LF.hpp>
#include <libmaus2/lf/MultiRankCacheLF.hpp>
#include <libmaus2/random/Random.hpp>
#include <iostream>

int bwttestdecodespeed(::libmaus2::util::ArgParser const & arg)
{
	std::string const infn = arg[0];
	uint64_t const n = libmaus2::huffman::RLDecoder::getLength(infn);
	libmaus2::huffman::RLDecoder dec(std::vector<std::string>(1,infn));
	std::string const inputtype = arg.uniqueArgPresent("T") ? arg["T"] : "bytestream";
	libmaus2::autoarray::AutoArray<uint8_t> BB(n);
	std::string const isaname = libmaus2::util::OutputFileNameTools::clipOff(infn,".bwt") + ".isa";

	if ( inputtype == "bytestream" )
	{
		// read inverse sampled suffix array
		libmaus2::aio::InputStreamInstance::unique_ptr_type isaISI(new libmaus2::aio::InputStreamInstance(isaname));
		libmaus2::fm::SampledISA<libmaus2::lf::LF>::readUnsignedInt(*isaISI); // isa rate
		libmaus2::autoarray::AutoArray<uint64_t> Aisa = libmaus2::fm::SampledISA<libmaus2::lf::LF>::readArray64(*isaISI);
		isaISI.reset();

		// decode BWT from run length encoding to byte array
		std::pair<int64_t,uint64_t> P;
		uint64_t z = 0;
		while ( (P = dec.decodeRun()).first >= 0 )
		{
			if ( P.first >= 256 )
			{
				::libmaus2::exception::LibMausException se;
				se.getStream() << "unsupported character value " << P.first << " for inputtype=bytestream" << std::endl;
				se.finish();
				throw se;
			}
			for ( uint64_t i = 0; i < P.second; ++i )
				BB[z++] = P.first;
		}

		// compute maximum symbol
		int64_t maxsym = -1;
		for ( uint64_t i = 0; i < BB.size(); ++i )
			maxsym = std::max(maxsym,static_cast<int64_t>(BB[i]));
		// compute alphabet size
		uint64_t const alsize = static_cast<uint64_t>(maxsym + 1);

		for ( uint64_t tpar = 1 ; tpar <= 8; ++tpar )
		{
			// clone BWT
			libmaus2::autoarray::AutoArray<uint8_t> B = BB.clone();
			// compute rank dictionary (one bit vector per symbol)
			libmaus2::lf::MultiRankCacheLF LF(B.begin(),B.size(),alsize);

			// generate random data and read it back to remove B and LF from cache
			{
				libmaus2::autoarray::AutoArray<uint8_t> RANDin(128*1024*1024);
				for ( uint64_t i = 0; i < RANDin.size(); ++i )
					RANDin[i] = libmaus2::random::Random::rand8();
				libmaus2::autoarray::AutoArray<uint8_t> RANDout = RANDin.clone();
			}

			uint64_t const step = (Aisa.size() + tpar - 1)/tpar;
			uint64_t const par = (Aisa.size() + step - 1)/step;
			uint64_t const tsteps = std::min((n + par - 1)/par,static_cast<uint64_t>(128ull*1024ull*1024ull));

			libmaus2::autoarray::AutoArray<uint64_t> R(par);
			for ( uint64_t i = 0; i < par; ++i )
				R[i] = Aisa[i * step];

			libmaus2::timing::RealTimeClock rtc; rtc.start();
			for ( uint64_t i = 0; i < tsteps; ++i )
				for ( unsigned int j = 0; j < par; ++j )
					R[j] = LF.step(B[R[j]],R[j]);

			double const t = rtc.getElapsedSeconds();
			std::cout << tpar << "\t" << t/par << " " << (par * tsteps)/t << std::endl;
		}

		return EXIT_SUCCESS;
	}
	else
	{
		::libmaus2::exception::LibMausException se;
		se.getStream() << "input type " << inputtype << " is currently not supported" << std::endl;
		se.finish();
		throw se;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.size() < 1 )
		{
			std::cerr << "usage: " << arg.progname << " <in.bwt>" << std::endl;
			return EXIT_FAILURE;
		}

		return bwttestdecodespeed(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
