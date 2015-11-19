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
#include <libmaus2/aio/BufferedOutput.hpp>
#include <libmaus2/aio/InputStreamFactoryContainer.hpp>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/sorting/MergingReadBack.hpp>
#include <libmaus2/util/Histogram.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		std::string const isain = arginfo.getUnparsedRestArg(0);
		std::string const isaout = arginfo.getUnparsedRestArg(1);

		std::string const tmpdir = arginfo.getUnparsedValue("tmpdir",arginfo.getCurDir());
		std::string const tmpout = tmpdir + "/" + arginfo.getDefaultTmpFileName() + "_sort_isa.tmp";

		libmaus2::util::TempFileRemovalContainer::addTempFile(tmpout);

		// std::string const tmpout = isaout + ".tmp";
		uint64_t const sortbufsize = arginfo.getValueUnsignedNumeric<uint64_t>("sortbufsize",16*1024*1024);
		uint64_t const inbufsize = arginfo.getValueUnsignedNumeric<uint64_t>("inbufsize",8*1024);
		bool const verbose = arginfo.getValue<int>("verbose",1);

		libmaus2::aio::OutputStream::unique_ptr_type Ptmp(libmaus2::aio::OutputStreamFactoryContainer::constructUnique(tmpout));
		libmaus2::aio::SortingBufferedOutput< std::pair<uint64_t,uint64_t> >::unique_ptr_type Psortout(
			new libmaus2::aio::SortingBufferedOutput< std::pair<uint64_t,uint64_t> >(*Ptmp,sortbufsize));

		libmaus2::aio::InputStream::unique_ptr_type Pin(libmaus2::aio::InputStreamFactoryContainer::constructUnique(isain));
		libmaus2::aio::SynchronousGenericInput<uint64_t> Sin(*Pin,inbufsize);
		uint64_t v0 = 0, v1 = 0;
		uint64_t c_in = 0;

		while ( Sin.peekNext(v0) )
		{
			bool const ok0 = Sin.getNext(v0);
			assert ( ok0 );
			bool const ok1 = Sin.getNext(v1);
			if ( ! ok1 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Uneven number of words in input file." << std::endl;
				lme.finish();
				throw lme;
			}
			Psortout->put(std::pair<uint64_t,uint64_t>(v1,v0));
			c_in += 1;

			if ( verbose && (c_in & (1024*1024-1)) == 0 )
				std::cerr << "[V] in " << c_in/(1024*1024) << std::endl;
		}
		Psortout->flush();
		std::vector<uint64_t> const blocksizes = Psortout->getBlockSizes();
		Psortout.reset();
		Ptmp->flush();
		Ptmp.reset();

		// libmaus2::aio::InputStream::unique_ptr_type Ptmpin(libmaus2::aio::InputStreamFactoryContainer::constructUnique(tmpout));
		libmaus2::sorting::MergingReadBack< std::pair<uint64_t,uint64_t> >::unique_ptr_type Min(
			new libmaus2::sorting::MergingReadBack< std::pair<uint64_t,uint64_t> >(tmpout,blocksizes)
		);

		std::pair<uint64_t,uint64_t> P;
		libmaus2::aio::OutputStream::unique_ptr_type Pout(libmaus2::aio::OutputStreamFactoryContainer::constructUnique(isaout));
		uint64_t c_out = 0;
		bool gotfirst = false;
		bool gotfirstdif = false;
		bool difconsistent = true;
		uint64_t firstpos = std::numeric_limits<uint64_t>::max();
		uint64_t prevpos = std::numeric_limits<uint64_t>::max();
		uint64_t firstdif = std::numeric_limits<uint64_t>::max();
		libmaus2::util::Histogram hist;
		while ( Min->getNext(P) )
		{
			libmaus2::util::NumberSerialisation::serialiseNumber(*Pout,P.first);
			libmaus2::util::NumberSerialisation::serialiseNumber(*Pout,P.second);

			c_out += 1;

			if ( verbose && (c_out & (1024*1024-1)) == 0 )
				std::cerr << "[V] out " << static_cast<double>(c_out) / c_in << " first pos " << firstpos << " first dif " << firstdif << " consistent " << difconsistent << std::endl;

			if ( ! gotfirst )
			{
				gotfirst = true;
				firstpos = P.first;
			}
			else
			{
				assert ( prevpos != std::numeric_limits<uint64_t>::max() );
				assert ( P.first > prevpos );
				uint64_t const dif = P.first - prevpos;

				hist(dif);

				if ( ! gotfirstdif )
				{
					gotfirstdif = true;
					firstdif = dif;
				}
				else
				{
					if ( dif != firstdif )
						difconsistent = false;
				}
			}
			prevpos = P.first;
		}

		if ( verbose )
			std::cerr << "[V] out " << static_cast<double>(c_out) / c_in << " first pos " << firstpos << " first dif " << firstdif << " consistent " << difconsistent << std::endl;


		Min.reset();
		Pout->flush();
		Pout.reset();

		hist.print(std::cerr);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
