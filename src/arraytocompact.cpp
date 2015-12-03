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
#include <iostream>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/serialize/Serialize.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/aio/SynchronousGenericInput.hpp>
#include <libmaus2/math/numbits.hpp>
#include <libmaus2/bitio/CompactArrayWriter.hpp>
#include <libmaus2/bitio/CompactArray.hpp>

/**
 * convert a sampled suffix array or sampled inverse suffix array to a CompactArray based data structure
 **/
int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const fn = arg[0];
		std::string const outfn = fn + ".compact";

		uint64_t samplingrate;
		libmaus2::aio::InputStreamInstance::unique_ptr_type PISI(new libmaus2::aio::InputStreamInstance(fn));
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*PISI,&samplingrate);
		uint64_t n;
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*PISI,&n);

		uint64_t maxval = 0;
		libmaus2::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type Sin(new libmaus2::aio::SynchronousGenericInput<uint64_t>(*PISI,8*1024));
		for ( uint64_t i = 0; i < n; ++i )
		{
			uint64_t v;
			bool const ok = Sin->getNext(v);
			if ( ! ok )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unexpected EOF while reading " << fn << std::endl;
				lme.finish();
				throw lme;
			}
			maxval = std::max(maxval,v);

			if ( (i % (32*1024*1024)) == 0 )
				std::cerr << "[V] " << static_cast<double>(i)/n << std::endl;
		}
		unsigned int const b = libmaus2::math::numbits(maxval);

		Sin.reset();
		PISI.reset();

		std::cerr << "[V] writing array of length " << n << " with " << b << " bits per number" << std::endl;

		libmaus2::aio::InputStreamInstance::unique_ptr_type tPISI(new libmaus2::aio::InputStreamInstance(fn));
		PISI = UNIQUE_PTR_MOVE(tPISI);
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*PISI,&samplingrate);
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*PISI,&n);
		libmaus2::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type tSin(new libmaus2::aio::SynchronousGenericInput<uint64_t>(*PISI,8*1024));
		Sin = UNIQUE_PTR_MOVE(tSin);

		libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(outfn));
		libmaus2::util::NumberSerialisation::serialiseNumber(*OSI,samplingrate);
		libmaus2::bitio::CompactArrayWriter::unique_ptr_type CAW(new libmaus2::bitio::CompactArrayWriter(*OSI,n,b));
		for ( uint64_t i = 0; i < n; ++i )
		{
			uint64_t v;
			bool const ok = Sin->getNext(v);
			if ( ! ok )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unexpected EOF while reading " << fn << std::endl;
				lme.finish();
				throw lme;
			}
			CAW->put(v);
		}
		CAW.reset();
		OSI->flush();
		OSI.reset();

		std::cerr << "[V] checking...";
		libmaus2::aio::InputStreamInstance::unique_ptr_type RISI(new libmaus2::aio::InputStreamInstance(outfn));
		uint64_t const rsamplingrate = libmaus2::util::NumberSerialisation::deserialiseNumber(*RISI);
		assert ( samplingrate == rsamplingrate );
		libmaus2::bitio::CompactArray CA(*RISI);

		libmaus2::aio::InputStreamInstance::unique_ptr_type DISI(new libmaus2::aio::InputStreamInstance(fn));
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*DISI,&samplingrate);
		libmaus2::serialize::Serialize<uint64_t>::deserialize(*DISI,&n);
		assert ( n == CA.size() );

		libmaus2::aio::SynchronousGenericInput<uint64_t>::unique_ptr_type Srin(new libmaus2::aio::SynchronousGenericInput<uint64_t>(*DISI,8*1024));

		for ( uint64_t i = 0; i < n; ++i )
		{
			uint64_t v;
			bool const ok = Srin->getNext(v);
			if ( ! ok )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unexpected EOF while reading " << fn << std::endl;
				lme.finish();
				throw lme;
			}
			bool const vok = (CA[i] == v);

			if ( ! vok )
			{
				std::cerr << "[D] failure i=" << i << " v=" << v << " CA[i]=" << CA[i] << std::endl;
			}
			assert ( CA[i] == v );
		}
		std::cerr << "done." << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
