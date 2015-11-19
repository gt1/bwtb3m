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
#include <libmaus2/fm/MausFmToBwaConversion.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <iostream>

void bwtb3mtobwa(::libmaus2::util::ArgInfo const & arginfo)
{
	std::string const infn = arginfo.getRestArg<std::string>(0);
	std::string const outbwt = arginfo.getRestArg<std::string>(1);
	std::string const outsa = arginfo.getRestArg<std::string>(2);

	::libmaus2::fm::MausFmToBwaConversion::rewrite(infn,outbwt,outsa);
}

int main(int argc, char * argv[])
{
	try
	{
		::libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( arginfo.helpRequested() || arginfo.restargs.size() < 3 )
		{
			::libmaus2::exception::LibMausException se;
			std::ostream & str = se.getStream();
			str << "usage: " << argv[0] << " <in.bwt> <out.bwt> <out.sa>" << std::endl;
			str << std::endl;

			se.finish();
			throw se;
		}

		bwtb3mtobwa(arginfo);

		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
