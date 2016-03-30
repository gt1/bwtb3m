/**
    bwtb3m
    Copyright (C) 2009-2016 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

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
#include <libmaus2/fastx/FastAToCompact4BigBandBiDir.hpp>

int fagzToCompact4BigBandBiDir(libmaus2::util::ArgInfo const & arginfo)
{
	int const verbose = arginfo.getValue<int>("verbose",1);
	std::ostream * logstr = 0;
	if ( verbose )
		logstr = &(std::cerr);

	std::vector<std::string> inputfilenames;
	inputfilenames = arginfo.restargs;

	if ( arginfo.hasArg("inputfilenames") )
	{
		std::string const inf = arginfo.getUnparsedValue("inputfilenames",std::string());
		libmaus2::aio::InputStreamInstance Pinf(inf);
		while ( Pinf )
		{
			std::string line;
			std::getline(Pinf,line);
			if ( line.size() )
				inputfilenames.push_back(line);
		}
	}

	std::string const outputfilename = arginfo.getUnparsedValue("outputfilename",std::string());

	return libmaus2::fastx::FastAToCompact4BigBandBiDir::fastaToCompact4BigBandBiDir(inputfilenames,logstr,true /* single strand */,outputfilename);
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		return fagzToCompact4BigBandBiDir(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
