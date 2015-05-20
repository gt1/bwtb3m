/*
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
*/
#include <config.h>

#include <libmaus2/aio/CheckedInputStream.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/StreamFastQReader.hpp>
#include <libmaus2/lf/ImpCompactHuffmanWaveletLF.hpp>
#include <libmaus2/lz/BufferedGzipStream.hpp>
#include <libmaus2/suffixtree/CompressedSuffixTree.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/util/ToUpperTable.hpp>

void evaluateAcc(std::ostream & out, std::vector<uint64_t> const & acccols, std::vector<uint64_t> vacc)
{
	std::vector < uint64_t > acccnt(acccols.size());
	std::sort(vacc.begin(),vacc.end());
	
	uint64_t l = 0;
	while ( l != vacc.size() )
	{
		uint64_t h = l;
		while ( h != vacc.size() && vacc[h] == vacc[l] )
			++h;
			
		for ( uint64_t i = 0; i < acccols.size(); ++i )
			if ( vacc[l] >= acccols[i] )
				acccnt[i] += (h-l);
			
		l = h;
	}
	
	for ( uint64_t i = 0; i < acccnt.size(); ++i )
		out << "\t" << acccnt[i];
}

/*
 * wavelet tree based method based on backward search
 */
template<typename _lf_type>
int probeScanDnaHwt(
	libmaus2::util::ArgInfo const & arginfo, unsigned int const probelen, std::string const suffix,
	std::vector<uint64_t> const & acccols
)
{
	// load index for reference
	typedef _lf_type lf_type;
	std::string const prefix = arginfo.getRestArg<std::string>(0);
	std::string const hwtname = prefix+suffix;
	std::cerr << "[V] loading index...";
	lf_type LF(hwtname);
	std::cerr << "done." << std::endl;

	// open queries
	libmaus2::fastx::StreamFastAReaderWrapper queriesIn(std::cin);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

	// read next query
	while ( queriesIn.getNextPatternUnlocked(pattern) )
	{
		// skip if query is too short
		if ( pattern.spattern.size() < probelen )
		{
			std::cout << pattern.sid << "\t" << pattern.spattern << "\n";
		}
		else
		{
			std::cout << pattern.sid << "\t" << pattern.spattern << "\t";
			
			std::vector<uint64_t> vacc;
			
			// look for probelen-mers
			for ( uint64_t i = 0; i < pattern.spattern.size()-probelen+1; ++i )
			{
				// forward
				std::string fquery = pattern.spattern.substr(i,probelen);
				// reverse complement
				std::string rquery = libmaus2::fastx::reverseComplementUnmapped(fquery);

				// suffix array intervals
				std::pair<uint64_t,uint64_t> fnode(0ull,LF.getN());
				std::pair<uint64_t,uint64_t> rnode(0ull,LF.getN());
				// backward search
				for ( uint64_t j = 0; j < probelen; ++j )
				{
					fnode = LF.step(libmaus2::fastx::mapChar(fquery[probelen-j-1])+1,fnode.first,fnode.second);
					rnode = LF.step(libmaus2::fastx::mapChar(rquery[probelen-j-1])+1,rnode.first,rnode.second);
				}
				
				uint64_t lcnt;
		
				if ( fquery == rquery )
				{
					lcnt = fnode.second-fnode.first;
					std::cout << lcnt << "(" << fnode.second-fnode.first << ")";
				}
				else
				{
					lcnt = (fnode.second-fnode.first)+(rnode.second-rnode.first);
					std::cout << lcnt << "(" << fnode.second-fnode.first << "," << rnode.second-rnode.first << ")";				
				}
				
				if ( i + 1 < pattern.spattern.size()-probelen+1 )
					std::cout.put(';');
					
				vacc.push_back(lcnt);
			}

			evaluateAcc(std::cout,acccols,vacc);
			
			std::cout << "\n";		
		}
	}	
	
	return EXIT_SUCCESS;
}

/*
 * kmer word based method
 */
template<typename uint_type>
int probeScanDna(
	libmaus2::util::ArgInfo const & arginfo, unsigned int const probelen,
	std::vector<uint64_t> const & acccols
)
{
	// this should be non empty
	assert ( probelen );
	// bits per k-mer
	unsigned int const probesymbits = 3 * probelen;
	// lookup table key bits
	unsigned int const lookupbits = std::min(static_cast<unsigned int>(arginfo.getValue<unsigned int>("lookupbits",24)),probesymbits);
	// rest bits
	unsigned int const probeshiftbits = probesymbits - lookupbits;

	// bit mask for k-mers
	uint_type mask = 0;
	
	for ( uint64_t i = 0; i < probelen; ++i )
	{
		mask <<= 3;
		mask |= 7;
	}

	std::cerr << "[V] reading queries...";
	libmaus2::fastx::StreamFastAReaderWrapper queriesIn(std::cin);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
	std::vector < std::pair<std::string,std::string> > queries;

	std::vector<uint_type> kmers;
	while ( queriesIn.getNextPatternUnlocked(pattern) )
	{
		queries.push_back(std::pair<std::string,std::string>(pattern.sid,pattern.spattern));
		
		std::string & pat = pattern.spattern;

		// unify symbols
		for ( uint64_t i = 0; i < pat.size(); ++i )
			pat[i] = libmaus2::fastx::remapChar(libmaus2::fastx::mapChar(pat[i]));
		
		if ( pat.size() < probelen )
			continue;

		// forward k-mers
		uint_type m = 0;
		for ( std::string::size_type i = 0; i < probelen-1; ++i )
		{
			m <<= 3;
			m |= libmaus2::fastx::mapChar(pat[i]);
		}

		for ( std::string::size_type i = probelen-1; i < pat.size(); ++i )
		{
			m <<= 3;
			m &= mask;
			m |= libmaus2::fastx::mapChar(pat[i]);

			kmers.push_back(m);
		}
		
		pat = libmaus2::fastx::reverseComplementUnmapped(pat);

		// reverse complement k-mers
		uint_type r = 0;
		for ( std::string::size_type i = 0; i < probelen-1; ++i )
		{
			r <<= 3;
			r |= libmaus2::fastx::mapChar(pat[i]);
		}

		for ( std::string::size_type i = probelen-1; i < pat.size(); ++i )
		{
			r <<= 3;
			r &= mask;
			r |= libmaus2::fastx::mapChar(pat[i]);

			kmers.push_back(r);
		}

		pat = libmaus2::fastx::reverseComplementUnmapped(pat);
	}
	std::cerr << "done." << std::endl;
	
	std::cerr << "[V] sorting kmers...";
	std::sort(kmers.begin(),kmers.end());
	std::cerr << "done." << std::endl;

	std::cerr << "[V] reducing to unique kmers...";
	kmers.resize ( std::unique (kmers.begin(), kmers.end()) - kmers.begin() );
	std::cerr << "done, number is " << kmers.size() << std::endl;

	std::cerr << "[V] producing lookup table for " << lookupbits << " bits...";
	libmaus2::autoarray::AutoArray< std::pair<uint32_t,uint32_t> > Alookup( 1ull << lookupbits );
	uint64_t l = 0;
	while ( l != kmers.size() )
	{
		uint64_t h = l;
		while ( h != kmers.size() && (kmers[h]>>probeshiftbits)==(kmers[l]>>probeshiftbits) )
			++h;
			
		assert ( (kmers[l]>>probeshiftbits) < Alookup.size() );
		
		Alookup[ kmers[l]>>probeshiftbits ] = std::pair<uint32_t,uint32_t>(l,h);
			
		l = h;
	}
	std::cerr << "done." << std::endl;

	// open reference file
	libmaus2::aio::CheckedInputStream gzCIS(arginfo.restargs.at(0));
	libmaus2::lz::BufferedGzipStream gzin(gzCIS);
	libmaus2::fastx::StreamFastAReaderWrapper refin(gzin);
	libmaus2::util::ToUpperTable toup;

	// kmer counts
	std::vector<uint64_t> cnts(kmers.size());
	while ( refin.getNextPatternUnlocked(pattern) )
	{
		std::cerr << "[V] matching " << pattern.sid << "...";
		
		// convert reference sequence to upper case
		std::string & text = pattern.spattern;
		for ( uint64_t i = 0; i < text.size(); ++i )
			text[i] = text[i] = static_cast<char>((toup(static_cast<unsigned char>(text[i]))));	
		
		// do nothing if reference sequence is shorter than k-mer length
		if ( text.size() < probelen )
			continue;

		uint_type v = 0;
		for ( uint64_t i = 0; i < probelen-1; ++i )
		{
			v <<= 3;
			v |= libmaus2::fastx::mapChar(text[i]);
		}
		for ( uint64_t i = probelen-1; i < text.size() ; ++i )
		{
			v <<= 3;
			v &= mask;		
			v |= libmaus2::fastx::mapChar(text[i]);
			
			uint64_t const il =  Alookup[ v >> probeshiftbits ].first;
			uint64_t const ih = Alookup[ v >> probeshiftbits ].second;
			
			if ( ih != il )
			{
				typename std::vector<uint_type>::const_iterator low = kmers.begin() + il;
				typename std::vector<uint_type>::const_iterator high = kmers.begin() + ih;
				
				std::pair<
					typename std::vector<uint_type>::const_iterator, 
					typename std::vector<uint_type>::const_iterator > eq = std::equal_range(low,high,v);
				
				if ( eq.second != eq.first )
				{
					assert ( *(eq.first) == v );
					cnts[eq.first-kmers.begin()]++;
				}
			}
			
			if ( ((i+1) & (1024*1024-1)) == 0 )
				std::cerr << "(" << (i+1)/static_cast<double>(text.size()) << ")";
		}

		std::cerr << "done." << std::endl;
	}

	// print k-mer counts for each query
	for ( uint64_t q = 0; q < queries.size(); ++q )
	{
		std::string const & id = queries[q].first;
		std::string & pat = queries[q].second;
		
		if ( pat.size() < probelen )
		{
			std::cout << id << "\t" << pat << "\n";
		}
		else
		{
			std::vector<uint_type> lkmersf;
			std::vector<uint_type> lkmersr;
			
			uint_type m = 0;
			for ( std::string::size_type i = 0; i < probelen-1; ++i )
			{
				m <<= 3;
				m |= libmaus2::fastx::mapChar(pat[i]);
			}

			for ( std::string::size_type i = probelen-1; i < pat.size(); ++i )
			{
				m <<= 3;
				m &= mask;
				m |= libmaus2::fastx::mapChar(pat[i]);

				lkmersf.push_back(m);
			}
			
			pat = libmaus2::fastx::reverseComplementUnmapped(pat);

			uint_type r = 0;
			for ( std::string::size_type i = 0; i < probelen-1; ++i )
			{
				r <<= 3;
				r |= libmaus2::fastx::mapChar(pat[i]);
			}

			for ( std::string::size_type i = probelen-1; i < pat.size(); ++i )
			{
				r <<= 3;
				r &= mask;
				r |= libmaus2::fastx::mapChar(pat[i]);

				lkmersr.push_back(r);
			}

			pat = libmaus2::fastx::reverseComplementUnmapped(pat);
		
			std::reverse(lkmersr.begin(),lkmersr.end());

			std::ostringstream ostr;

			std::vector<uint64_t> vacc;
			
			ostr << id << "\t" << pat << "\t";
			
			for ( uint64_t i = 0; i < lkmersf.size(); ++i )
			{
				std::pair<
					typename std::vector<uint_type>::const_iterator, 
					typename std::vector<uint_type>::const_iterator > eqf = std::equal_range(kmers.begin(),kmers.end(),lkmersf[i]);
				
				uint64_t fcnt = 0;
				if ( eqf.first != eqf.second )
				{
					assert ( eqf.second-eqf.first == 1 );
					fcnt = cnts[eqf.first - kmers.begin()];
				}
				
				std::pair<
					typename std::vector<uint_type>::const_iterator, 
					typename std::vector<uint_type>::const_iterator > eqr = std::equal_range(kmers.begin(),kmers.end(),lkmersr[i]);	
					
				uint64_t rcnt = 0;
				if ( eqr.first != eqr.second )
				{
					assert ( eqr.second-eqr.first == 1 );					
					rcnt = cnts[eqr.first - kmers.begin()];
				}
				
				uint64_t lcnt = (lkmersf[i] != lkmersr[i]) ? (fcnt + rcnt) : fcnt;
				
				ostr << lcnt;
				
				ostr << "(";
				
				if ( lkmersf[i] != lkmersr[i] )
				{
					ostr << fcnt << "," << rcnt;
				}
				else
				{
					ostr << fcnt;
				}
				
				ostr << ")";
				
				if ( i+1 < lkmersf.size() )
					ostr << ";";

				vacc.push_back(lcnt);
			}


			evaluateAcc(ostr,acccols,vacc);
			
			std::cout << ostr.str() << "\n";
		}
	}
	
	std::cout.flush();
	                                                                   		
	return EXIT_SUCCESS;
}

static std::vector<uint64_t> parseAccCols(std::string const & scols)
{
	std::deque<std::string> const stokens = libmaus2::util::stringFunctions::tokenize(scols,std::string(","));
	std::vector<uint64_t> acccols;
	
	for ( uint64_t i = 0; i < stokens.size(); ++i )
	{
		std::string const stoken = stokens[i];

		uint64_t v = 0;

		for ( uint64_t j = 0; j < stoken.size(); ++j )
		{
			if ( !isdigit(static_cast<unsigned char>(stoken[j])) )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Cannot parse " << stoken << " as a number." << std::endl;
				lme.finish();
				throw lme;
			}
			else
			{
				v *= 10;
				v += static_cast<unsigned char>(stoken[j]) - '0';
			}
		}
		
		acccols.push_back(v);
	}
	
	return acccols;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc, argv);
		
		if ( arginfo.restargs.size() < 1 )
		{
			std::cerr << "[E] usage: " << argv[0] << " <ref_fa.gz>" << std::endl;
			return EXIT_FAILURE;
		}

		uint64_t const probelen = arginfo.getValue<unsigned int>("probelen",20);
		std::string const mode = arginfo.getValue<std::string>("mode","words");
		std::vector<uint64_t> acccols = parseAccCols(arginfo.getUnparsedValue("acccols",std::string("2,4,10,20,50,100")));
		
		std::cerr << "[V] accumulation columns ";
		for ( uint64_t i = 0; i < acccols.size(); ++i )
			std::cerr << acccols[i] << ((i+1<acccols.size())?",":"");
		std::cerr << std::endl;
	
		if ( mode == "words" )
		{
			if ( 3*probelen <= 64 )
				return probeScanDna<uint64_t>(arginfo,probelen,acccols);
			#if defined(LIBMAUS_HAVE_UNSIGNED_INT128)
			else if ( 3*probelen <= 128 )
				return probeScanDna<libmaus2::uint128_t>(arginfo,probelen,acccols);
			#endif
			else
			{
				std::cerr << "[E] Cannot handle length of probe." << std::endl;
				return EXIT_FAILURE;
			}
		}
		else if ( mode == "rlhwt" )
		{
			return probeScanDnaHwt<libmaus2::lf::ImpCompactRLHuffmanWaveletLF>(arginfo,probelen,".rlhwt",acccols);
		}
		else if ( mode == "hwt" )
		{
			return probeScanDnaHwt<libmaus2::lf::ImpCompactHuffmanWaveletLF>(arginfo,probelen,".hwt",acccols);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
