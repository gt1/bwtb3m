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
#include <libmaus2/fm/BidirectionalDnaIndexTemplate.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/StreamFastQReader.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>

template<typename pattern_type>
std::string getQuality(pattern_type const & pattern)
{
	return pattern.quality;
}

template<>
std::string getQuality(libmaus2::fastx::Pattern const & pattern)
{
	return std::string(pattern.spattern.size(),'H');
}

template<typename _reader_type>
int allHamDna(libmaus2::util::ArgInfo const & arginfo)
{
	typedef _reader_type reader_type;

	libmaus2::aio::CheckedInputStream CIS(arginfo.restargs[1]);

	std::istringstream maxmisistr(arginfo.restargs[2]);
	uint64_t maxmis;
	maxmisistr >> maxmis;
	if ( ! maxmisistr )
	{
		std::cerr << "[E] cannot parse " << arginfo.restargs[2] << " as maximum number of mismatches" << std::endl;
		return EXIT_FAILURE;
	}

	std::cerr << "[V] loading index...";
	libmaus2::fm::BidirectionalDnaIndexImpCompactHuffmanWaveletTree index(arginfo.restargs[0]);
	std::cerr << "done." << std::endl;

	std::string const & basename = index.basename;
	std::string const fainame = basename + ".fa.fai";

	if ( ! libmaus2::util::GetFileSize::fileExists(fainame) )
	{
		std::cerr << "[E] file " << fainame << " does not exist." << std::endl;
		return EXIT_FAILURE;
	}

	std::vector<libmaus2::bambam::Chromosome> chromosomes;
	{
		libmaus2::aio::CheckedInputStream CIS(fainame);

		while ( CIS )
		{
			std::string line;
			std::getline(CIS,line);

			if ( line.size() )
			{
				std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(line,std::string("\t"));
				if ( tokens.size() == 5 )
				{
					std::string const & name = tokens[0];
					std::string const & slen = tokens[1];
					std::istringstream slenistr(slen);
					uint64_t len = 0;
					slenistr >> len;
					if ( ! slenistr )
					{
						std::cerr << "[E] failed to parse line " << line << " in " << fainame << std::endl;
						return EXIT_FAILURE;
					}

					chromosomes.push_back(
						libmaus2::bambam::Chromosome(name,len)
					);
				}
			}
		}
	}

	reader_type FAin(CIS);
	typename reader_type::pattern_type pat;
	std::vector< std::pair<uint64_t, libmaus2::fm::BidirectionalIndexInterval > > VBI;
	std::vector<uint64_t> const seqstart = index.getSeqStartPositions();

	if ( seqstart.size() % 2 )
	{
		std::cerr << "[E] Number of terminator symbols in index is not even, index is broken." << std::endl;
		return EXIT_FAILURE;
	}

	uint64_t const numseq = seqstart.size()/2;

	if ( numseq != chromosomes.size() )
	{
		std::cerr << "[E] error: index is inconsistent with fai file" << std::endl;
		return EXIT_FAILURE;
	}

	std::vector<uint64_t> seqlen(numseq);
	for ( uint64_t i = 0; i < seqlen.size(); ++i )
	{
		seqlen[i] = seqstart[2*i+1]-(seqstart[2*i]+1);

		if ( seqlen[i] != chromosomes[i].getLength() )
		{
			std::cerr << "[E] error: index is inconsistent with fai file, length of sequence " << chromosomes[i].getNameString() << " is wrong" << std::endl;
			return EXIT_FAILURE;
		}
	}

	libmaus2::bambam::BamHeader header;
	for ( uint64_t i = 0; i < chromosomes.size(); ++i )
		header.addChromosome(
			chromosomes[i].getNameString(),chromosomes[i].getLength()
		);
	header.initSetup();

	// add PG line to header
	header.text = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		header.text,
		"allHamDna", // ID
		"allHamDna", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(header.text).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN
	);

	libmaus2::bambam::BamWriter wr(std::cout,header);
	::libmaus2::bambam::MdStringComputationContext mdcontext;

	while ( FAin.getNextPatternUnlocked(pat) )
	{
		VBI.clear();
		index.hammingSearchRecUnmapped(
			pat.spattern,
			maxmis,
			std::numeric_limits<uint64_t>::max(),
			VBI
		);

		std::string const & query = pat.spattern;
		std::string const quality = getQuality(pat);
		std::ostringstream cigarstr;
		cigarstr << query.size() << "M";
		std::string const cigar = cigarstr.str();
		bool second = false;

		for ( uint64_t i = 0; i < VBI.size(); ++i )
		{
			uint64_t const mismatches = VBI[i].first;
			libmaus2::fm::BidirectionalIndexInterval const & interval = VBI[i].second;

			for ( uint64_t j = 0; j < interval.siz; ++j )
			{
				uint64_t const rank_f = interval.spf+j;
				uint64_t const rank_r = interval.spr+j;

				uint64_t const fullpos_f = (*(index.SA))[rank_f];
				uint64_t const fullpos_r = (*(index.SA))[rank_r];

				std::vector<uint64_t>::const_iterator ita_f = std::lower_bound(seqstart.begin(),seqstart.end(),fullpos_f);
				if ( ita_f == seqstart.end() || *ita_f != fullpos_f ) --ita_f;
				std::vector<uint64_t>::const_iterator ita_r = std::lower_bound(seqstart.begin(),seqstart.end(),fullpos_r);
				if ( ita_r == seqstart.end() || *ita_r != fullpos_r ) --ita_r;

				uint64_t const prerefid_f = ita_f - seqstart.begin();
				uint64_t const prerefid_r = ita_r - seqstart.begin();

				if ( (prerefid_f & 1) == 0 )
				{
					uint64_t const refid_f = prerefid_f/2;
					uint64_t const pos_f = fullpos_f - seqstart[prerefid_f];

					wr.encodeAlignment(
						pat.sid,refid_f,pos_f,0,
						second ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY : 0, // flags
						cigar,-1,-1,0,query,quality
					);
					wr.putAuxNumber("NM",'i',mismatches);

					std::string const text = index.getTextUnmapped(rank_f,pat.spattern.size());
					std::pair<uint8_t const *,uint64_t> algn = wr.getAlignment();
					::libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(
						algn.first,algn.second,mdcontext,text.begin(),false);

					wr.putAuxString("MD", mdcontext.md.begin());

					wr.commit();

					#if 0
					bool const strand = true;
					std::cerr << "refid=" << refid_f << " pos=" << pos_f << " strand=" << strand << std::endl;
					std::cerr << query << std::endl;
					std::cerr << text << std::endl;

					for ( uint64_t k = 0; k < text.size(); ++k )
						if ( query[k] != text[k] )
							std::cerr.put('-');
						else
							std::cerr.put(' ');
					std::cerr.put('\n');
					#endif

					second = true;
				}
				if ( (prerefid_r & 1) == 0 )
				{
					uint64_t const refid_r = prerefid_r/2;
					uint64_t const pos_r = fullpos_r - seqstart[prerefid_r];
					std::string rquality = quality;
					std::reverse(rquality.begin(),rquality.end());

					wr.encodeAlignment(
						pat.sid,refid_r,pos_r,0,
						libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE |
						(second ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY : 0), // flags
						cigar,-1,-1,0,libmaus2::fastx::reverseComplementUnmapped(query),rquality
					);
					wr.putAuxNumber("NM",'i',mismatches);

					std::string const pretext = index.getTextUnmapped(rank_r,pat.spattern.size());
					std::string const text = libmaus2::fastx::reverseComplementUnmapped(pretext);
					std::pair<uint8_t const *,uint64_t> algn = wr.getAlignment();
					::libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(
						algn.first,algn.second,mdcontext,pretext.begin(),false);

					wr.putAuxString("MD", mdcontext.md.begin());

					wr.commit();

					#if 0
					bool const strand = false;
					std::cerr << "refid=" << refid_r << " pos=" << pos_r << " strand=" << strand << std::endl;
					std::cerr << pat.spattern << std::endl;
					std::cerr << text << std::endl;


					for ( uint64_t k = 0; k < text.size(); ++k )
						if ( query[k] != text[k] )
							std::cerr.put('-');
						else
							std::cerr.put(' ');
					std::cerr.put('\n');
					#endif

					second = true;
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc, argv);

		if ( arginfo.restargs.size() < 3 )
		{
			std::cerr << "[E] usage: " << argv[0] << " <index.hwt> <patterns.fa> <maxmis>" << std::endl;
			return EXIT_FAILURE;
		}

		libmaus2::aio::CheckedInputStream CIS(arginfo.restargs[1]);
		bool fastq;
		if ( CIS.peek() >= 0 && CIS.peek() == '>' )
			fastq = false;
		else if ( CIS.peek() >= 0 && CIS.peek() == '@' )
			fastq = true;
		else
		{
			std::cerr << "[E] cannot detect format of pattern file." << std::endl;
			return EXIT_FAILURE;
		}
		CIS.close();

		if ( fastq )
			return allHamDna<libmaus2::fastx::StreamFastQReaderWrapper>(arginfo);
		else
			return allHamDna<libmaus2::fastx::StreamFastAReaderWrapper>(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
