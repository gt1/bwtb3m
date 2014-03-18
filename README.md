bwtb3m
======

Burrows Wheeler Transform By Balanced Block Merging
---------------------------------------------------

bwtb3m is a set of tools for indexing data.

Source
------

The bwtb3m source code is hosted on github:

	git@github.com:gt1/bwtb3m.git

Compilation of bwtb3m
---------------------

bwtb3m needs libmaus [https://github.com/gt1/libmaus] . When libmaus
is installed in ${LIBMAUSPREFIX} then bwtb3m can be compiled and
installed in ${HOME}/bwtb3m using

	- autoreconf -i -f
	- ./configure --with-libmaus=${LIBMAUSPREFIX} \
		--prefix=${HOME}/bwtb3m
	- make install

Generating an index for BWA
---------------------------

The following script can be used to generate an index for BWA. It expects 4 arguments:

 - the path to the bwa program
 - the path to the bwtb3m program
 - the path to the bwtb3mtobwa program
 - the file name of the input FastA file to be indexed

	#! /bin/bash
	if [ $# -lt 4 ] ; then
		echo "usage: ${SHELL} $0 /path/to/bwa /path/to/bwtb3m /path/to/bwtb3mtobwa <in.fa>"
		exit 1
	fi

	BWA="$1"
	BWTB3M="$2"
	BWTB3MTOBWA="$3"
	INPUT="$4"

	if [ ! -e "${BWA}" ] ; then
		echo "File ${BWA} does not exist"
		exit 1
	fi
	if [ ! -e "${BWTB3M}" ] ; then
		echo "File ${BWTB3M} does not exist"
		exit 1
	fi
	if [ ! -e "${BWTB3MTOBWA}" ] ; then
		echo "File ${BWTB3MTOBWA} does not exist"
		exit 1
	fi

	"${BWA}" fa2pac "${INPUT}"
	"${BWTB3M}" inputtype=pacterm outputfilename=${INPUT}.pac.bwt ${INPUT}.pac
	"${BWTB3MTOBWA}" ${INPUT}.pac.bwt ${INPUT}.bwt ${INPUT}.sa
	"${BWA}" bwtupdate "${INPUT}.bwt"
	rm -f "${INPUT}".pac.*
	"${BWA}" fa2pac -f "${INPUT}"
