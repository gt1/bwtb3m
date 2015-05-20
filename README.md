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

bwtb3m needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then bwtb3m can be compiled and
installed in ${HOME}/bwtb3m using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/bwtb3m
	- make install

Generating an index for BWA
---------------------------

The following script can be used to generate an index for BWA. It expects 5 arguments:

 - the path to the bwa program
 - the path to the bwtb3m program
 - the path to the bwtb3mtobwa program
 - the amount of memory for bwtb3m in giga bytes
 - the file name of the input FastA file to be indexed

```
#! /bin/bash
if [ $# -lt 5 ] ; then
	echo "usage: ${SHELL} $0 /path/to/bwa /path/to/bwtb3m /path/to/bwtb3mtobwa <mem/GB> <in.fa>"
	exit 1
fi

BWA="$1"
BWTB3M="$2"
BWTB3MTOBWA="$3"
MEM="$4"
INPUT="$5"

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
"${BWTB3M}" inputtype=pacterm mem="${MEM}g" outputfilename=${INPUT}.pac.bwt ${INPUT}.pac
"${BWTB3MTOBWA}" ${INPUT}.pac.bwt ${INPUT}.bwt ${INPUT}.sa
"${BWA}" bwtupdate "${INPUT}.bwt"
rm -f "${INPUT}".pac.*
"${BWA}" fa2pac -f "${INPUT}"
```
