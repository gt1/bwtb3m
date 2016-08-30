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

Calling bwtb3m and options
--------------------------

usage: src/bwtb3m [options] <inputfile>

options:

```
* inputtype=[<bytestream>] (bytestream,compactstream,pac,pacterm,lz4,utf-8)
* outputfilename=[<bwtb3m_myers-mac-8.local_35627_1472591133.bwt>] (name of output .bwt file)
* sasamplingrate=[32] sampling rate for sampled suffix array
* isasamplingrate=[262144] sampling rate for sampled inverse suffix array
* mem=[2147483648] memory target (suffixes k,m and g are accepted)
* numthreads=[8] number of threads
* bwtonly=[0] compute BWT only (no sampled suffix array and reverse)
* tmpprefix=[bwtb3m_myers-mac-8.local_35627_1472591133] (prefix for tmp files)
* sparsetmpprefix=[tmpprefix] (prefix for sparse gap tmp files)
* copyinputtomemory=[0] (copy input file to memory)
* largelcpthres=[16384] (large LCP value threshold)
* verbose=[0] (verbosity level)
```

Output
------

bwtb3m computes the BWT of the given input file. Note that this refers to
the BWT as in the original definition (see M. Burrows and D. J. Wheeler: A Block-sorting Lossless
Data Compression Algorithm, Digital Research Report 124), in particular no
(implicit or explicit) terminator symbol is used. The input string is
considered as circular for the sake of comparisons. The output .bwt file can
be decoded using the bwtb3mdecoderl program or read using the class
libmaus2::huffman::RLDecoder in libmaus2.

If bwtonly=0, then bwtb3m computes a sampled suffix array by loading the final BWT to memory in
the form of a Huffman shaped wavelet tree. If bwtonly=1, then the program
only computes a hint file with suffix .preisa, which can be used for calling
the bwtcomputessa program to compute a sampled suffix array and sampled
inverse suffix array in external memory without loading the BWT into memory.

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
