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
------------------------

bwtb3m needs libmaus [https://github.com/gt1/libmaus] . When libmaus
is installed in ${LIBMAUSPREFIX} then bwtb3m can be compiled and
installed in ${HOME}/bwtb3m using

	- autoreconf -i -f
	- ./configure --with-libmaus=${LIBMAUSPREFIX} \
		--prefix=${HOME}/bwtb3m
	- make install
