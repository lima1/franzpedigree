FRANz 2.0.0
===========

A pedigree (family tree) reconstruction tool for natural populations.

Key Features:
-------------

* Allele frequency analyses:
     * Heterozygosities 
     * Polymorphic Information Content (PIC) 
     * Exact Hardy Weinberg Test
     * Test for Null Alleles

* Exclusion probabilities: Estimate the power of your marker suite. 

* FRANz is flexible, you CAN (but don't have to!) incorporate available prior
  knowledge:

     * Known sub-pedigrees (e.g. mother-offspring arcs) or fullsib groups
     * The sex and years of birth and death of individuals
     * Sampling locations (pairwise distances or coordinates)
     * Sampling rate (number of unsampled males and females)

  Most other tools assume that you know all or most of these.

* User friendly: FRANz uses this prior knowledge to generate the lists of 
  candidate parents internally. Once you have your data in a FRANz input
  file, you are almost done. A GUI for input file generation is available on
  our website.

* Can handle data with typing errors and missing alleles.

* Estimation of the number of unsampled candidate mothers and fathers.

* FRANz is fast: It was written for data sets with thousands of individuals.
  It supports modern CPUs with multiple cores.

* Available for Windows, Mac and Linux.

* Open Source: Obtain FRANz under the GPL Version 3.



Installation
============

Source Package
--------------

If possible, try to compile FRANz. Supporting the different versions of the
major operating systems is difficult, but the source code should always work.
If you have trouble compiling the code on Linux or Mac, please contact us.
Download the latest stable source code from SourceForge.net:

http://sourceforge.net/projects/franzpedigree/files/

Alternatively, get the newest code from github:

https://github.com/lima1/franzpedigree/tarball/master


Then unzip:

    tar xvfz lima1-franzpedigree-30eaf35.tar.gz
    cd lima1-franzpedigree-30eaf35


Then compile:

    autoreconf (only necessary if you use the latest code from github)
    ./configure 
    make check  (optional, may take a while)
    make install

For multicore CPUs (requires GCC >=4.2 or ICC).

    ./configure --enable-openmp
    make check (optional)
    make install



Mac OS X
========

For Leopard or newer, you can either download the installer or you can compile the
source package (this works with older Mac OS X Versions as well).

To compile the source package, you need the Apple Developer Tools. To enable
OpenMP, you have to compile FRANz with gcc-4.2 (or newer) or Intel's ICC. GCC
4.2 is available on Leopard:

    ./configure --enable-openmp CC=gcc-4.2
    make check (optional)
    make install

Make sure that you use gcc or icc for OpenMP; on newer Mac versions, you have
to install the command line tools in XCode 4.3 or newer.

On Tiger if you don't want to install a recent compiler (e.g. from
finkproject.org):

    ./configure
    make check (optional)
    make install

Windows
=======

If you want to install FRANz on Windows, you can use the installer available
from our website. Currently, only the single CPU version is available for
Windows. Linux and Mac are still the recommended platforms.

To compile the source package, you will need something like Cygwin (this
should give you the multi-core Version) or MinGW (there are experimental
releases with OpenMP support). Make sure that FRANz is compiled with -DWIN32.
Otherwise at least the progressbar will not work.


Additional configure Flags
==========================

    --enable-openmp      OpenMP support
    --enable-debug       Enables debug code and compiles with -g -O0
    --enable-assertions  Enables assertions. Warning: will make the program
                         really slow

Tutorial
========

Once you have your data in the input format, analyzing it should be easy.
There are some Perl scripts in extras/input to convert Excel (saved as
CSV-file) or migrate data. There is also an user friendly CSV to FRANz
converter available on the FRANz website.

Also check the manpage 
  
    man FRANz


Known Problems
==============

* Allele frequencies with partially missing genotypes. CERVUS ignores the rare
  cases with partially missing genotypes (e.g. 123.?) in the allele frequency 
  calculation and we followed this strategy. This can cause  problems for rare 
  alleles.


Troubleshooting
===============

Most of the problems are related to wrongly formatted input files. If you get
a segmentation fault, please if possible try to re-compile with --enable-debug.
Then start FRANz in the GNU debugger:

    $ gdb lima1-franzpedigree-30eaf35/scr/FRANz

Then run the program with the same set of paramaters:

    GNU gdb 6.3.50-20050815 (Apple version gdb-1515) (Sat Jan 15 08:33:48 UTC
    2011)
    Copyright 2004 Free Software Foundation, Inc.
    GDB is free software, covered by the GNU General Public License, and you are
    welcome to change it and/or distribute copies of it under certain conditions.
    Type "show copying" to see the conditions.
    There is absolutely no warranty for GDB.  Type "show warranty" for details.
    This GDB was configured as "x86_64-apple-darwin"...Reading symbols for shared
    libraries .. done

    (gdb) run input.dat --femrepro 1:10 --malerepro 1:10 --N 10 

After the segmentation fault, type backtrace and send us the output with your
bugreport:

    Program received signal EXC_BAD_ACCESS, Could not access memory.
    Reason: KERN_INVALID_ADDRESS at address: 0x0000000000000000
    0x00007fff86bddc00 in strlen ()
    (gdb) backtrace 



Feedback
========

Any comments, questions, critics or suggestions are gratefully received.  So
please don't hesitate to contact us! We would be happy to help. Your feedback
will help us improving this software.

markus@bioinf.uni-leipzig.de

Copyright
=========

Copyright 2008-2010. University of Leipzig (all code except parts of
src/hwe.c, src/exact.c, src/tap*, libdir/dcmt and src/uthash.h. See code for
details).

Some parts of dag.c are adapted from code presented in the book Algorithms in
C, Part 5 and parts of utils.c are taken from Numerical Recipes in C++. 
