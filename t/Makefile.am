# binaries

LDADD = ../libdir/dcmt/lib/libdcmt.a
#AUTOMAKE_OPTIONS = subdir-objects
AM_CFLAGS = -I../libdir/dcmt/lib

check_PROGRAMS = dag_tests dag_tests2 io_tests prob_tests ll_tests utils_tests hwe_tests mcmc_tests cervus_tests mlrelate_tests

commonSOURCES = tap.c  $(top_srcdir)/src/vtprogressbar.c $(top_srcdir)/src/dag.c $(top_srcdir)/src/listdag.c $(top_srcdir)/src/dataio.c $(top_srcdir)/src/freq.c $(top_srcdir)/src/listparentage.c $(top_srcdir)/src/options.c $(top_srcdir)/src/prob.c $(top_srcdir)/src/hwe.c $(top_srcdir)/src/utils.c  $(top_srcdir)/src/null.c $(top_srcdir)/src/ibd.c $(top_srcdir)/src/lod.c $(top_srcdir)/src/pedigree.c $(top_srcdir)/src/genotype.c $(top_srcdir)/src/pvalue.c $(top_srcdir)/src/exact.c




commonHEADERS = tap.h $(top_srcdir)/src/macros.h $(top_srcdir)/src/global.h $(top_srcdir)/src/ll.h $(top_srcdir)/src/vtprogressbar.h $(top_srcdir)/src/dag.h $(top_srcdir)/src/listdag.h $(top_srcdir)/src/listparentage.h $(top_srcdir)/src/freq.h $(top_srcdir)/src/dataio.h $(top_srcdir)/src/options.h $(top_srcdir)/src/prob.h  $(top_srcdir)/src/hwe.h $(top_srcdir)/src/utils.h  $(top_srcdir)/src/null.h $(top_srcdir)/src/ibd.h  $(top_srcdir)/src/lod.h  $(top_srcdir)/src/pedigree.h  $(top_srcdir)/src/genotype.h $(top_srcdir)/src/pvalue.h $(top_srcdir)/src/exact.h





io_tests_SOURCES = io_tests.c $(commonSOURCES)\
				   $(commonHEADERS)

prob_tests_SOURCES = prob_tests.c $(top_srcdir)/src/mcmc.c $(top_srcdir)/src/sim.c $(commonSOURCES)\
				   $(top_srcdir)/src/mcmc.h $(top_srcdir)/src/sim.h $(commonHEADERS)

mcmc_tests_SOURCES = mcmc_tests.c $(top_srcdir)/src/mcmc.c $(commonSOURCES)\
				   $(top_srcdir)/src/mcmc.h $(commonHEADERS)

cervus_tests_SOURCES = cervus_tests.c $(top_srcdir)/src/sim.c $(commonSOURCES)\
					   $(top_srcdir)/src/sim.h $(commonHEADERS)

mlrelate_tests_SOURCES = mlrelate_tests.c $(top_srcdir)/src/sim.c $(commonSOURCES)\
					   $(top_srcdir)/src/sim.h $(commonHEADERS)

hwe_tests_SOURCES = hwe_tests.c $(commonSOURCES) \
				    $(commonHEADERS)

dag_tests_SOURCES = dag_tests.c $(commonSOURCES)\
                  $(top_srcdir)/src/pedigree.h $(commonHEADERS)

dag_tests2_SOURCES = $(dag_tests_SOURCES)

ll_tests_SOURCES = ll_tests.c tap.c $(top_srcdir)/src/listdag.c $(top_srcdir)/src/vtprogressbar.c \
				   tap.h $(top_srcdir)/src/ll.h $(top_srcdir)/src/listdag.h $(top_srcdir)/src/vtprogressbar.h 

utils_tests_SOURCES = utils_tests.c tap.c $(top_srcdir)/src/utils.c $(top_srcdir)/src/pvalue.c  $(top_srcdir)/src/vtprogressbar.c \
				   tap.h $(top_srcdir)/src/utils.h $(top_srcdir)/src/global.h $(top_srcdir)/src/macros.h  $(top_srcdir)/src/pvalue.h $(top_srcdir)/src/vtprogressbar.h 


dag_tests_CFLAGS = -g -I $(top_srcdir)/src
dag_tests2_CFLAGS = -g -I $(top_srcdir)/src -DADJLIST

ll_tests_CFLAGS  = $(dag_tests_CFLAGS)

utils_tests_CFLAGS  = $(dag_tests_CFLAGS)
hwe_tests_CFLAGS  = $(dag_tests_CFLAGS)

io_tests_CFLAGS = $(dag_tests_CFLAGS) -DFILENAME='"$(top_srcdir)/t/testdata/EsCastel.dat"' \
				  -DFILENAMEMD='"$(top_srcdir)/examples/simpsons/simpsons.dat"'\
				  -DFILENAMEMOTHERS='"$(top_srcdir)/examples/simpsons/simpsons.mothers"'\
				  -DFILENAMEPOPS='"$(top_srcdir)/t/testdata/123.dat"'\
				  -DFILENAMEPOPSDIST='"$(top_srcdir)/t/testdata/123.dist"'\
				  -DFILENAMEPOPSCOORD='"$(top_srcdir)/t/testdata/123.coord"'\
				  -DFILENAMEPMDIST='"$(top_srcdir)/t/testdata/penaeus_monodon.dist"'\
				  -DFILENAMEPM='"$(top_srcdir)/t/testdata/penaeus_monodon.dat"'


mcmc_tests_CFLAGS = $(dag_tests_CFLAGS) -DFILENAMEMD='"$(top_srcdir)/examples/simpsons/simpsons.dat"'

prob_tests_CFLAGS = $(io_tests_CFLAGS)

cervus_tests_CFLAGS = $(dag_tests_CFLAGS)  -DFILENAME='"$(top_srcdir)/examples/penaeus_monodon/input/penaeus_monodon.dat"'\
   -DFILENAMEPED='"$(top_srcdir)/examples/penaeus_monodon/input/penaeus_monodon.mothers"'

mlrelate_tests_CFLAGS = $(cervus_tests_CFLAGS)
TESTS = dag_tests dag_tests2 io_tests prob_tests ll_tests utils_tests hwe_tests mcmc_tests cervus_tests mlrelate_tests
