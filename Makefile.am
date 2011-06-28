ACLOCAL_AMFLAGS = -I m4

SUBDIRS= libdir examples man src t extras

EXTRA_DIST = $(top_srcdir)/t/testdata/*.dat $(top_srcdir)/t/testdata/*.dist $(top_srcdir)/t/testdata/*.coord $(top_srcdir)/extras/input
LCOV = lcov
GENHTML = genhtml

cov: clean-cov
	mkdir $(top_builddir)/t/coverage
	$(LCOV) -c -o $(top_builddir)/t/coverage/franz.info.tmp -d $(top_srcdir)/t
	$(LCOV) -r $(top_builddir)/t/coverage/franz.info.tmp -o $(top_builddir)/t/coverage/franz.info 
	rm $(top_builddir)/t/coverage/franz.info.tmp
	$(GENHTML) -t "franz" -o $(top_builddir)/t/coverage --legend $(top_builddir)/t/coverage/franz.info

clean-cov:
	rm -rf $(top_builddir)/t/coverage

maintainer-clean-local:
	rm -f $(top_builddir)/t/*.gcda
	rm -f $(top_builddir)/t/*.gcno
	rm -f $(top_builddir)/src/*.gcda
	rm -f $(top_builddir)/src/*.gcno
	rm -f $(top_builddir)/libdir/dcmt/lib/*.gcda
	rm -f $(top_builddir)/libdir/dcmt/lib/*.gcno
