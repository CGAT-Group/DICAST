library_includedir=@includedir@/@PACKAGE@-@VERSION@
library_include_HEADERS = progressBar.h
lib_LTLIBRARIES = libProgressBar.la
libProgressBar_la_SOURCES = progressBar.h progressBar.cpp
libProgressBar_la_configdir = @libdir@/@PACKAGE@-@VERSION@/include
libProgressBar_la_config_DATA = config.h
libProgressBar_la_LDFLAGS = -version-info @LIB_VERSION@ -release @VERSION@ $(AM_LDFLAGS)

noinst_PROGRAMS=testProgressBar
testProgressBar_SOURCES=testPB.cpp
testProgressBar_LDADD=libProgressBar.la
testProgressBar_LDFLAGS=-lProgressBar $(AM_LDFLAGS)

install-data-hook:
	$(LN_S) -f @PACKAGE@-@VERSION@ $(DESTDIR)@includedir@/@PACKAGE@

uninstall-hook:
	test "`readlink $(DESTDIR)@includedir@/@PACKAGE@`" != @PACKAGE@-@VERSION@ \
        || $(RM) -f $(DESTDIR)@includedir@/@PACKAGE@

MOSTLYCLEANFILES=*~
CLEANFILES=*~
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
