SUBDIRS = $(PACKAGE)

if BUILD_DOC
  SUBDIRS += doc
endif

DISTSUBDIRS = $(PACKAGE)
EXTRA_DIST=$(PACKAGE).lsm.in $(PACKAGE).spec.in control.in AUTHORS COPYING ChangeLog INSTALL NEWS README

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@-@PACKAGE_VERSION@.pc

@PACKAGE@-@PACKAGE_VERSION@.pc: @PACKAGE@.pc
	cp "$^" "$@"

deb: @PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb

rpm: rpm-bin rpm-src
rpm-bin: @PACKAGE@-@PACKAGE_VERSION@-1.@ARCH@.rpm
rpm-src: @PACKAGE@-@PACKAGE_VERSION@-1.src.rpm

ACLOCAL_AMFLAGS=-I m4

@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb:
	@if [ "@DEB_BUILDER@" == "" ]; then \
	  echo "You need to install a debian package builder and to reconfigure this package.";\
	  exit 1;\
	fi
	make dist
	$(am__remove_distdir)
	GZIP=$(GZIP_ENV) gunzip -c $(distdir).tar.gz | $(AMTAR) xf -
	mkdir $(distdir)/_build
	mkdir -p $(distdir)/_inst/_deb/$(prefix)
	chmod a-w $(distdir)
	dc_install_base=`$(am__cd) $(distdir)/_inst/_deb/$(prefix) && pwd | sed -e 's,^[^:\\/]:[\\/],/,'` \
	  && dc_destdir="$${TMPDIR-/tmp}/am-dc-$$$$/" \
	  && cd $(distdir)/_build \
	  && localedir=$(datadir)/locale \
	  && HACKED_PREFIX=$(prefix) ../configure --srcdir=.. --prefix="$$dc_install_base" \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) localedir=$$localedir \
	  && $(MAKE) $(AM_MAKEFLAGS) dvi \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && cd ../_inst/_deb \
	  && find . -type f -exec md5sum {} \; | sed "s,\.\/,,g" > ../md5sums \
          && mkdir DEBIAN \
	  && mv ../md5sums DEBIAN/ \
	  && cp ../../_build/control DEBIAN/ \
	  && cp ../../_build/postinst DEBIAN/ \
	  && cp ../../_build/prerm DEBIAN/ \
	  && sed "s/@INSTALLED_SIZE@/`du -ck . | tail -1 | cut -f1`/" ../../_build/control > DEBIAN/control \
	  && cd ../ \
	  && dpkg -b _deb $@ \
	  && mv $@ ../../ \
	  && rm -f "$$dc_install_base/DEBIAN/control" \
	  && rm -f "$$dc_install_base/DEBIAN/prerm" \
	  && rm -f "$$dc_install_base/DEBIAN/postinst" \
	  && cd ../_build \
	  && $(MAKE) $(AM_MAKEFLAGS) uninstall \
	  && $(MAKE) $(AM_MAKEFLAGS) distuninstallcheck_dir="$$dc_install_base" \
	        distuninstallcheck \
	  && chmod -R a-w "$$dc_install_base" \
	  && ({ \
	       (cd ../.. && $(mkinstalldirs) "$$dc_destdir") \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" install \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" uninstall \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" \
	            distuninstallcheck_dir="$$dc_destdir" distuninstallcheck; \
	      } || { rm -rf "$$dc_destdir"; exit 1; }) \
	  && rm -rf "$$dc_destdir" \
	  && $(MAKE) $(AM_MAKEFLAGS) distcleancheck
	$(am__remove_distdir)
	@echo "$@ is ready for distribution" | \
	  sed 'h;s/./=/g;p;x;p;x'

rpm-all:
	@if [ "@RPM_BUILDER@" == "" ]; then \
	  echo "You need to install a rpm package builder and to reconfigure this package.";\
	  exit 1;\
	fi
	make dist
	$(am__remove_distdir)
	mkdir $(distdir)
	mkdir -p $(distdir)/BUILD
	mkdir -p $(distdir)/RPMS
	mkdir -p $(distdir)/SOURCES
	mkdir -p $(distdir)/SPECS
	mkdir -p $(distdir)/SRPMS
	cp $(distdir).tar.gz $(distdir)/SOURCES

@PACKAGE@-@VERSION@-1.@ARCH@.rpm : @PACKAGE@.spec
	make rpm-all
	@RPM_BUILDER@ --clean -bb --sign --quiet $<
	mv $(distdir)/RPMS/@ARCH@/$@ .
	$(am__remove_distdir)
	@echo "$@ is ready for distribution" | \
	  sed 'h;s/./=/g;p;x;p;x'

@PACKAGE@-@VERSION@-1.src.rpm : @PACKAGE@.spec
	make rpm-all
	@RPM_BUILDER@ --clean -bs --sign --quiet $<
	mv $(distdir)/SRPMS/$@ .
	$(am__remove_distdir)
	@echo "$@ is ready for distribution" | \
	  sed 'h;s/./=/g;p;x;p;x'


install-data-hook:
	$(LN_S) -f @PACKAGE@-@VERSION@.pc $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc

uninstall-hook:
	test "`readlink $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc`" != @PACKAGE@-@VERSION@.pc \
        || $(RM) $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc

MOSTLYCLEANFILES=*~
CLEANFILES=*~
DISTCLEANFILES=*~ control $(distdir).tar.gz
DISTCLEANFILES+=@PACKAGE@*-@PACKAGE_VERSION@.pc
DISTCLEANFILES+=@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb
DISTCLEANFILES+=@PACKAGE@-@PACKAGE_VERSION@-1.@ARCH@.rpm
DISTCLEANFILES+=@PACKAGE@-@PACKAGE_VERSION@-1.src.rpm
MAINTAINERCLEANFILES=*~
