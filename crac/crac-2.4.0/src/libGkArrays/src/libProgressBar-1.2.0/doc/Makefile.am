#  build and install the .info pages
info_TEXINFOS  = libProgressBar.texinfo
libProgressBar_TEXINFOS = lgpl.texinfo

#  install the man pages
man_MANS = libProgressBar.3

HTMLS  = html
DOCDIR = $(datadir)/$(PACKAGE)-$(VERSION)/
DOCS   = $(HTMLS) $(PSS) $(PDFS)

AM_MAKEINFOFLAGS = --enable-encoding
TEXI2DVI         = texi2dvi -E

all-local: info ps pdf html

install-data-local: $(DOCS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(DOCDIR)
	@for f in $^; do \
	  echo "test -f $$f && $(INSTALL_DATA) $$f $(DESTDIR)$(DOCDIR) || cp -R $$f $(DESTDIR)$(DOCDIR)"; \
	  test -f $$f && $(INSTALL_DATA) $$f $(DESTDIR)$(DOCDIR) || cp -R $$f $(DESTDIR)$(DOCDIR); \
	done

uninstall-local: $(DOCS)
	@$(NORMAL_UNINSTALL)
	@for f in $^; do \
	  echo "rm -f $(DESTDIR)$(DOCDIR)$$f"; \
	  rm -rf $(DESTDIR)$(DOCDIR)$$f; \
	done
	@rmdir $(DESTDIR)$(DOCDIR) 2> /dev/null && \
	  echo "$(DESTDIR)$(DOCDIR) is removed." || \
	  echo "$(DESTDIR)$(DOCDIR) is not empty or doesn't exist. So not removed."

mostlyclean-local:
	rm -rf $(DOCS)

html: libProgressBar.texinfo
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) $(MAKEINFOFLAGS) -I $(srcdir) \
	 --html --number-sections `test -f '$<' || echo '$(srcdir)/'`$< -o $@

EXTRA_DIST=$(man_MANS)

MOSTLYCLEANFILES=$(man_MANS) *~
CLEANFILES=*~
DISTCLEANFILES=*~ $(PACKAGE).dvi
MAINTAINERCLEANFILES=*~
