###############################################################################
#                                                                             #
#  Copyright © 2010-2013 -- IRB/INSERM                                        #
#                           (Institut de Recherches en Biothérapie /          #
#                           Institut National de la Santé et de la Recherche  #
#                           Médicale)                                         #
#                           LIFL/INRIA                                        #
#                           (Laboratoire d'Informatique Fondamentale de       #
#                           Lille / Institut National de Recherche en         #
#                           Informatique et Automatique)                      #
#                           LIRMM/CNRS                                        #
#                           (Laboratoire d'Informatique, de Robotique et de   #
#                           Microélectronique de Montpellier /                #
#                           Centre National de la Recherche Scientifique)     #
#                           LITIS                                             #
#                           (Laboratoire d'Informatique, du Traitement de     #
#                           l'Information et des Systèmes).                   #
#                                                                             #
#                                                                             #
#  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              #
#                   Mikaël SALSON    <mikael.salson@lifl.fr>                  #
#                   Thérèse COMMES   <commes@univ-montp2.fr>                  #
#                   Éric RIVALS      <eric.rivals@lirmm.fr>                   #
#                                                                             #
#  Programmeurs                                                               #
#      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              #
#                   Mikaël SALSON    <mikael.salson@lifl.fr>                  #
#  with additional contribution for the packaging of:	                      #
#                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               #
#                                                                             #
#  Contact:         CRAC list        <crac-bugs@lists.gforge.inria.fr>        #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce fichier fait partie du programme CRAC.                                  #
#                                                                             #
#  Crac est un outil d'analyse de données de RNA-Seq provenant des NGS.       #
#                                                                             #
#  Ce logiciel est régi  par la licence CeCILL  soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur  #
#  le site "http://www.cecill.info".                                          #
#                                                                             #
#  En contrepartie de l'accessibilité au code source et des droits de copie,  #
#  de modification et de redistribution accordés par cette licence, il n'est  #
#  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  #
#  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  #
#  titulaire des droits patrimoniaux et les concédants successifs.            #
#                                                                             #
#  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  #
#  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  #
#  développement  et à la reproduction du  logiciel par  l'utilisateur étant  #
#  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  #
#  manipuler et qui le réserve donc à des développeurs et des professionnels  #
#  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  #
#  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  #
#  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  #
#  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  #
#  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         #
#                                                                             #
#  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  #
#  pris connaissance  de la licence CeCILL,  et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This File is part of the CRAC program.                                     #
#                                                                             #
#  Crac is a tool to analyse RNA-Seq data provided by NGS.                    #
#                                                                             #
#  This software is governed by the CeCILL license under French law and       #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL     #
#  license as circulated by CEA, CNRS and INRIA at the following URL          #
#  "http://www.cecill.info".                                                  #
#                                                                             #
#  As a counterpart to the access to the source code and rights to copy,      #
#  modify and redistribute granted by the license, users are provided only    #
#  with a limited warranty and the software's author, the holder of the       #
#  economic rights, and the successive licensors have only limited            #
#  liability.                                                                 #
#                                                                             #
#  In this respect, the user's attention is drawn to the risks associated     #
#  with loading, using, modifying and/or developing or reproducing the        #
#  software by the user in light of its specific status of free software,     #
#  that may mean that it is complicated to manipulate, and that also          #
#  therefore means that it is reserved for developers and experienced         #
#  professionals having in-depth computer knowledge. Users are therefore      #
#  encouraged to load and test the software's suitability as regards their    #
#  requirements in conditions enabling the security of their systems and/or   #
#  data to be ensured and, more generally, to use and operate it in the same  #
#  conditions as regards security.                                            #
#                                                                             #
#  The fact that you are presently reading this means that you have had       #
#  knowledge of the CeCILL license and that you accept its terms.             #
#                                                                             #
###############################################################################


SUBDIRS = src man scripts test_f

#  uncomment the following if crac requires the math library
#crac_LDADD=-lm

EXTRA_DIST= @PACKAGE@.lsm.in @PACKAGE@.spec.in control.in AUTHORS COPYING ChangeLog INSTALL NEWS README
# dist-hook:
# 	rm -rf `find $(distdir) -type d -name .svn`

valgrindcheck: 
	LAUNCHER='eval valgrind --tool=memcheck --suppressions=$(abs_top_srcdir)/test_f/valgrind.supp --leak-check=full --show-reachable=yes --undef-value-errors=yes --track-origins=yes --xml=yes --xml-file=$$(mktemp --tmpdir=. --suffix=.valgrind)' make check


deb: @PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb



always-build:

@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb: always-build
	@if [ "@DEB_BUILDER@" == "" ]; then \
	  echo "You need to install a debian package builder and to reconfigure this package.";\
	  exit 1;\
	fi
	$(MAKE) dist
	$(am__remove_distdir)
	GZIP=$(GZIP_ENV) gunzip -c $(distdir).tar.gz | $(AMTAR) xf -
	mkdir $(distdir)/_build
	mkdir -p $(distdir)/_inst/_deb/$(prefix)
	chmod a-w $(distdir)
	dc_install_base=`$(am__cd) $(distdir)/_inst/_deb/$(prefix) && pwd | sed -e 's,^[^:\\/]:[\\/],/,'` \
	  && dc_destdir="$${TMPDIR-/tmp}/am-dc-$$$$/" \
	  && cd $(distdir)/_build \
	  && localedir=$(datadir)/locale \
	  && ../configure --srcdir=.. --prefix="$$dc_install_base" \
            --enable-static --disable-shared $(DEB_OPT) \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) localedir=$$localedir \
	  && $(MAKE) $(AM_MAKEFLAGS) dvi \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && cd src/libGkArrays \
	  && $(MAKE) $(AM_MAKEFLAGS) uninstall \
	  && cd ../../ \
	  && cd ../_inst/_deb \
	  && find . -depth -type d -empty -exec rmdir {} \; \
	  && find . -type f -exec md5sum {} \; | sed "s,\.\/,,g" > ../md5sums \
          && mkdir DEBIAN \
	  && mv ../md5sums DEBIAN/ \
	  && sed "s,^\(Installed-Size:.*\)0$$,\1$$(du -sk $$dc_install_base/ | cut -f1)," ../../_build/control > DEBIAN/control \
	  && cd ../ \
	  && dpkg -b _deb $@ \
	  && mv $@ ../../ \
	  && cd ../_build/src/libGkArrays \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && cd ../../ \
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

rpm: always-build
	@if [ "@RPM_BUILDER@" == "" ]; then \
	  echo "You need to install a rpm package builder and to reconfigure this package.";\
	  exit 1;\
	fi
	$(MAKE) dist
	$(am__remove_distdir)
	- rm -rf $(distdir)/_build/rpm/BUILD/
	- rm -rf $(distdir)/_build/rpm/RPMS/
	- rm -rf $(distdir)/_build/rpm/SRPMS/
	- rm -rf $(distdir)/_build/rpm/SOURCES/
	mkdir -p $(distdir)/_build/rpm/SOURCES
	cp $(distdir).tar.gz $(distdir)/_build/rpm/SOURCES/
	@RPM_BUILDER@ --quiet --nodeps --define="_topdir "$$(pwd)"/$(distdir)/_build/rpm" \
	-bb crac.spec
	mv $(distdir)/_build/rpm/RPMS/*/*.rpm .

MOSTLYCLEANFILES=*~
CLEANFILES=*~
DISTCLEANFILES=*~ control $(distdir).tar.gz
DISTCLEANFILES+=@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb *.rpm
MAINTAINERCLEANFILES=*~

ACLOCAL_AMFLAGS = -I config/m4
