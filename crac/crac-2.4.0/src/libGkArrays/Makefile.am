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
#                   Thierry LECROQ   <thierry.lecroq@univ-rouen.fr>           #
#                   Martine LÉONARD  <Martine.Leonard@univ-rouen.fr>          #
#                   Éric RIVALS      <eric.rivals@lirmm.fr>                   #
#                                                                             #
#  Programmeurs                                                               #
#      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              #
#                   Mikaël SALSON    <mikael.salson@lifl.fr>                  #
#  with additional contribution for the packaging of:	                      #
#                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               #
#                                                                             #
#  Contact:         Gk-Arrays list   <crac-gkarrays@lists.gforge.inria.fr>    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce fichier fait partie de la librairie Gk-arrays.                          #
#                                                                             #
#  La librairie Gk-arrays  a  pour objectif d'indexer de grands ensembles de  #
#  lectures de séquences issues du séquençage haut-débit.                     #
#                                                                             #
#  Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA sur  #
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
#  pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This File is part of the Gk-arrays library.                                #
#                                                                             #
#  The Gk-arrays library aims at indexing k-factors from a huge set of        #
#  sequencing reads.                                                          #
#                                                                             #
#  This software is governed by the CeCILL-C license under French law and     #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL-C   #
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
#  knowledge of the CeCILL-C license and that you accept its terms.           #
#                                                                             #
###############################################################################
#
# $Id: Makefile.am 778 2011-07-11 08:05:38Z salson $
#
###############################################################################

SUBDIRS=
SUBDIRS += src
if BUILD_DOC
  SUBDIRS += doc
endif
SUBDIRS += tests 

EXTRA_DIST=control.in AUTHORS COPYING ChangeLog INSTALL NEWS README

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@-@PACKAGE_VERSION@.pc

@PACKAGE@-@PACKAGE_VERSION@.pc: @PACKAGE@.pc
	cp "$^" "$@"

ACLOCAL_AMFLAGS = -I config/m4

deb: @PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb

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
	  && HACKED_PREFIX=$(prefix) ../configure --srcdir=.. --prefix="$$dc_install_base" CC="$(CC)" CXX="$(CXX)" \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) localedir=$$localedir \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && cd ../_inst/_deb \
	  && find . -type f -exec md5sum {} \; | sed "s,\.\/,,g" > ../md5sums \
          && mkdir DEBIAN \
	  && mv ../md5sums DEBIAN/ \
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
	      } \
	      || { rm -rf "$$dc_destdir"; exit 1; }) \
	  && rm -rf "$$dc_destdir" \
	  && $(MAKE) $(AM_MAKEFLAGS) distcleancheck
	$(am__remove_distdir)
	@echo "$@ is ready for distribution" | \
	  sed 'h;s/./=/g;p;x;p;x'


install-data-hook:
	$(LN_S) -f @PACKAGE@-@VERSION@.pc $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc

uninstall-hook:
	test "`readlink $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc`" != @PACKAGE@-@VERSION@.pc \
        || $(RM) $(DESTDIR)@libdir@/pkgconfig/@PACKAGE@.pc

MOSTLYCLEANFILES=*~
CLEANFILES=*~ *.gcda *.gcno
DISTCLEANFILES=*~ control $(distdir).tar.gz
DISTCLEANFILES+=@PACKAGE@*-@PACKAGE_VERSION@.pc
DISTCLEANFILES+=@PACKAGE@_@PACKAGE_VERSION@_@ARCH@.deb
MAINTAINERCLEANFILES=*~
