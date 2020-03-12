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
# $Id: Makefile.am 763 2011-07-10 20:48:29Z salson $
#
###############################################################################

SUBDIRS = kseq sux

if INCLUDED_PROGRESSBAR
  SUBDIRS += libProgressBar-1.2.0
endif

KSEQ_DIR=kseq
KSEQ_HDRS=$(KSEQ_DIR)/kseq.h

SUX_DIR=sux
SUX_SRCS=$(SUX_DIR)/rank9b.cpp
SUX_HDRS=\
	$(SUX_DIR)/macros.h \
	$(SUX_DIR)/popcount.h \
	$(SUX_DIR)/rank9b.h

AM_LDFLAGS = -lpthread 
AM_CPPFLAGS = -I@abs_top_srcdir@
# if INCLUDED_DIVSUFSORT
#   AM_LDFLAGS += -L@abs_top_builddir@/libdivsufsort/lib/.libs
#   AM_CPPFLAGS += -I@abs_top_builddir@/libdivsufsort/include
# endif


includedir = @includedir@/@PACKAGE@-@VERSION@
libdir=@libdir@
include_HEADERS = gkArrays.h solArray.h gkArraysTypes.h gkArraysBuilder.h gkSABuilder.h gkCFPSBuilder.h gkIFABuilder.h readsReader.h bitvector.hxx largebinvect.hxx bitvector.cxx largebinvect.cxx binvect.hxx
noinst_HEADERS = utils.h basic.h const.h $(KSEQ_HDRS) $(SUX_HDRS)

lib_LTLIBRARIES = libGkArrays.la
libGkArrays_ladir = @abs_srcdir@
libGkArrays_la_SOURCES = gkArrays.cpp solArray.cpp utils.cpp gkArraysBuilder.cpp gkSABuilder.cpp gkCFPSBuilder.cpp gkIFABuilder.cpp  readsReader.cpp largebinvect.cxx bitvector.cxx $(SUX_SRCS)
libGkArrays_la_configdir = @libdir@/@PACKAGE@-@VERSION@/include
libGkArrays_la_config_DATA = @abs_top_builddir@/config/config.h
libGkArrays_la_LDFLAGS = -version-info @LIB_VERSION@ -release @VERSION@ $(AM_LDFLAGS) # @LIBDIVSUFSORT@

if INCLUDED_PROGRESSBAR
  AM_LDFLAGS += -L@abs_builddir@/libProgressBar-1.2.0/libProgressBar/.libs -lProgressBar
  AM_CPPFLAGS += -I@abs_srcdir@/libProgressBar-1.2.0/
endif

noinst_PROGRAMS = buildTables
buildTables_SOURCES = buildTables.cpp
buildTables_LDADD=libGkArrays.la
buildTables_LDFLAGS = -lGkArrays $(AM_LDFLAGS)

install-data-hook:
	$(LN_S) -f @PACKAGE@-@VERSION@ $(DESTDIR)@includedir@/@PACKAGE@

uninstall-hook:
	test "`readlink $(DESTDIR)@includedir@/@PACKAGE@`" != @PACKAGE@-@VERSION@ \
        || $(RM) -f $(DESTDIR)@includedir@/@PACKAGE@

MOSTLYCLEANFILES=*~
CLEANFILES=*~ *.gcda *.gcno
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
