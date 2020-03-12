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


bin_PROGRAMS=crac

CRAC_LINKS = crac-index

all-local: $(bin_PROGRAMS) $(CRAC_LINKS)


SUBDIRS = libReadsInfo libSSA libgzstream libReadsReader
if INCLUDED_GKARRAYS
  SUBDIRS += libGkArrays
endif
if INCLUDED_JELLYFISH
  SUBDIRS += libJellyfish-2.0
endif

crac_SOURCES =                    \
  basic.h const.h types.h         \
  Classifier.cpp Classifier.h	    \
  classifyTags.cpp classifyTags.h \
	ReadIndex.h                     \
	GkReadIndex.h										\
	Read.h												  \
	GkReadIndex.cpp GkReadIndex.h   \
	JellyReadIndex.cpp JellyReadIndex.h   \
  crac.cpp                        \
  Parameters.cpp Parameters.h     \
  SupportBreak.cpp SupportBreak.h \
  BreakList.cpp BreakList.h \
  CandidatBreak.cpp CandidatBreak.h \
  Support.cpp Support.h           \
  locateServer.cpp locateServer.h \
  utils.cpp utils.h \
  Bitset.cpp Bitset.h
crac_LDADD = @abs_builddir@/libSSA/libSSA.a @abs_builddir@/libReadsInfo/libReadsInfo.a @abs_builddir@/libgzstream/libgzstream.a @abs_builddir@/libReadsReader/libReadsReader.a 

#crac_index_SOURCES = @abs_srcdir@/libSSA/cracIndex.cpp
#crac_index_LDADD = @abs_builddir@/libSSA/libSSA.a

crac-index: @abs_builddir@/libSSA/crac-index
	$(LN_S) -f $^ $@

crac-emt: @abs_builddir@/libSSA/crac-emt
	$(LN_S) -f $^ $@


AM_LDFLAGS = -lm -lpthread -lhts -lSSA -lReadsInfo -lReadsReader -lgzstream -L@abs_builddir@/libSSA -L@abs_builddir@/libReadsInfo -L@abs_builddir@/libgzstream -L@abs_builddir@/libReadsReader  #-pg # for Gprof
AM_CPPFLAGS = -I@abs_top_srcdir@ -I@abs_srcdir@ -I@abs_srcdir@/libSSA -I@abs_srcdir@/libSSA/karkkainen_bwt -I@abs_srcdir@/libgzstream 

if INCLUDED_JELLYFISH
  AM_LDFLAGS += -L@abs_builddir@/libJellyfish-2.0/.libs
  AM_CPPFLAGS += -I@abs_srcdir@/libJellyfish-2.0/include
  AM_CXXFLAGS = -std=c++0x
endif

if INCLUDED_GKARRAYS
  AM_LDFLAGS += -L@abs_builddir@/libGkArrays/src/.libs
  AM_CPPFLAGS += -I@abs_srcdir@/libGkArrays/src/
if INCLUDED_PROGRESSBAR
  AM_LDFLAGS += -L@abs_builddir@/libGkArrays/src/libProgressBar-1.2.0/libProgressBar/.libs -lProgressBar
  AM_CPPFLAGS += -I@abs_srcdir@/libGkArrays/src/libProgressBar-1.2.0/
endif
endif

crac_LDFLAGS =

if HAVE_LIBGKARRAYS
  crac_LDFLAGS += -lGkArrays
endif
if HAVE_LIBJELLYFISH
  crac_LDFLAGS += -ljellyfish-2.0
endif

crac_LDFLAGS += $(AM_LDFLAGS)

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(EXTRA_PROGRAMS) $(CRAC_LINKS)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
