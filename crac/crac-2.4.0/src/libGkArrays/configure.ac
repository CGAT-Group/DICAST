dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2010-2013 -- IRB/INSERM
dnl                          (Institut de Recherches en Biothérapie /
dnl                          Institut National de la Santé et de la Recherche
dnl                          Médicale)
dnl                          LIFL/INRIA
dnl                          (Laboratoire d'Informatique Fondamentale de
dnl                          Lille / Institut National de Recherche en
dnl                          Informatique et Automatique)
dnl                          LIRMM/CNRS
dnl                          (Laboratoire d'Informatique, de Robotique et de
dnl                          Microélectronique de Montpellier /
dnl                          Centre National de la Recherche Scientifique)
dnl                          LITIS
dnl                          (Laboratoire d'Informatique, du Traitement de
dnl                          l'Information et des Systèmes).
dnl
dnl Auteurs/Authors:
dnl                   Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
dnl                   Mikaël SALSON    <mikael.salson@lifl.fr>
dnl                   Thierry LECROQ   <thierry.lecroq@univ-rouen.fr>
dnl                   Martine LÉONARD  <Martine.Leonard@univ-rouen.fr>
dnl                   Éric RIVALS      <eric.rivals@lirmm.fr>
dnl
dnl Programmeurs
dnl      /Progammers:
dnl                   Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
dnl                   Mikaël SALSON    <mikael.salson@lifl.fr>
dnl  with additional contribution for the packaging of:
dnl                   Alban MANCHERON  <alban.mancheron@lirmm.fr>
dnl
dnl  Contact:         Gk-Arrays list   <crac-gkarrays@lists.gforge.inria.fr>
dnl -------------------------------------------------------------------------
dnl
dnl Ce fichier fait partie de la librairie Gk-arrays.
dnl
dnl La librairie Gk-arrays  a  pour objectif d'indexer de grands ensembles de
dnl lectures de séquences issues du séquençage haut-débit.
dnl
dnl Ce logiciel est régi par la licence CeCILL-C soumise au droit français et
dnl respectant les principes  de diffusion des logiciels libres.  Vous pouvez
dnl utiliser, modifier et/ou redistribuer ce programme sous les conditions de
dnl la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA sur
dnl le site "http://www.cecill.info".
dnl
dnl En contrepartie de l'accessibilité au code source et des droits de copie,
dnl de modification et de redistribution accordés par cette licence, il n'est
dnl offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
dnl seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le
dnl titulaire des droits patrimoniaux et les concédants successifs.
dnl
dnl À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques
dnl associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au
dnl développement  et à la reproduction du  logiciel par  l'utilisateur étant
dnl donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à
dnl manipuler et qui le réserve donc à des développeurs et des professionnels
dnl avertis  possédant  des  connaissances  informatiques  approfondies.  Les
dnl utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du
dnl logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la
dnl sécurité de leurs systêmes et ou de leurs données et,  plus généralement,
dnl à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
dnl
dnl Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez
dnl pris connaissance de la licence CeCILL-C, et que vous en avez accepté les
dnl termes.
dnl
dnl -------------------------------------------------------------------------
dnl
dnl This File is part of the Gk-arrays library.
dnl
dnl The Gk-arrays library aims at indexing k-factors from a huge set of
dnl sequencing reads
dnl
dnl This software is governed by the CeCILL-C license under French law and
dnl abiding by the rules of distribution of free software. You can use,
dnl modify and/ or redistribute the software under the terms of the CeCILL-C
dnl license as circulated by CEA, CNRS and INRIA at the following URL
dnl "http://www.cecill.info".
dnl
dnl As a counterpart to the access to the source code and rights to copy,
dnl modify and redistribute granted by the license, users are provided only
dnl with a limited warranty and the software's author, the holder of the
dnl economic rights, and the successive licensors have only limited
dnl liability.
dnl
dnl In this respect, the user's attention is drawn to the risks associated
dnl with loading, using, modifying and/or developing or reproducing the
dnl software by the user in light of its specific status of free software,
dnl that may mean that it is complicated to manipulate, and that also
dnl therefore means that it is reserved for developers and experienced
dnl professionals having in-depth computer knowledge. Users are therefore
dnl encouraged to load and test the software's suitability as regards their
dnl requirements in conditions enabling the security of their systems and/or
dnl data to be ensured and, more generally, to use and operate it in the same
dnl conditions as regards security.
dnl
dnl The fact that you are presently reading this means that you have had
dnl knowledge of the CeCILL-C license and that you accept its terms.

dnl $Id: configure.ac 763 2011-07-10 20:48:29Z salson $

AC_PREREQ(2.64)

m4_define([VERSION_MAJOR], [2])
m4_define([VERSION_MINOR], [1])
m4_define([VERSION_MICRO], [1])

AC_INIT([libGkArrays],
        [VERSION_MAJOR.VERSION_MINOR.VERSION_MICRO],
        [crac-gkarrays@lists.gforge.inria.fr],
        [libGkArrays],
        [http://crac.gforge.inria.fr/gkarrays])
PACKAGE_DESCRIPTION="Library for building and using GkArrays"
AC_SUBST([PACKAGE_DESCRIPTION])

AC_CONFIG_SRCDIR([src/gkArrays.h])
AC_REVISION([$Revision: 763 $])
AC_COPYRIGHT([
*******************************************************************************
*                                                                             *
*  Copyright(C) 2010-2013 --IRB/INSERM                                        *
*                           (Institut de Recherches en Biothérapie /          *
*                           Institut National de la Santé et de la Recherche  *
*                           Médicale)                                         *
*                           LIFL/INRIA                                        *
*                           (Laboratoire d'Informatique Fondamentale de       *
*                           Lille / Institut National de Recherche en         *
*                           Informatique et Automatique)                      *
*                           LIRMM/CNRS                                        *
*                           (Laboratoire d'Informatique, de Robotique et de   *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique)     *
*                           LITIS                                             *
*                           (Laboratoire d'Informatique, du Traitement de     *
*                           l'Information et des Systèmes).                   *
*                                                                             *
*  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*                   Thierry LECROQ   <thierry.lecroq@univ-rouen.fr>           *
*                   Martine LÉONARD  <Martine.Leonard@univ-rouen.fr>          *
*                   Éric RIVALS      <eric.rivals@lirmm.fr>                   *
*                                                                             *
*  Programmeurs                                                               *
*      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*  with additional contribution:                                              *
*                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               *
*                                                                             *
*  Contact:         Gk-Arrays list   <crac-gkarrays@lists.gforge.inria.fr>    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce fichier fait partie de la librairie Gk-arrays.                          *
*                                                                             *
*  La librairie Gk-arrays  a  pour objectif d'indexer de grands ensembles de  *
*  lectures de séquences issues du séquençage haut-débit.                     *
*                                                                             *
*  Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This File is part of the Gk-arrays library.                                *
*                                                                             *
*  The Gk-arrays library aims at indexing k-factors from a huge set of        *
*  sequencing reads.                                                          *
*                                                                             *
*  This software is governed by the CeCILL-C license under French law and     *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL-C   *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
******************************************************************************/
])


dnl The most recent interface number that this library implements.
LIB_CURRENT=2
dnl The implementation number of the current interface.
LIB_REVISION=0
dnl The difference between the newest and oldest interfaces that this library implements.
LIB_AGE=1
dnl LIB_VERSION is not the release number.
dnl Please read libtool documentation about versionning
LIB_VERSION=$LIB_CURRENT:$LIB_REVISION:$LIB_AGE
dnl Ending numbers of library will be $((LIB_CURRENT-LIB_AGE)).$LIB_AGE.$LIB_REVISION

AC_SUBST([VERSION])
AC_SUBST([LIB_VERSION])
AC_SUBST([LIB_VERSION_EXT], ["VERSION_MAJOR"_"VERSION_MINOR"_"VERSION_MICRO"])

AS_IF([test x"$CFLAGS" == "x"], [CFLAGS=""])
AS_IF([test x"$CPPFLAGS" == "x"], [CPPFLAGS=""])
AS_IF([test x"$CXXFLAGS" == "x"], [CXXFLAGS=""])
AS_IF([test x"$LDFLAGS" == "x"], [LDFLAGS=""])

AC_CONFIG_AUX_DIR([config])
dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config/config.h])
dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([config/m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_CHECK_PROG([TEXI2DVI], [texi2dvi], [texi2dvi])
AC_CHECK_PROG([MAKEINFO], [makeinfo], [makeinfo])
AC_CHECK_PROG([PDFTEX], [pdftex], [pdftex])
AC_CHECK_PROG([KPSEWHICH], [kpsewhich], [kpsewhich])
BUILD_DOC=false
AS_IF([test "x$TEXI2DVI" != "xtexi2dvi" || test "x$PDFTEX" != "xpdftex" || test "x$MAKEINFO" != "xmakeinfo"],
      [BUILD_DOC=false
       AC_MSG_WARN([The 'texi2dvi' or 'pdftex' or 'makeinfo' program could not be found. This isn't critical, but it means that you won't be able to create the documentation and thus to make releases.])])

AS_IF([test "x$KPSEWHICH" != "xkpsewhich" || test "$(kpsewhich  -show-path=fmt)" == "/nonesuch"],
      [BUILD_DOC=false
       AC_MSG_WARN([The 'texi2dvi' or 'pdftex' is missing a correct latex installation. This isn't critical, but it means that you won't be able to create the documentation and thus to make releases.])])
AM_CONDITIONAL([BUILD_DOC], [$BUILD_DOC])

LT_INIT
LT_PATH_LD

ISODATE=`date +%Y-%m-%d`
AC_SUBST([ISODATE])

ARCH=`uname -m | sed s/i.86/i386/ | sed s/x86_64/amd64/`
AC_SUBST([ARCH])

AC_PATH_PROG([DEB_BUILDER], [dpkg])

dnl Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h])

# Checks for library functions.
AC_C_INLINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday memset select])
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UNSIGNED_LONG_LONG_INT
dnl AC_CHECK_TYPE([uint],
dnl               [],
dnl               [AC_DEFINE_UNQUOTED([uint], [unsigned int],
dnl                [Define to `unsigned int' if <sys/types.h> does not define.])])

AC_CHECK_LIB([z], [gzread], ,
             [AC_MSG_ERROR([zlib not found, see http://www.zlib.net])])
AC_CHECK_LIB([pthread], [pthread_create], ,
             [AC_MSG_ERROR([pthread library not found])])


CPPFLAGS="$CPPFLAGS -Wall -ansi -funroll-loops"

hacked_prefix="${HACKED_PREFIX:-$prefix}"
if test "x$hacked_prefix" == "xNONE"; then
  hacked_prefix=$ac_default_prefix
fi
AC_SUBST([hacked_prefix])

dnl Additional Scripts for debian packages
AC_SUBST_FILE([postrm_script])
postrm_script=$srcdir/postrm
AC_SUBST_FILE([postinst_script])
postinst_script=$srcdir/postinst

dnl Add option --enable-debug (default to --disable-debug)
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging informations])],
               [DEBUG=1
                AC_DEFINE([OUTPUT_NB_OCC_FACTORS],[1],[Output the factors' number of occurences])],
              [DEBUG=0])

AC_ARG_ENABLE([coverage],
              [AS_HELP_STRING([--enable-coverage],
                              [enable coverage information])],
              [COVERAGE=1],
              [COVERAGE=0])

dnl Add option --without-included-ProgressBar (default --with-included-ProgressBar=check)
AC_ARG_WITH([included-ProgressBar],
            [AS_HELP_STRING([--without-included-ProgressBar],
                            [don't use the ProgressBar library included here])],
            [],
            [with_included_ProgressBar="check"])

dnl Add option --with-libProgressBar-prefix=PATH (default --without-libProgressBar-prefix)
AC_ARG_WITH([libProgressBar-prefix],
            [AS_HELP_STRING([--with-libProgressBar-prefix=DIR],
                            [search for libProgressBar in DIR/include/libProgressBar and DIR/lib])],
            [LDFLAGS="$LDFLAGS -L$with_libProgressBar_prefix/lib"
             CPPFLAGS="$CPPFLAGS -I$with_libProgressBar_prefix/include/libProgressBar"],
            [with_libProgressBar_prefix=""])

dnl Check if we need included ProgressBar library.
PB_OK=0
AS_IF([test "x$with_included_ProgressBar" == "xcheck"],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [with_included_ProgressBar="no"],
             [with_included_ProgressBar="yes"])])

AS_IF([test "x$with_included_ProgressBar" == "xyes"],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-ProgressBar and --with-libProgressBar-prefix options])],
             [AC_CONFIG_SUBDIRS([src/libProgressBar-1.2.0])
              PB_OK=1])],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [PB_OK=1],
             [AC_CHECK_LIB([ProgressBar],
                           [libProgressBarVersion],
                           [PB_OK=1
                            ProgressBar_LIBS=-lProgressBar])])]) dnl Override the default behavior of AC_CHECK_LIB (we don't want -lProgressBar in LIBS).

AS_IF([test "$PB_OK" == "1"], [AC_DEFINE([HAVE_LIBPROGRESSBAR], [1], [Define to 1 in order to Use libProgressBar])])

dnl Set INCLUDED_PROGRESSBAR conditionnal for automake
AM_CONDITIONAL([INCLUDED_PROGRESSBAR], [test "x$with_included_ProgressBar" == xyes])

if test $COVERAGE -eq 1 || test $DEBUG -eq 1; then
   CPPFLAGS="$CPPFLAGS -g -O0"
   if test $COVERAGE -eq 1; then
      CPPFLAGS="$CPPFLAGS -fprofile-arcs -ftest-coverage"
      LDFLAGS="-fprofile-arcs -lgcov"
   fi
else
   CPPFLAGS="$CPPFLAGS -O4 -DNDEBUG"
fi

LDFLAGS="$LDFLAGS $ProgressBar_LIBS"
CPPFLAGS="$CPPFLAGS $ProgressBar_CFLAGS"

dnl Output files and make debian scripts executable
AC_CONFIG_FILES([\
  postinst prerm control libGkArrays.pc \
  Makefile src/Makefile tests/Makefile \
  src/sux/Makefile src/kseq/Makefile \
  doc/Makefile doc/libGkArrays.3 doc/package.texi \
], [ \
  chmod 755 postinst && chmod 755 prerm \
])
AC_OUTPUT
