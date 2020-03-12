
dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2016 -- LIRMM/CNRS                                      
dnl                          (Laboratoire d'Informatique, de Robotique et de 
dnl                          Microélectronique de Montpellier /              
dnl                          Centre National de la Recherche Scientifique)   
dnl                          LIFL/INRIA                                      
dnl                          (Laboratoire d'Informatique Fondamentale de     
dnl                          Lille / Institut National de Recherche en       
dnl                          Informatique et Automatique)                    
dnl                          LITIS                                           
dnl                          (Laboratoire d'Informatique, du Traitement de   
dnl                          l'Information et des Systèmes).                 
dnl                                                                          
dnl Copyright © 2011-2016 -- IRB/INSERM                                      
dnl                          (Institut de Recherches en Biothérapie /        
dnl                          Institut National de la Santé et de la Recherche
dnl                          Médicale).                                      
dnl                                                                          
dnl Copyright © 2015-2016 -- AxLR/SATT                                       
dnl                          (Lanquedoc Roussilon /                          
dnl                           Societe d'Acceleration de Transfert de         
dnl                           Technologie).	                                
dnl                                                                          
dnl Programmeurs/Progammers:                                                 
dnl                   Nicolas PHILIPPE <nphilippe.resear@gmail.com>          
dnl                   Mikaël SALSON    <mikael.salson@lifl.fr>               
dnl                   Jérôme Audoux    <jerome.audoux@gmail.com>             
dnl  with additional contribution for the packaging of:	                    
dnl                   Alban MANCHERON  <alban.mancheron@lirmm.fr>            
dnl                                                                          
dnl  Contact:         CRAC list   <crac-bugs@lists.gforge.inria.fr>          
dnl  Paper:           CRAC: An integrated RNA-Seq read analysis              
dnl                   Philippe N., Salson M., Commes T., Rivals E.           
dnl                   Genome Biology 2013; 14:R30.                           
dnl                                                                          
dnl -------------------------------------------------------------------------
dnl                                                                          
dnl  This File is part of the CRAC program.                                  
dnl                                                                          
dnl  This program is free software: you can redistribute it and/or modify    
dnl  it under the terms of the GNU General Public License as published by    
dnl  the Free Software Foundation, either version 3 of the License, or (at   
dnl  your option) any later version.  This program is distributed in the     
dnl  hope that it will be useful, but WITHOUT ANY WARRANTY; without even     
dnl  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR     
dnl  PURPOSE.  See the GNU General Public License for more details.  You     
dnl  should have received a copy of the GNU General Public License along     
dnl  with this program.  If not, see <http://www.gnu.org/licenses/>.         

dnl $Id$

dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
m4_define([V_MAJOR], [2])
m4_define([V_MINOR], [4])
m4_define([V_MICRO], [0])

AC_INIT(
  [crac],
  [V_MAJOR.V_MINOR.V_MICRO],
  [crac-bugs@lists.gforge.inria.fr],
  [crac],
  [http://crac.gforge.inria.fr/])

PACKAGE_DESCRIPTION="Crac is a tool to analyse RNA-Seq data provided by NGS"
AC_SUBST([PACKAGE_DESCRIPTION])

AC_CONFIG_SRCDIR([src/crac.cpp])
AC_REVISION([$Revision$])
AC_COPYRIGHT([
/******************************************************************************
*  Copyright © 2009-2016 -- LIRMM/CNRS                                        *
*                           (Laboratoire d'Informatique, de Robotique et de   *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique)     *
*                           LIFL/INRIA                                        *
*                           (Laboratoire d'Informatique Fondamentale de       *
*                           Lille / Institut National de Recherche en         *
*                           Informatique et Automatique)                      *
*                           LITIS                                             *
*                           (Laboratoire d'Informatique, du Traitement de     *
*                           l'Information et des Systèmes).                   *
*                                                                             *
*  Copyright © 2011-2016 -- IRB/INSERM                                        *
*                           (Institut de Recherches en Biothérapie /          *
*                           Institut National de la Santé et de la Recherche  *
*                           Médicale).                                        *
*                                                                             *
*  Copyright © 2015-2016 -- AxLR/SATT                                         *
*                           (Lanquedoc Roussilon /                            *
*                            Societe d'Acceleration de Transfert de           *
*                            Technologie).	                                  *
*                                                                             *
*  Programmeurs/Progammers:                                                   *
*                    Nicolas PHILIPPE <nphilippe.resear@gmail.com>            *
*                    Mikaël SALSON    <mikael.salson@lifl.fr>                 *
*                    Jérôme Audoux    <jerome.audoux@gmail.com>               *
*   with additional contribution for the packaging of:	                      *
*                    Alban MANCHERON  <alban.mancheron@lirmm.fr>              *
*                                                                             *
*   Contact:         CRAC list   <crac-bugs@lists.gforge.inria.fr>            *
*   Paper:           CRAC: An integrated RNA-Seq read analysis                *
*                    Philippe N., Salson M., Commes T., Rivals E.             *
*                    Genome Biology 2013; 14:R30.                             *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*   This File is part of the CRAC program.                                    *
*                                                                             *
*   This program is free software: you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation, either version 3 of the License, or (at     *
*   your option) any later version.  This program is distributed in the       *
*   hope that it will be useful, but WITHOUT ANY WARRANTY; without even       *
*   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR       *
*   PURPOSE.  See the GNU General Public License for more details.  You       *
*   should have received a copy of the GNU General Public License along       *
*   with this program.  If not, see <http://www.gnu.org/licenses/>.           *
*                                                                             *
******************************************************************************/
])

AS_IF([test x"$CFLAGS" == "x"], [CFLAGS=""])
AS_IF([test x"$CPPFLAGS" == "x"], [CPPFLAGS=""])
AS_IF([test x"$CXXFLAGS" == "x"], [CXXFLAGS=""])
AS_IF([test x"$LDFLAGS" == "x"], [LDFLAGS=""])

AC_CONFIG_AUX_DIR([config])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config/config.h])

dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE()

AC_CANONICAL_HOST
AC_LANG([C++])

AC_CONFIG_MACRO_DIR([config/m4])

dnl Checks for programs.
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_LN_S

AC_SUBST([VERSION])

dnl Disable shared library construction by default (must appear before LT_INIT)
AC_DISABLE_SHARED

LT_INIT
LT_PATH_LD

ISODATE=`date +%Y-%m-%d`
AC_SUBST([ISODATE])

ARCH=`uname -m | sed s/i.86/i386/ | sed s/x86_64/amd64/`
AC_SUBST([ARCH])

AC_PATH_PROG([DEB_BUILDER], [dpkg])
AC_PATH_PROG([RPM_BUILDER], [rpmbuild])

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [pthread_join])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_MAJOR
AC_CHECK_HEADERS([\
  fcntl.h inttypes.h libintl.h limits.h locale.h \
  memory.h stdint.h stdlib.h string.h strings.h \
  sys/file.h sys/param.h sys/time.h unistd.h utime.h \
])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

AX_CXX_COMPILE_STDCXX_11()

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([\
  floor gettimeofday memchr memmove memset \
  mkfifo pow select sqrt strrchr strstr strtol \
])


CPPFLAGS="$CPPFLAGS -Wall -ansi -funroll-loops"

dnl Add option --enable-debug (default to --disable-debug)
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging informations])],
              [CPPFLAGS="$CPPFLAGS -g -O0"],
              [CPPFLAGS="$CPPFLAGS -O4 -DNDEBUG"])

dnl Add option --without-included-GkArrays (default --with-included-GkArrays=check)
AC_ARG_WITH([included-GkArrays],
            [AS_HELP_STRING([--without-included-GkArrays],
                            [don't use the GkArrays library included here])],
            [],
            [with_included_GkArrays="check"])

dnl Add option --with-libGkArrays-prefix=PATH (default --without-libGkArrays-prefix)
AC_ARG_WITH([libGkArrays-prefix],
            [AS_HELP_STRING([--with-libGkArrays-prefix=DIR],
                            [search for libGkArrays in DIR/include/libGkArrays and DIR/lib])],
            [LDFLAGS="$LDFLAGS -L$with_libGkArrays_prefix/lib"
             CPPFLAGS="$CPPFLAGS -I$with_libGkArrays_prefix/include/libGkArrays"],
            [with_libGkArrays_prefix=""])

dnl Check if we need included GkArrays library.
AS_IF([test "x$with_included_GkArrays" == "xcheck"],
      [AS_IF([test "x$with_libGkArrays_prefix" != "x"],
             [with_included_GkArrays="no"],
             [with_included_GkArrays="yes"])])

AS_IF([test "x$with_included_GkArrays" == "xyes"],
      [AS_IF([test "x$with_libGkArrays_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-GkArrays and --with-libGkArrays-prefix options])])])

dnl Add option --without-included-jellyfish (default --with-included-jellyfish=check)
AC_ARG_WITH([included-jellyfish],
            [AS_HELP_STRING([--without-included-jellyfish],
                            [don't use the jellyfish library included here])],
            [],
            [with_included_jellyfish="check"])

dnl Add option --with-libjellyfish-prefix=PATH (default --without-libjellyfish-prefix)
AC_ARG_WITH([libjellyfish-prefix],
            [AS_HELP_STRING([--with-libjellyfish-prefix=DIR],
                            [search for libjellyfish in DIR/include/jellyfish and DIR/lib])],
            [LDFLAGS="$LDFLAGS -L$with_libjellyfish_prefix/lib"
             CPPFLAGS="$CPPFLAGS -I$with_libjellyfish_prefix/include/jellyfish"],
            [with_libjellyfish_prefix=""])

dnl Check if we need included jellyfish library.
AS_IF([test "x$with_included_jellyfish" == "xcheck"],
      [AS_IF([test "x$with_libjellyfish_prefix" != "x"],
             [with_included_jellyfish="no"],
             [with_included_jellyfish="yes"])])

AS_IF([test "x$with_included_jellyfish" == "xyes"],
      [AS_IF([test "x$with_libjellyfish_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-jellyfish and --with-libjellyfish-prefix options])])])

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
             [with_included_ProgressBar="$with_included_GkArrays"])])

AS_IF([test "x$with_included_ProgressBar" == "xyes"],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-ProgressBar and --with-libProgressBar-prefix options])],
             [AS_IF([test "x$with_included_GkArrays" != "xyes"],
                    [AC_MSG_FAILURE([You can't use both --with-included-ProgressBar and --without-included-GkArrays options])],
                    [PB_OK=1])
             ])],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [PB_OK=1
              LDFLAGS="$LDFLAGS -L$with_libProgressBar_prefix -L$with_libProgressBar_prefix/lib"
              CPPFLAGS="$CPPFLAGS -I$with_libProgressBar_prefix -I$with_libProgressBar_prefix/include"],
             [AC_CHECK_LIB([ProgressBar],
                           [libProgressBarVersion],
                           [PB_OK=1
                            ProgressBar_LIBS=-lProgressBar])])]) dnl Override the default behavior of AC_CHECK_LIB (we don't want -lProgressBar in LIBS).

AS_IF([test "$PB_OK" == "1"], [AC_DEFINE([HAVE_LIBPROGRESSBAR], [1], [Define to 1 in order to Use libProgressBar])])

dnl Set INCLUDED_PROGRESSBAR conditionnal for automake
AM_CONDITIONAL([INCLUDED_PROGRESSBAR], [test "x$with_included_ProgressBar" == "xyes"])

AC_CHECK_LIB([z], [gzread], ,
             [AC_MSG_ERROR([zlib not found, see http://www.zlib.net])])

AC_CHECK_LIB([hts], [hts_version], ,
             [AC_MSG_ERROR([htslib not found, see http://www.htslib.org/])])

GK_OK=1
AS_IF([test "x$with_included_GkArrays" != "xyes"],
      [AC_CHECK_LIB([GkArrays],
                    [libGkArraysVersion],
                    [:],dnl Override the default behavior of AC_CHECK_LIB (we don't want -lGkArrays in LIBS).
                    [GK_OK=0])],
      [AC_CONFIG_SUBDIRS([src/libGkArrays])])

HAS_READ_INDEX=0

AS_IF([test "$GK_OK" == "1"], 
      [AC_DEFINE([HAVE_LIBGKARRAYS],
                 [1],
                 [Define to 1 in order to Use libGkArrays])
       HAS_READ_INDEX=1],
      [AC_MSG_WARN([Library GkArrays not found])])

JF_OK=1
AS_IF([test "x$with_included_jellyfish" != "xyes"],
      [AC_CHECK_LIB([jellyfish-2.0],
                    [libJellyfishVersion],
                    [:],dnl Override the default behavior of AC_CHECK_LIB (we don't want -ljellyfish-2.0 in LIBS).
                    [JF_OK=0])],
      [AC_CONFIG_SUBDIRS([src/libJellyfish-2.0])])

AS_IF([test "$JF_OK" == "1"], 
      [AC_DEFINE([HAVE_LIBJELLYFISH],
                 [1],
                 [Define to 1 in order to Use libjellyfish-2.0])
       HAS_READ_INDEX=1],
      [AC_MSG_WARN([Library jellyfish-2.0 not found])])

dnl Check if we have at least on read index
AS_IF([test "$HAS_READ_INDEX" == "0"], [AC_MSG_FAILURE([No read index library found, please use included libraries])])


dnl Set INCLUDED_GKARRAYS conditionnal for automake
AM_CONDITIONAL([INCLUDED_GKARRAYS], [test "x$with_included_GkArrays" == xyes])
AM_CONDITIONAL([HAVE_LIBGKARRAYS], [test "$GK_OK" == "1"])

dnl Set INCLUDED_JELLYFISH conditionnal for automake
AM_CONDITIONAL([INCLUDED_JELLYFISH], [test "x$with_included_jellyfish" == xyes])
AM_CONDITIONAL([HAVE_LIBJELLYFISH], [test "$JF_OK" == "1"])

dnl Forward options given to current configure script to libGkArrays configure script if needed
AS_IF([test "x$enable_shared" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-shared"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-shared"])])
AS_IF([test "x$enable_static" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-static"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-static"])])

dnl Forward FLAGS to libGkArrays/libProgressBar configure scripts if needed
AS_VAR_APPEND([ac_configure_args], [" CPPFLAGS=\"$CPPFLAGS\""])
AS_VAR_APPEND([ac_configure_args], [" CFLAGS=\"$CFLAGS\""])
AS_VAR_APPEND([ac_configure_args], [" CXXFLAGS=\"$CXXFLAGS\""])

LDFLAGS="$LDFLAGS $ProgressBar_LIBS"
CPPFLAGS="$CPPFLAGS $ProgressBar_CFLAGS"

AC_CONFIG_FILES([\
  control crac.spec                             \
  Makefile src/Makefile scripts/Makefile        \
  src/libReadsInfo/Makefile src/libSSA/Makefile \
  src/libgzstream/Makefile                      \
  src/libReadsReader/Makefile                   \
  src/libReadsReader/kseq/Makefile              \
  man/Makefile man/crac.1                       \
  scripts/crac-client                           \
  test_f/Makefile test_f/bug/Makefile           \
  test_f/causes/Makefile test_f/cigar/Makefile  \
  test_f/splices_adjust/Makefile
], [ \
  test -f scripts/crac-client && chmod 755 scripts/crac-client ; \
])

dnl Configure jellyfish with its own ./configure file
dnl AC_CONFIG_SUBDIRS([src/libJellyfish/Makefile])

AC_OUTPUT
