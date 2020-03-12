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


bin_PROGRAMS = crac-index
crac_index_SOURCES = cracIndex.cpp cracIndex.h

AM_CPPFLAGS = -I@abs_top_srcdir@ -I@abs_srcdir@/karkkainen_bwt -I@abs_top_srcdir@/src/libGkArrays/src
AM_LDFLAGS = -lm -lpthread -lSSA -L@abs_builddir@/
LDADD = libSSA.a

KARKKAINEN_BWT =                     \
  karkkainen_bwt/blockwise_sa.h      \
  karkkainen_bwt/zbox.h              \
  karkkainen_bwt/assert_helpers.h    \
  karkkainen_bwt/random_source.h     \
  karkkainen_bwt/binary_sa_search.h  \
  karkkainen_bwt/multikey_qsort.h    \
  karkkainen_bwt/auto_array.h        \
  karkkainen_bwt/diff_sample.h       \
  karkkainen_bwt/alphabet.h          \
  karkkainen_bwt/timer.h

SEQAN_INDEX =                            \
  seqan/index.h                          \
  seqan/index/index_qgram_find.h         \
  seqan/index/index_sa_lss.h             \
  seqan/index/index_generated_forwards.h \
  seqan/index/shape_base.h               \
  seqan/index/shape_predefined.h         \
  seqan/index/index_manual_forwards.h    \
  seqan/index/shape_gapped.h             \
  seqan/index/index_find.h               \
  seqan/index/index_base.h

SEQAN_SEQUENCE =                               \
  seqan/sequence.h                             \
  seqan/sequence/string_pointer.h              \
  seqan/sequence/sequence_forwards.h           \
  seqan/sequence/string_alloc.h                \
  seqan/sequence/string_base.h                 \
  seqan/sequence/sequence_shortcuts.h          \
  seqan/sequence/sequence_generated_forwards.h \
  seqan/sequence/sequence_multiple.h           \
  seqan/sequence/string_packed.h               \
  seqan/sequence/lexical.h                     \
  seqan/sequence/string_value_expand.h         \
  seqan/sequence/std_string.h                  \
  seqan/sequence/string_array.h                \
  seqan/sequence/segment_prefix.h              \
  seqan/sequence/segment_base.h                \
  seqan/sequence/segment_infix.h               \
  seqan/sequence/segment_suffix.h              \
  seqan/sequence/sequence_interface.h          \
  seqan/sequence/string_cstyle.h               \
  seqan/sequence/string_stack.h

SEQAN_FILE =                           \
  seqan/file.h                         \
  seqan/file/file_format_fasta_align.h \
  seqan/file/file_format_raw.h         \
  seqan/file/file_format.h             \
  seqan/file/file_format_guess.h       \
  seqan/file/cstream.h                 \
  seqan/file/file_base.h               \
  seqan/file/file_format_fasta.h       \
  seqan/file/file_generated_forwards.h \
  seqan/file/file_forwards.h           \
  seqan/file/file_array.h              \
  seqan/file/chunk_collector.h         \
  seqan/file/stream.h                  \
  seqan/file/file_format_cgviz.h       \
  seqan/file/stream_algorithms.h       \
  seqan/file/meta.h

SEQAN_FIND =                           \
  seqan/find.h                         \
  seqan/find/find_score.h              \
  seqan/find/find_set_horspool.h       \
  seqan/find/find_generated_forwards.h \
  seqan/find/find_base.h               \
  seqan/find/find_multi.h              \
  seqan/find/find_horspool.h           \
  seqan/find/find_multiple_shiftand.h

SEQAN_BASIC =                              \
  seqan/basic.h                            \
  seqan/basic/basic_operator.h             \
  seqan/basic/basic_allocator_multipool.h  \
  seqan/basic/basic_compare.h              \
  seqan/basic/basic_allocator_singlepool.h \
  seqan/basic/basic_allocator_to_std.h     \
  seqan/basic/basic_alphabet_interface2.h  \
  seqan/basic/basic_pointer.h              \
  seqan/basic/basic_alphabet_simple_tabs.h \
  seqan/basic/basic_iterator.h             \
  seqan/basic/basic_converter.h            \
  seqan/basic/basic_alphabet_interface.h   \
  seqan/basic/basic_metaprogramming.h      \
  seqan/basic/basic_profile.h              \
  seqan/basic/basic_proxy.h                \
  seqan/basic/basic_alphabet_trait_basic.h \
  seqan/basic/basic_allocator_simple.h     \
  seqan/basic/basic_alphabet_simple.h      \
  seqan/basic/basic_host.h                 \
  seqan/basic/basic_counted_ptr.h          \
  seqan/basic/basic_allocator_interface.h  \
  seqan/basic/basic_iterator_position.h    \
  seqan/basic/basic_iterator_adaptor.h     \
  seqan/basic/basic_iterator_adapt_std.h   \
  seqan/basic/basic_iterator_simple.h      \
  seqan/basic/basic_forwards.h             \
  seqan/basic/basic_debug.h                \
  seqan/basic/basic_generated_forwards.h   \
  seqan/basic/basic_transport.h            \
  seqan/basic/basic_type.h                 \
  seqan/basic/basic_holder.h               \
  seqan/basic/basic_aggregates.h           \
  seqan/basic/basic_definition.h           \
  seqan/basic/basic_volatile_ptr.h         \
  seqan/basic/basic_iterator_base.h        \
  seqan/basic/basic_tag.h

SEQAN =             \
  seqan/platform.h  \
  $(SEQAN_INDEX)    \
  $(SEQAN_SEQUENCE) \
  $(SEQAN_FILE)     \
  $(SEQAN_FIND)     \
  $(SEQAN_BASIC)

noinst_LIBRARIES = libSSA.a
libSSA_a_SOURCES =                        \
  basic.h	                          \
  bitarray.cpp bitarray.h                 \
  bitrankw32int.cpp bitrankw32int.h       \
  chrPosition.cpp chrPosition.h           \
  HuffAlphabetRank.cpp HuffAlphabetRank.h \
  Huffman_Codes.cpp Huffman_Codes.h       \
  locateOnGenome.cpp locateOnGenome.h     \
  locateTags.cpp locateTags.h     \
  SSA.cpp SSA.h                           \
  utils.cpp utils.h                       \
  GenomeInfo.cpp GenomeInfo.h             \
  $(SEQAN) $(KARKKAINEN_BWT) $(KSEQ)

MOSTLYCLEANFILES=*~
CLEANFILES=*~ $(EXTRA_PROGRAMS)
DISTCLEANFILES=*~
MAINTAINERCLEANFILES=*~
