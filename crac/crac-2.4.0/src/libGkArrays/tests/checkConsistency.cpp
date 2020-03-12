/******************************************************************************
*                                                                             *
*  Copyright © 2010-2013 -- IRB/INSERM                                        *
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
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  *
*  développement  et à la reproduction du  logiciel par  l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  *
*  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  *
*  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  *
*  pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  *
*  termes.                                                                    *
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
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL-C license and that you accept its terms.           *
*                                                                             *
******************************************************************************/

#undef NDEBUG
#include "checkConsistency.h"
#include <gkArrays.h>
#include <utils.h>
#include <cstring>
#include <iostream>
#include "testing.h"

using namespace std;
using namespace gkarrays;

void checkConsistency(const char *file, uint k, uint read_length, bool use_bitvector,
                      bool stranded, uint nb_threads) {

#undef TAP_ADDITIONAL_INFOS
#define TAP_ADDITIONAL_INFOS "file = " << file << ", read_length = " << read_length <<", k = " << k <<", stranded = " << stranded << ", bit_vector = " << use_bitvector <<", nb_threads = " << nb_threads 

  RECORD_TAP_TEST(TEST_SORTED_STRING, "Check that strings are sorted in gkSA");
  RECORD_TAP_TEST(TEST_GKCFA, "Check GkCFA values");
  RECORD_TAP_TEST(TEST_GKISA, "Check GkIFA values");
  RECORD_TAP_TEST(TEST_GKCFALENGTH, "Check GkCFA length");
  RECORD_TAP_TEST(TEST_READING_READS, "Compare reading sequences from file and retrieving them from the structure");
  RECORD_TAP_TEST(TEST_READ_LENGTH, "Check read length");
  RECORD_TAP_TEST(TEST_EPOTN, "Check endPosOfTagNum");
  RECORD_TAP_TEST(TEST_SPOTN, "Check startPosOfTagNum");
  RECORD_TAP_TEST(TEST_SQPOTN, "Check startQPosOfTagNum");
  RECORD_TAP_TEST(TEST_SUPPORT_LENGTH, "Check support length");
  RECORD_TAP_TEST(TEST_TAGNUM, "Check read number");
  RECORD_TAP_TEST(TEST_TAGPOS, "Check tag position computation");
  RECORD_TAP_TEST(TEST_NBOCCS_TAGS, "Check number of occurrences among reads");
  RECORD_TAP_TEST(TEST_SUPPORT, "Check support");
  RECORD_TAP_TEST(TEST_KMER, "Check k-mer extraction from the structure");
  RECORD_TAP_TEST(TEST_FIND_MYSELF, "Check I find myself among the occurrences");
  RECORD_TAP_TEST(TEST_NBTAGS_FOUND, "Check number of occurrences among reads (with multiplicity");
  RECORD_TAP_TEST(TEST_PPOSITION, "Check isPposition is true when needed");
  RECORD_TAP_TEST(TEST_PPOS_TO_QPOS, "Check conversion of Pposition in Qposition");
  RECORD_TAP_TEST(TEST_NOT_PPOSITION, "Check isPposition is false when needed");
  RECORD_TAP_TEST(TEST_NB_READS, "Check number of reads");

  gkArrays *gk = new gkArrays((char *)file, k, use_bitvector, read_length, stranded,
                              nb_threads);
  readIterator *it = gk->getReads()->begin();
  uint nb_reads = 0;
  uint current_length;
  uintSA cumulative_length = 0;
  uintSA cumulativePpos = 0;

  // Testing GkISA, GkSA and GkCFPS
  char *previous_factor = gk->getTextFactor(gk->getGkSA(0), k);
  uintSA nb_distinct_factors = 0;
  uintSA nb_factors_current_bucket = 1;
  for (uintSA i = 1; i < gk->getGkSALength(); i++) {
    char *current_factor = gk->getTextFactor(gk->getGkSA(i), k);
    int compar = strcmp(previous_factor, current_factor);
    TAP_TEST(compar <= 0, TEST_SORTED_STRING,
             "i = " << i << ", gkSA[i] = " << gk->getGkSA(i) << ", we should have "
             << previous_factor << "<=" << current_factor);
 
    if (compar < 0) {
      TAP_TEST(gk->getGkCFA(nb_distinct_factors) == nb_factors_current_bucket,
               TEST_GKCFA,
               "i = " << i << ", gkSA[i] = " << gk->getGkSA(i) << ", we should have "
               << gk->getGkCFA(nb_distinct_factors) << "==" 
               << nb_factors_current_bucket);
      nb_distinct_factors++;
      nb_factors_current_bucket = 1;
    } else 
      nb_factors_current_bucket++;

    TAP_TEST(gk->getGkISA(gk->convertPposToQpos(gk->getGkSA(i))) == nb_distinct_factors, 
             TEST_GKISA,
             "gk->getGkISA(gk->convertPposToQpos(gk->getGkSA(i))) = " << gk->getGkISA(gk->convertPposToQpos(gk->getGkSA(i))) << ", nb_distinct_factors = " << nb_distinct_factors);
    delete [] previous_factor;
    previous_factor = current_factor;
  }
  TAP_TEST(gk->getGkCFALength() == nb_distinct_factors+1, 
           TEST_GKCFALENGTH,
           "|GkCFA| = " << gk->getGkCFALength() << ", nb_distinct_factors = " << nb_distinct_factors);

  while (!it->isFinished()) {
    char *current_read = gk->getTag(nb_reads);
    char *it_read = it->getSequence();
    current_length = gk->getTagLength(nb_reads);
    if (read_length) {
      dnaFiltration(it_read, read_length);
      TAP_TEST(strncmp(current_read, it_read, read_length) == 0,
               TEST_READING_READS,
               "Read n°" << nb_reads << " read " << it_read
               << ", in GkArrays: " << current_read);
      TAP_TEST(current_length == read_length,
               TEST_READ_LENGTH,
               "read_length (parameter) = " << read_length 
               << " current_length (given by gkA) = " 
               << current_length);
    } else {
      dnaFiltration(it_read, strlen(it_read));
      TAP_TEST(strcmp(current_read, it_read) == 0, 
               TEST_READING_READS,
               "current_read = " << current_read <<", it->getSequence() = " << it_read);
      TAP_TEST(strlen(current_read) == current_length,
               TEST_READ_LENGTH,
               "length of read is " << strlen(current_read) <<", should have been "
               << current_length << " (according to gkA)");
    }
    TAP_TEST(gk->getEndPosOfTagNum(nb_reads) == cumulative_length + current_length - 1,
             TEST_EPOTN,
             "nb_reads = " << nb_reads 
             << ", gk->getEndPosOfTagNum(nb_reads) = " 
             << gk->getEndPosOfTagNum(nb_reads)
             << "cumulative_length + current_length - 1 = " 
             << cumulative_length + current_length - 1);
    TAP_TEST(gk->getStartPosOfTagNum(nb_reads) == cumulative_length,
             TEST_SPOTN,
             "nb_reads = " << nb_reads 
             << ", gk->getStartPosOfTagNum(nb_reads) = " 
             << gk->getStartPosOfTagNum(nb_reads)
             << "cumulative_length = " << cumulative_length);
    TAP_TEST(gk->getStartQPosOfTagNum(nb_reads) == cumulativePpos,
             TEST_SQPOTN,
             "nb_reads = " << nb_reads 
             << ", gk->getStartQPosOfTagNum(nb_reads) = " 
             << gk->getStartQPosOfTagNum(nb_reads)
             << "cumulativePpos = " << cumulativePpos);
             
    TAP_TEST(gk->getSupportLength(nb_reads) == current_length - k + 1,
             TEST_SUPPORT_LENGTH,
             "nb_reads = " << nb_reads << " support length: "
             << gk->getSupportLength(nb_reads) << ", current_length - k + 1 = "
             << current_length - k + 1);

    for (uint i = 0; i < current_length; i++) {
      for (uint j = 1; j < current_length - i; j++) {
        char *factor = gk->getTagFactor(nb_reads, i, j);
        assert(strncmp(factor, &current_read[i], j) == 0);
        delete [] factor;
      }
    }

    for (uint i = 0; i < current_length; i++) {
      std::pair<uint, uint> tagNumPos = gk->getTagNumAndPosFromAbsolutePos(cumulative_length + i);
      TAP_TEST(gk->getTagNum(cumulative_length + i) == nb_reads, 
               TEST_TAGNUM,
               "cumulative_length = " << cumulative_length << ", i = " << i <<", nb_reads = " << nb_reads <<", gk->getTagNum(cumulative_length + i) = " << gk->getTagNum(cumulative_length + i));
      TAP_TEST(i == tagNumPos.second,
               TEST_TAGPOS,
               "i = " << i <<", tagNumPos.second = " << tagNumPos.second);
      TAP_TEST(nb_reads == tagNumPos.first, TEST_TAGNUM,
               "nb_reads = " << nb_reads <<", tagNumPos.first = " << tagNumPos.first);
    }

    uint *support = gk->getSupport(nb_reads);

    for (uint i = 0; i < current_length - k + 1; i++) {
      uint nb_tags = gk->getNbTagsWithFactor(nb_reads, i, false);
      char *kmer = gk->getTagFactor(nb_reads, i, k);
      uint *tags_num = gk->getTagNumWithFactor(nb_reads, i);
      uint nb_occs = gk->getNbTagsWithFactor(nb_reads, i, true);
      pair<uint, uint> *tags = gk->getTagsWithFactor(nb_reads, i);
      bool found_me = false;
      uint nb_found = 0;

      TAP_TEST(nb_occs >= nb_tags, TEST_NBOCCS_TAGS,
               "i = " << i << "nb_occs = " << nb_occs <<", nb_tags = " << nb_tags);
      TAP_TEST(support[i] == nb_occs, TEST_SUPPORT,
               "i = " << i << ", support[i] = " << support[i] 
               << ", nb_occs = " << nb_occs);

      for (uint j = 0; j < nb_occs; j++) {
        if (tags[j].first == nb_reads)
          found_me = true;

        // Search this tag in tags_num
        uint p;
        for (p = 0; p < nb_tags && tags_num[p] != tags[j].first ; p++) {}
        if (p < nb_tags) {
          // Count the read found
          nb_found++;
          // And change it in the array so that we don't count it twice.
          tags_num[p] = ~0;
        }

        char *current_kmer = gk->getTagFactor(tags[j].first, tags[j].second, k);
        TAP_TEST(strcmp(current_kmer, kmer) == 0, 
                 TEST_KMER,
                 "current_kmer = " << current_kmer << ", kmer = " << kmer << ", j = " << j << ", nb_occs = " << nb_occs);
        delete [] current_kmer;
      }

      TAP_TEST(found_me, TEST_FIND_MYSELF,
               "i = " << i);
      TAP_TEST(nb_found == nb_tags, 
               TEST_NBTAGS_FOUND,
               "nb_found = " << nb_found <<", nb_tags = " << nb_tags );
      TAP_TEST(gk->isPposition(cumulative_length + i),
               TEST_PPOSITION,
               "i = " << i <<", cumulative_length = " << cumulative_length);
      TAP_TEST(gk->convertPposToQpos(cumulative_length + i) == cumulativePpos + i,
               TEST_PPOS_TO_QPOS,
               "i = " << i <<", cumulative_length = " << cumulative_length
               <<", cumulativePpos = " << cumulativePpos <<
               "gk->convertPposToQpos(cumulative_length + i) = "
               << gk->convertPposToQpos(cumulative_length + i));
      delete [] kmer;
      free(tags_num);
      delete [] tags;
    }

    delete [] support;
    
    for (uint i = current_length - k + 1; i < current_length; i++) {
      TAP_TEST(! gk->isPposition(cumulative_length+i), TEST_NOT_PPOSITION,
               "cumulative_length = " << cumulative_length << ", i = " << i);
    }

    cumulative_length += current_length;
    cumulativePpos += current_length - k + 1;
    nb_reads++;
    ++(*it);

    delete [] current_read;
  }

  delete it;

  TAP_TEST(nb_reads == gk->getNbTags(),
           TEST_NB_READS,
           "nb_reads = " << nb_reads <<", gk->getNbTags() = "
           << gk->getNbTags());

  delete gk;
}
