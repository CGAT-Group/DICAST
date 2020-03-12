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
#include <iostream>
#include <cstdlib>
#include <ctime>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include "compareVariableFixed.h"
#include <gkArrays.h>
#include <gkArraysTypes.h>
#include <utils.h>
#include <sys/time.h>
#include "testing.h"

#ifndef ullong
#define ullong unsigned long long
#endif

using namespace std;
using namespace gkarrays;

void compareVariableFixed(const char *file, uint read_length, uint k, bool stranded,
                          bool bit_vector, uint nb_threads) {

  RECORD_TAP_TEST(TEST_NBTAGS, "Check number of tags");
  RECORD_TAP_TEST(TEST_TAG_NUM, "Check that tag numbers are consistent");
  RECORD_TAP_TEST(TEST_TAG_LENGTH, "Check tag lengths");
  RECORD_TAP_TEST(TEST_GKSALENGTH, "Check GkSA length");
  RECORD_TAP_TEST(TEST_GKSA, "Check GkSA values");
  RECORD_TAP_TEST(TEST_GKISA, "Check GkISA values");
  RECORD_TAP_TEST(TEST_GKCFALENGTH, "Check GkCFA length");
  RECORD_TAP_TEST(TEST_GKCFA, "Check GkCFA values");

#undef TAP_ADDITIONAL_INFOS
#define TAP_ADDITIONAL_INFOS "file = " << file << ", read_length = " << read_length <<", k = " << k <<", stranded = " << stranded << ", bit_vector = " << bit_vector <<", nb_threads = " << nb_threads 

  gkArrays *gkFixed = new gkArrays((char *)file, k, bit_vector, read_length, stranded,
                                   nb_threads);
  gkArrays *gkVariable = new gkArrays((char *)file, k, bit_vector, 0, stranded, nb_threads);

  TAP_TEST(gkFixed->getNbTags() == gkVariable->getNbTags(),
           TEST_NBTAGS,
           "Nb fixed tags (" << gkFixed->getNbTags() << ") != Nb variable tags ("
           << gkVariable->getNbTags() << ")");

  uint nbtags = gkFixed->getNbTags();

  for (uintSA i = 0; i < nbtags; i++) {
    for (uintSA j = 0; j < (uintSA) read_length; j++)
      TAP_TEST(i == gkFixed->getTagNum(i * read_length + j) && i == gkVariable->getTagNum(i * read_length + j ), 
               TEST_TAG_NUM,
               "i = " << i << ", j = " << j << ", num_fixed = " << gkFixed->getTagNum(i * read_length + j) << ", num_var = " << gkVariable->getTagNum(i * read_length + j));

    TAP_TEST(gkFixed->getTagLength(i) == (uintSA)read_length && gkVariable->getTagLength(i) == (uintSA)read_length, 
             TEST_TAG_LENGTH,
             "i = " << i << " fixed tag length: " << gkFixed->getTagLength(i) <<", variable tag length: " << gkVariable->getTagLength(i));
  }


  TAP_TEST(gkFixed->getGkSALength() == gkVariable->getGkSALength(), 
           TEST_GKSALENGTH,
           "fixed gkSA length = " << gkFixed->getGkSALength() << ", variable gkSA length = " << gkVariable->getGkSALength());

  uintSA max = gkFixed->getGkSALength();

  for (uintSA i = 0; i < max; i++) {
    TAP_TEST(gkFixed->getGkSA(i) == gkVariable->getGkSA(i), 
             TEST_GKSA,
             "i = " << i);
    TAP_TEST(gkFixed->getGkISA(i) == gkVariable->getGkISA(i), 
             TEST_GKISA,
             "i = " << i << ", GkISA[i] (fixed) = " << gkFixed->getGkISA(i)
             << ", GkISA[i] (variable) = " << gkVariable->getGkISA(i));
  }
  
  TAP_TEST(gkFixed->getGkCFALength() == gkVariable->getGkCFALength(), 
           TEST_GKCFALENGTH,
           "fixed GkCFA length = " << gkFixed->getGkCFALength() << ", variable GkCFA length = " << gkVariable->getGkCFALength());

  max = gkFixed->getGkCFALength();
  
  for (uintSA i = 0; i < max; i++) {
    TAP_TEST(gkFixed->getGkCFA(i) == gkVariable->getGkCFA(i), 
             TEST_GKCFA,
             "i = " << i);
  }
  
  delete gkFixed;
  delete gkVariable;
}
