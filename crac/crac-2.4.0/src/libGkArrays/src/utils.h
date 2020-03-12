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
*                   Jérôme AUDOUX    <jerome.audoux@etud.univ-montp2.fr>      *
*  with additional contribution for the packaging of:	                      *
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

/*
 * =============================================
 *
 * $Id: utils.h 397 2011-04-14 08:47:49Z nphilippe $
 *
 * =============================================
 */
 
#ifndef UTILS_H
#define UTILS_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <gkArraysTypes.h>
#include <cstdlib>
#include <stdexcept>
#include <cassert>

#ifndef NDEBUG
#define PRINT_VAR(v) std::cerr << #v << " = " << v << std::endl;
#define ASSERT_MSG(condition, message) \
  if (!(condition)) { std::cerr << message << std::endl; assert(condition); }
#else
#define ASSERT_MSG(condition, message) ((void)0)
#define PRINT_VAR(v) ((void)0)

#endif

namespace gkarrays {

  template <typename T>
  inline T *Malloc() {
    T *ptr;
    if (!(ptr = (T *) malloc(sizeof(T)))) {
      throw std::bad_alloc();
    }
    return ptr;
  }

  template <typename T>
  inline T *Calloc(size_t nb) {
    T *ptr;
    if (!(ptr = (T *) calloc(nb, sizeof(T)))) {
      throw std::bad_alloc();
    }
    return ptr;
  }

  template <typename T>
  inline T *Realloc(T *ptr, size_t nb) {
    if (!(ptr = (T *) realloc((void *) ptr, nb * sizeof(T)))) {
      throw std::bad_alloc();
    }
    return ptr;
  }

  int comparUint(const void *a1, const void* a2);

  int comparUintSA(const void *a1, const void* a2);

  uint convNuc(char nuc);

  char intToNuc(uint c);

  /**
   * Modify the parameter so that every lowercase nucleotide
   * is changed to an uppercase nucleotide
   * and every unknown letter (other than A, C, G or T) is changed to A.
   * The second parameter corresponds to the length of the DNA string
   * to be treated.
   */
  void dnaFiltration(char *, uintSA);

  uintSA DNAtoInt(char *dna, uintSA dna_length);

  void intToDNA(uint64 code, uint dna_length, char *dna);
  
  /**
   * @param pos: position in the string
   * @param dna: DNA string, where each nucleotide is encoded on 2 bits
   * @param length: length of the factor to extract
   * @return a numeric representation of the DNA factor
   * @pre length >= 4 && length <= 32
   */
  /**
   * @param pos: position in the string
   * @param dna: DNA string, where each nucleotide is encoded on 2 bits
   * @param length: length of the factor to extract
   * @return a numeric representation of the DNA factor
   * @pre length >= 4 && length <= 32
   */
  uint64 factorsToInt(uintSA pos, uchar *dna, uint length);

  /**
   * @param pos: position in the string
   * @param dna: DNA string, where each nucleotide is encoded on 2 bits
   * @param length: length of the factor to extract
   * @return a numeric representation of the DNA factor such that
   *         a factor and its reverse complement share the same value.
   * @pre length >= 4 && length <= 32
   */
  uint64 factorsToIntNoStranded(uintSA pos, uchar *dna, uint length);

  /**
   * @param factor: a numeric representation of a DNA factor
   * @param length: length of the DNA factor
   * @return the numeric representation of the DNA factor's revcomp.
   */
  uint64 intRevcomp(uint64 factor,  uint length);

  /**
   * @return the integer corresponding to the complementary nucleotide
   * @post intToNuc(complementaryNuc(convNuc(c))) == c'
   *       where c is a nucleotide and c' its complementary
   */
  inline uintSA complementaryNuc(int i){
    return i ^ 3;
  }

}
#endif

// Local Variables:
// mode:c++
// End:
