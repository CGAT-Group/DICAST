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
 
#ifndef SOL_ARRAY_H
#define SOL_ARRAY_H

#include <stdexcept>
#include "gkArraysTypes.h"
#include "bitvector.hxx"

using namespace std;

namespace gkarrays {

  typedef enum {SMALL_ARRAY, LARGE_ARRAY, OPTIMAL_ARRAY} array_type;

  /**
   * Small Or Large (SOL) Array storing integers. A SOL array can either
   * store 32-bit, 64-bit integers, or fixed length integers on a different number
   * of bits, depending on the user choice.
   * This class allows to deal with the same variable whatever the type of
   * integers we are storing.
   */
  class SolArray {
  public:
    array_type type;
    uintSA *large;
    uint *small;
    BitVector<uintSA, uintSA> *bv;
    uintSA len;

    /**
     * Default constructor
     * @post getType() == SMALL_ARRAY && length() == 0
     */
    SolArray();

    /**
     * @param nbElements: number of elements to be stored in the array.
     * @param t: type of the array to be constructed.
     *           This specifies if we need to build a small array (32-bit integers)
     *           a large array (64-bit integers) or an optimal array which uses the
     *           optimal number of bits in memory (but which is longer)
     */
    SolArray(uintSA nbElements, array_type t);

    /**
     * @param nbElements: number of elements to be stored in the array.
     * @param max: maximal value stored in the array
     * @param t: type of the array to be constructed.
     *           This specifies if we need to build a small array (32-bit integers)
     *           a large array (64-bit integers) or an optimal array which uses the
     *           optimal number of bits in memory (but which is longer)
     */
    SolArray(uintSA nbElements, uintSA max, array_type t);

    /**
     * @param ptr: pointer to an already allocated memory
     * @param t: type of the array to be built
     * @param nbElements: Number of elements that can be stored in the array
     * @pre t == SMALL_ARRAY || t == LARGE_ARRAY (doestn't work for bit vectors -- yet?)
     */
    SolArray(void *ptr, array_type t, uintSA nbElements) throw (invalid_argument);

    /**
     * Copy constructor.
     * Warning, the arrays are not copied, only the pointers are (to avoid
     * huge memory consumption).
     */
    SolArray(SolArray &sa);

    ~SolArray();

    /**
     * @param nbElements: number of elements to be stored in the array.
     * @param max: maximal value stored in the array.
     * @param t: type of the array to be constructed
     * @pre The previous arrays must have been deleted if needed.
     * @post length() == nbElements && getType() == t
     */
    void init(uintSA nbElements, uintSA max, array_type t) ;

    /**
     * @param pos: position in the array
     * @pre pos >= 0 && pos < length()
     * @return the value in the array at position pos
     */
    uintSA get(uintSA pos);

    /**
     * @return the type of the SolArray that was constructed .
     */
    array_type getType();

    /**
     * @return The number of elements we can store in the array.
     */
    uintSA length();

    /**
     * @param nbElements: number of elements the array must be reallocated to
     * @pre getType() != OPTIMAL_ARRAY
     * @post length() == nbElements AND the memory has been reallocated.
     */
    void realloc(uintSA nbElements) throw (logic_error);
    
    /**
     * @param pos: position in the array
     * @param value: value to be stored.
     * @pre pos >= 0 && pos < length()
     * @post get(pos) == value
     */
    void set(uintSA pos, uintSA value);

    /**
     * A shortcut to the get method.
     */
    uintSA operator[](uintSA pos);
  };


  inline uintSA SolArray::get(uintSA pos) {
    switch(type) {
    case SMALL_ARRAY:
      return small[pos];
    case LARGE_ARRAY:
      return large[pos];
    default:
      return bv->GetValue(pos);
    }
  }

  inline array_type SolArray::getType() {
    return type;
  }

  inline uintSA SolArray::length() {
    return len;
  }

  inline void SolArray::set(uintSA pos, uintSA value) {
    switch(type) {
    case SMALL_ARRAY:
      small[pos] = value;
      break;
    case LARGE_ARRAY:
      large[pos] = value;
      break;
    case OPTIMAL_ARRAY:
      bv->SetValue(pos, value);
      break;
    }
  }

  inline uintSA SolArray::operator[](uintSA pos) {
    return get(pos);
  }
}
#endif

// Local Variables:
// mode:c++
// End:
