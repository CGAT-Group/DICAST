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
*      /Progammers:                                                           *
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

#ifndef __LARGEBINVECT_HXX__
#define __LARGEBINVECT_HXX__

#include "binvect.hxx"

template <typename CountWordSize = unsigned char>
class LargeBinVect {
private:
  unsigned char *data;
  CountWordSize length;

  // Get the data internal array length
  CountWordSize GetDataArrayLength() const;
  // Get the last unsigned char mask for the internal array
  unsigned char GetDataArrayEndingMask() const;

public:
  // Constructors
  LargeBinVect(const CountWordSize length = 0, const bool value = 0) throw (std::bad_alloc);
  LargeBinVect(const char *binary_string) throw (BinVectOverflow, std::bad_alloc);
  LargeBinVect(const LargeBinVect<CountWordSize> &v) throw (std::bad_alloc);

  // Destructors
  ~LargeBinVect();

  // Operators
  inline LargeBinVect<CountWordSize> &operator=(const LargeBinVect<CountWordSize> &v) throw (std::bad_alloc);

  // Bitwise NOT
  inline LargeBinVect<CountWordSize> operator~() const;
  // Bitwise XOR
  inline LargeBinVect<CountWordSize> operator^(const LargeBinVect<CountWordSize> &v) const;
  // Bitwise OR
  inline LargeBinVect<CountWordSize> operator|(const LargeBinVect<CountWordSize> &v) const;
  // Bitwise AND
  inline LargeBinVect<CountWordSize> operator&(const LargeBinVect<CountWordSize> &v) const;
  // Right to Left shift by 's' bits
  inline LargeBinVect<CountWordSize> operator<<(const CountWordSize s) const;
  // Left to Right shift by 's' bits
  inline LargeBinVect<CountWordSize> operator>>(const CountWordSize s) const;

  // Bitwise XOR with self assignment
  inline LargeBinVect<CountWordSize> &operator^=(const LargeBinVect<CountWordSize> &v);
  // Bitwise OR with self assignment
  inline LargeBinVect<CountWordSize> &operator|=(const LargeBinVect<CountWordSize> &v);
  // Bitwise AND with self assignment
  inline LargeBinVect<CountWordSize> &operator&=(const LargeBinVect<CountWordSize> &v);
  // Right to Left shift by 's' bits with self assignment
  inline LargeBinVect<CountWordSize> &operator<<=(const CountWordSize s);
  // Left to Right shift by 's' bits with self assignment
  inline LargeBinVect<CountWordSize> &operator>>=(const CountWordSize s);

  // Sub-LargeBinVector starting at position b (>= 0) and
  // of length lg (or to the end)
  inline LargeBinVect<CountWordSize> operator()(const CountWordSize b = 0, const CountWordSize lg = ~((CountWordSize)0)) const throw (std::bad_alloc, std::out_of_range);

  // Optimized already casted version
  inline unsigned long int GetValue(const CountWordSize b = 0, const CountWordSize lg = ~((CountWordSize)0)) const throw (std::bad_alloc, std::out_of_range);

  // [ReadOnly] Access bit at position i (>= 0, < length)
  inline bool operator[](const CountWordSize i) const throw (std::out_of_range);

  // Bitwise equality. False if *this and v doesn't have the same length.
  inline bool operator==(const LargeBinVect<CountWordSize> &v) const;
  // Bitwise inequality. True if *this and v doesn't have the same length.
  inline bool operator!=(const LargeBinVect<CountWordSize> &v) const;
  // less than lexicographic operator
  inline bool operator<(const LargeBinVect<CountWordSize> &v) const;
  // greater than lexicographic operator
  inline bool operator>(const LargeBinVect<CountWordSize> &v) const;
  // less than or equal to lexicographic operator
  inline bool operator<=(const LargeBinVect<CountWordSize> &v) const;
  // greater than or equal to lexicographic operator
  inline bool operator>=(const LargeBinVect<CountWordSize> &v) const;

  // equivalent to operator[]
  inline bool GetBit(const CountWordSize i) const throw (std::out_of_range);
  // Set Bit at position i (>= 0, < length) to value
  inline LargeBinVect<CountWordSize> &SetBit(const CountWordSize i, const bool value) throw (std::out_of_range);

  // Set Bits between positions i and j (>= 0, < length) to value
  inline LargeBinVect<CountWordSize> &SetBitRange(const CountWordSize i, const CountWordSize j, const long unsigned int value) throw (std::out_of_range);

  // Clear the bit vector and set all bits to value
  inline LargeBinVect<CountWordSize> &Clear(const bool value = 0);

  // set each bit of the vector to a random value
  inline LargeBinVect<CountWordSize> &RandomInit();

  // Resize the bit vector according to length and set all extra bits to value (if any)
  inline LargeBinVect<CountWordSize> &Resize(const CountWordSize length, const bool value = 0) throw (std::bad_alloc);

  // Append a bit to the vector and set its value to value
  inline LargeBinVect<CountWordSize> &Append(const bool value) throw (BinVectOverflow, std::bad_alloc);

  // Prepend a bit to the vector and set its value to value
  inline LargeBinVect<CountWordSize> &Prepend(const bool value) throw (BinVectOverflow, std::bad_alloc);

  // Concatenate v to the current binary vector
  inline LargeBinVect<CountWordSize> &operator+=(const LargeBinVect<CountWordSize> &v) throw (BinVectOverflow, std::bad_alloc);

  // Pre-Increment the current binary vector by 1
  inline LargeBinVect<CountWordSize> &operator++() throw (BinVectOverflow, std::bad_alloc);
  // Post-Increment the current binary vector by 1
  inline LargeBinVect<CountWordSize> operator++(int) throw (BinVectOverflow, std::bad_alloc);

  // Pre-Decrement the current binary vector by 1
  inline LargeBinVect<CountWordSize> &operator--() throw (BinVectOverflow);
  // Post-Decrement the current binary vector by 1
  inline LargeBinVect<CountWordSize> operator--(int) throw (BinVectOverflow, std::bad_alloc);

  // Split the bit vector at position i.
  // Current object is set to the prefix of length i and the remaining suffix is returned.
  inline LargeBinVect<CountWordSize> Split(const CountWordSize i) throw (std::out_of_range, std::bad_alloc);

  // Get longuest common prefix length.
  inline CountWordSize GetLonguestCommonPrefix(const LargeBinVect<CountWordSize> &v) const;
  // Get longuest common suffix length.
  inline CountWordSize GetLonguestCommonSuffix(const LargeBinVect<CountWordSize> &v) const;

  // Get bit vector length.
  inline CountWordSize GetLength() const;

  // Convert to ascii
  inline char *ToString() const;

  // Cast operators
  inline operator char() const;
  inline operator unsigned char() const;
  inline operator short int() const;
  inline operator unsigned short int() const;
  inline operator int() const;
  inline operator unsigned int() const;
  inline operator long int() const;
  inline operator unsigned long int() const;
};

template <typename CountWordSize>
std::ostream &operator<<(std::ostream &os, const LargeBinVect<CountWordSize> &v);
  
#include "largebinvect.cxx"
#endif
