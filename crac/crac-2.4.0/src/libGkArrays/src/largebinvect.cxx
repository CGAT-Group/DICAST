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

#ifndef __LARGEBINVECT_CXX__
#define __LARGEBINVECT_CXX__

#include "largebinvect.hxx"

using std::ostream;

template <typename CountWordSize>
CountWordSize LargeBinVect<CountWordSize>::GetDataArrayLength() const {
  CountWordSize w = length/8 + ((length % 8) > 0);
  return w;  
}

template <typename CountWordSize>
unsigned char LargeBinVect<CountWordSize>::GetDataArrayEndingMask() const {
  unsigned char mask = (length % 8 ? 255 << (8 - (length % 8)) : 255);
  return mask;
}

template <typename CountWordSize>
unsigned char *MyCalloc(const CountWordSize w, const int value = -1) throw (std::bad_alloc) {
  unsigned char *array;
  if (!w) {
    return NULL;
  }
  if ((array = (unsigned char *) malloc(w*sizeof(unsigned char))) == NULL) {
    throw std::bad_alloc();
  }
  if (value >= 0) {
    memset(array, value, w*sizeof(unsigned char));
  }
  return array;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::LargeBinVect(const CountWordSize length, const bool value) throw (std::bad_alloc)
  : length(length) {
  CountWordSize w = GetDataArrayLength();
  unsigned char mask = GetDataArrayEndingMask();
  data = MyCalloc(w, (value ? 255 : 0));
  if (w) {
    data[w - 1] |= ~mask;
  }
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::LargeBinVect(const char *binary_string) throw (BinVectOverflow, std::bad_alloc) {
  length = strlen(binary_string);
  if (length != (long unsigned int) length) {
    throw BinVectOverflow();
  }
  CountWordSize w = GetDataArrayLength();
  unsigned char mask = GetDataArrayEndingMask();
  data = MyCalloc(w, 0);
  for (CountWordSize i = 0; i < length; i++) {
    if (binary_string[i] == '1') {
      unsigned char m = 1 << (7 - (i%8));
      data[i/8] |= m;
    }
  }
  if (w) {
    data[w - 1] |= ~mask;
  }
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::LargeBinVect(const LargeBinVect<CountWordSize> &v) throw (std::bad_alloc)
  : length(v.length)  {
  CountWordSize w = v.GetDataArrayLength();
  data = MyCalloc(w);
  for (CountWordSize i = 0; i < w; i++) {
    data[i] = v.data[i];
  }
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::~LargeBinVect() {
  if (data) {
    free(data);
  }
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator=(const LargeBinVect<CountWordSize> &v) throw (std::bad_alloc) {
  if (this != &v) {
    if (data) {
      free(data);
    }
    length = v.length;
    CountWordSize w = v.GetDataArrayLength();
    data = MyCalloc(w);
    for (CountWordSize i = 0; i < w; i++) {
      data[i] = v.data[i];
    }
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator~() const {
  LargeBinVect<CountWordSize> bv(*this);
  CountWordSize w = bv.GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    bv.data[i] ^= 255;
  }
  unsigned char mask = GetDataArrayEndingMask();
  if (w) {
    bv.data[w - 1] |= ~mask;
  }
  return bv;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator^(const LargeBinVect<CountWordSize> &v) const {
  assert(v.length == length);
  LargeBinVect<CountWordSize> bv(*this);
  CountWordSize w = bv.GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    bv.data[i] ^= v.data[i];
  }
  unsigned char mask = GetDataArrayEndingMask();
  if (w) {
    bv.data[w - 1] |= ~mask;
  }
  return bv;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator|(const LargeBinVect<CountWordSize> &v) const {
  assert(v.length == length);
  LargeBinVect<CountWordSize> bv(*this);
  CountWordSize w = bv.GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    bv.data[i] |= v.data[i];
  }
  return bv;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator&(const LargeBinVect<CountWordSize> &v) const {
  assert(v.length == length);
  LargeBinVect<CountWordSize> bv(*this);
  CountWordSize w = bv.GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    bv.data[i] &= v.data[i];
  }
  return bv;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator<<(const CountWordSize s) const {
  LargeBinVect<CountWordSize> v(*this);
  v<<=s;
  return v;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator>>(const CountWordSize s) const {
  LargeBinVect<CountWordSize> v(*this);
  v>>=s;
  return v;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator^=(const LargeBinVect<CountWordSize> &v) {
  assert(v.length == length);
  CountWordSize w = GetDataArrayLength();
  unsigned char mask = GetDataArrayEndingMask();
  for (CountWordSize i = 0; i < w; i++) {
    data[i] ^= v.data[i];
  }
  if (w) {
    data[w - 1] |= ~mask;
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator|=(const LargeBinVect<CountWordSize> &v) {
  assert(v.length == length);
  CountWordSize w = GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    data[i] |= v.data[i];
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator&=(const LargeBinVect<CountWordSize> &v) {
  assert(v.length == length);
  CountWordSize w = GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    data[i] &= v.data[i];
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator<<=(const CountWordSize s) {
  CountWordSize w = GetDataArrayLength();
  unsigned char m = ~(255>>s);
  unsigned char mask = GetDataArrayEndingMask();
  for (CountWordSize i = 0; i < w - 1; i++) {
    if (i) {
      data[i - 1] |= (data[i] & m) >> (8-s);
    }
    data[i] <<= s;
  }
  if (w) {
    if (w - 1) {
      data[w - 2] |= (data[w - 1] & m & mask) >> (8-s);
    }
    data[w - 1] &= mask;
    data[w - 1] <<= s; 
    data[w - 1] |= ~mask;
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator>>=(const CountWordSize s) {
  CountWordSize w = GetDataArrayLength();
  unsigned char m = ~(255<<s);
  unsigned char mask = GetDataArrayEndingMask();
  if (!w) {
    return *this;
  }
  data[w - 1] &= mask;
  data[w - 1] >>= s;
  for (CountWordSize i = w - 2; (i >= 0) && (i < w - 1); i--) {
    data[i + 1] |= (data[i] & m) << (8-s);
    data[i] >>= s;
  }
  data[w - 1] |= ~mask;
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator()(const CountWordSize b, const CountWordSize lg) const  throw (std::bad_alloc, std::out_of_range) {
  if ((b < 0) || (b > length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "b=" << (unsigned long int) b
	<< " should be >= 0 and <= " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  CountWordSize l = (lg > length ? length : lg);
  if (b+l > length) {
    l = length - b;
  }
  assert(l <= length || l == 0);
  LargeBinVect<CountWordSize> bv(l);
  for (CountWordSize i = 0; i < l; i++) {
    bv.SetBit(i, GetBit(i+b));
  }
  return bv;
}


template <typename CountWordSize>
unsigned long int LargeBinVect<CountWordSize>::GetValue(const CountWordSize b, const CountWordSize lg) const throw (std::bad_alloc, std::out_of_range) {
  if ((b < 0) || (b > length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "b=" << (unsigned long int) b
	<< " should be >= 0 and <= " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  unsigned long int r = 0;
  CountWordSize l = b + (lg > sizeof(unsigned long int)*8 ? sizeof(unsigned long int)*8 : lg);
//   std::cerr << "1: b = " << (unsigned int) b
// 	    << ", l = " << (unsigned int) l
// 	    << std::endl;
  if ((l > length) || (l < b)) {
    l = length;
  }
//   std::cerr << "2: b = " << (unsigned int) b
// 	    << ", l = " << (unsigned int) l
// 	    << std::endl;
  assert(l <= length || l == 0);
  for (CountWordSize i = b; i < l; i++) {
    r <<= 1;
    r |= GetBit(i);
//     std::cerr << "i = " << (unsigned int) i
// 	      << ", GetBit(i) = " << GetBit(i)
// 	      << ", r = " << r
// 	      << std::endl;
  }
  return r;
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator[](const CountWordSize i) const throw (std::out_of_range) {
  if ((i < 0) || (i >= length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "i=" << (unsigned long int) i
	<< " should be >= 0 and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  return (128&(data[i/8]<<(i%8)));
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator==(const LargeBinVect<CountWordSize> &v) const {
  if (length != v.length) {
    return false;
  }
  const CountWordSize w = GetDataArrayLength();
  for (CountWordSize i = 0; i < w; i++) {
    if (data[i] ^ v.data[i]) { // True if they differ by at least one bit
      return false;
    }
  }
  return true;
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator!=(const LargeBinVect<CountWordSize> &v) const {
  return !(*this == v);
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator<(const LargeBinVect<CountWordSize> &v) const {
  const CountWordSize &lg = (length > v.length ? v.length : length);
  const CountWordSize w = lg/8;
  CountWordSize i;
  for (i = 0; i < w; i++) {
    if (data[i] ^ v.data[i]) { // True if they differ by at least one bit
      return (data[i] < v.data[i]);
    }
  }
  if (lg%8) {
    unsigned char mask = (length > v.length ? v.GetDataArrayEndingMask() : GetDataArrayEndingMask());
    if ((data[i]&mask) ^ (v.data[i]&mask)) { // True if they differ by at least one bit
     return ((data[i]&mask) < (v.data[i]&mask));
    }
  }
  return (length < v.length);
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator>(const LargeBinVect<CountWordSize> &v) const {
  return v < *this;
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator<=(const LargeBinVect<CountWordSize> &v) const {
  const CountWordSize &lg = (length > v.length ? v.length : length);
  const CountWordSize w = lg/8;
  CountWordSize i;
  for (i = 0; i < w; i++) {
    if (data[i] ^ v.data[i]) { // True if they differ by at least one bit
      return (data[i] < v.data[i]);
    }
  }
  if (lg%8) {
    unsigned char mask = (length > v.length ? v.GetDataArrayEndingMask() : GetDataArrayEndingMask());
    if ((data[i]&mask) ^ (v.data[i]&mask)) { // True if they differ by at least one bit
     return ((data[i]&mask) < (v.data[i]&mask));
    }
  }
  return (length <= v.length);
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::operator>=(const LargeBinVect<CountWordSize> &v) const {
  return v <= *this;
}

template <typename CountWordSize>
bool LargeBinVect<CountWordSize>::GetBit(const CountWordSize i) const throw (std::out_of_range) {
  if ((i < 0) || (i >= length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "i=" << (unsigned long int) i
	<< " should be >= 0 and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  return (data[i/8] & (1 << (7-(i%8))));
}


template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::SetBit(const CountWordSize i, const bool value) throw (std::out_of_range) {
  if ((i < 0) || (i >= length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "i=" << (unsigned long int) i
	<< " should be >= 0 and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  if (value) {
    data[i/8] |= (1 << (7-(i%8)));
  } else {
    data[i/8] &= ~(1 << (7-(i%8)));
  }
  return *this;
}

// Set Bits between positions i and j (>= 0, < length) to value
template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::SetBitRange(const CountWordSize i, const CountWordSize j, const long unsigned int value) throw (std::out_of_range) {
  if ((i < 0) || (i >= length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "i=" << (unsigned long int) i
	<< " should be >= 0 and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  if ((j < i) || (j >= length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "j=" << (unsigned long int) j
	<< " should be >= " << i << " and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  long unsigned int v = value;
//   cerr << "i / 8 = " << (unsigned int) i <<" / 8 = " << i / 8 << " et i % 8 = " << i % 8 << endl;
//   cerr << "j / 8 = " << (unsigned int) j <<" / 8 = " << j / 8 << " et j % 8 = " << j % 8 << endl;
  if ((j / 8) > (i / 8)) {
//     cerr << "cas 1" << endl;
    if ((j + 1) % 8) {
//       cerr << "cas 1A" << endl;
      unsigned char mask = ~(255 << (j % 8 + 1));
      unsigned char clear_mask = ~(255 << (7 - (j % 8)));
//       cerr << "clearing mask = " << bitset<8>(clear_mask) << ", setting mask = " << bitset<8>(mask) << endl;
//       cerr << "data["<<(unsigned int) j<<" / 8] = " << bitset<8>(data[j / 8]) << endl;
      data[j / 8] &= clear_mask;
//       cerr << "clearing => " << bitset<8>(data[j / 8]) << endl;
      data[j / 8] |= ((char) v & mask) << (7 - (j % 8));
//       cerr << "setting with " << bitset<8>(((char) v & mask)) << "<<" << 7 - (j % 8) << " => " << bitset<8>(data[j / 8]) << endl;
      v >>= (j % 8 + 1);
//       cerr << "shifting value by " << j % 8 + 1 << " bits => " << bitset<8>(v) << endl;
    }
    for (CountWordSize x = j / 8 - ((j + 1) % 8 != 0); x > i / 8; x--) {
//       cerr << "cas 1B for data[x = " << (unsigned int) x << "]" << endl;      
      data[x] = (char) v & 255;
      v >>= 8;
//       cerr << "shifting value by 8 bits => " << bitset<8>(v) << endl;
    }
//     cerr << "cas 1C" << endl;      
    unsigned char mask = 255 << (8 - (i % 8));
//     cerr << "clearing mask = " << bitset<8>(mask) << ", setting mask = " << bitset<8>(~mask) << endl;
//     cerr << "data[" << (unsigned int) i <<" / 8] = " << bitset<8>(data[i / 8]) << endl;
    data[i / 8] &= mask;
//     cerr << "clearing => " << bitset<8>(data[i / 8]) << endl;
    data[i / 8] |= (char) v & ~mask;
//     cerr << "setting => " << bitset<8>(data[i / 8]) << endl;
  } else {
//     cerr << "cas 2" << endl;
    if (i / 8 == length / 8) {
//       cerr << "cas 2A" << endl;
      unsigned char mask = ~(255 << (j - i + 1));
      unsigned char clear_mask = ~(mask << (8 - (length % 8)));
//       cerr << "clearing mask is " << bitset<8>(clear_mask) << endl;
      data[i / 8] &= clear_mask;
//       cerr << "clearing => " << bitset<8>(data[i / 8]) << endl;
//       cerr << "value setting mask = " << bitset<8>(mask) << endl;
//       cerr << "result of masking will be shifted by " << (8 - (length % 8)) << " positions" << endl;
      data[i / 8] |= (v & mask) << (8 - (length % 8));
//       cerr << "setting => " << bitset<8>(data[i / 8]) << endl;
    } else {
//       cerr << "cas 2B" << endl;
      unsigned char clearing_mask = ~((255 >> (i % 8)) & (255 << (7 - (j % 8))));
      unsigned char setting_mask = ~(255 << (j - i + 1));
//       cerr << "clearing mask = " << bitset<8>(clearing_mask) << ", setting mask = " << bitset<8>(setting_mask) << endl;
//       cerr << "data[" << (unsigned int) i <<" / 8] = " << bitset<8>(data[i / 8]) << endl;
      data[i / 8] &= clearing_mask;
//       cerr << "clearing => " << bitset<8>(data[i / 8]) << endl;
      data[i / 8] |= ((char) v & setting_mask) << (7 - (j % 8));
//       cerr << "setting => " << bitset<8>(data[i / 8]) << endl;
    }
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::Clear(const bool value) {
  unsigned char mask = GetDataArrayEndingMask();
  CountWordSize w = GetDataArrayLength();
  if (w) {
    memset(data, (value ? 255 : 0), w*sizeof(unsigned char));
    data[w - 1] |= ~mask;
  }
  return *this;
}


template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::RandomInit() {
  unsigned char mask = GetDataArrayEndingMask();
  CountWordSize w = GetDataArrayLength();
  if (w) {
    for (CountWordSize i = 0; i < w; i++) {
      data[i] = (unsigned char) random();
    }
    data[w - 1] |= ~mask;
  }
  return *this;
}


template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::Resize(const CountWordSize length, const bool value) throw (std::bad_alloc) {
  CountWordSize old_w = GetDataArrayLength();
  unsigned char mask = GetDataArrayEndingMask();
  this->length = length;
  CountWordSize w = GetDataArrayLength();

  if (!w) { // If length is 0 ;-)
    free(data);
    data = NULL;
    return *this;
  }

  if (old_w) {
    if (value) {
      data[old_w - 1] |= ~mask;
    } else {
      data[old_w - 1] &= mask;
    }
  }

  if ((data = (unsigned char *) realloc(data, w*sizeof(unsigned char))) == NULL) {
    throw std::bad_alloc();
  }
  if (w > old_w) {
    memset(data+old_w, (value ? 255 : 0), (w-old_w)*sizeof(unsigned char));
  }
  mask = GetDataArrayEndingMask();
  data[w - 1] |= ~mask;
  
  return *this;
  
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::Append(const bool value) throw (BinVectOverflow, std::bad_alloc) {
  if (length == CountWordSize(~0)) {
    throw BinVectOverflow();
  }
  return Resize(length+1, value);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::Prepend(const bool value) throw (BinVectOverflow, std::bad_alloc) {
  if (length == CountWordSize(~0)) {
    throw BinVectOverflow();
  }
  Resize(length+1) >>= 1;
  if (value) {
    data[0] |= 128;
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator+=(const LargeBinVect<CountWordSize> &v) throw (BinVectOverflow, std::bad_alloc) {
  if ((length + v.length < length) || (length + v.length < v.length)) {
    throw BinVectOverflow();
  }
  if (!v.length) {
    return *this;
  }
  CountWordSize old_w = GetDataArrayLength();
  unsigned char mask = GetDataArrayEndingMask();
  CountWordSize app_w = v.GetDataArrayLength();
  unsigned int shift = (length%8);
  Resize(length + v.length);
  if (shift) {
    assert(old_w);
    assert(app_w);
    old_w--;
    data[old_w] &= mask;
    mask = v.GetDataArrayEndingMask();
    v.data[app_w - 1] &= mask;
    for (CountWordSize i = 0; i < app_w - 1; i++) {
      data[old_w + i] |= v.data[i] >> shift;
      data[old_w + i + 1] |= v.data[i] << (8 - shift);
    }
    data[old_w + app_w - 1] |= v.data[app_w - 1] >> shift;
    if (old_w + app_w == GetDataArrayLength() - 1) {
      data[old_w + app_w] |= v.data[app_w - 1] << (8 - shift);
    }
  } else {
    for (CountWordSize i = 0; i < app_w; i++) {
      data[old_w + i] = v.data[i];
    }
    mask = GetDataArrayEndingMask();
    data[old_w + app_w - 1] |= ~mask;  
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator++() throw (BinVectOverflow, std::bad_alloc) {
  CountWordSize w = GetDataArrayLength();
  unsigned char mask = ~GetDataArrayEndingMask();
  ++mask;
  while (mask && w--) {
    unsigned char old = data[w];
    data[w] += mask;
    mask = (data[w] < old);
  }
  if (mask) {
    Prepend(true);
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator++(int) throw (BinVectOverflow, std::bad_alloc) {
  LargeBinVect<CountWordSize> ret(*this);
  ++(*this);
  return ret;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> &LargeBinVect<CountWordSize>::operator--() throw (BinVectOverflow) {
  CountWordSize w = GetDataArrayLength();
  unsigned char mask = ~GetDataArrayEndingMask();
  ++mask;
  while (mask && w--) {
    unsigned char old = data[w];
    data[w] -= mask;
    mask = (data[w] > old);
  }
  if (mask) {
    throw BinVectOverflow();
  }
  return *this;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::operator--(int) throw (BinVectOverflow, std::bad_alloc) {
  LargeBinVect<CountWordSize> ret(*this);
  --(*this);
  return ret;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize> LargeBinVect<CountWordSize>::Split(const CountWordSize i) throw (std::out_of_range, std::bad_alloc) {
  if ((i < 0) || (i > length)) {
    std::stringstream msg;
    msg << __FILE__ << ":"
	<< __PRETTY_FUNCTION__ << ":"
	<< __LINE__ << ":"
	<< "i=" << (unsigned long int) i
	<< " should be >= 0 and < " << (unsigned long int) length;
    throw std::out_of_range(msg.str());
  }
  if (i == length) {
    return LargeBinVect<CountWordSize>();
  } else {
    LargeBinVect<CountWordSize> bv = (*this)(i, length - i);
    Resize(i);
    return bv;
  }  
}

template <typename CountWordSize>
CountWordSize LargeBinVect<CountWordSize>::GetLonguestCommonPrefix(const LargeBinVect<CountWordSize> &v) const {
  const CountWordSize &lg = (length > v.length ? v.length : length);
  const CountWordSize w = lg/8;
  CountWordSize i;
  for (i = 0; i < w; i++) {
    unsigned char m = data[i] ^ v.data[i];
    CountWordSize p;
    if (m) { // True if they differ by at least one bit
      for (p = i*8; !(128 & m); p++) {
	m <<= 1;
      }
      return p;
    }
  }
  if (lg%8) {
    unsigned char mask = (length > v.length ? v.GetDataArrayEndingMask() : GetDataArrayEndingMask());
    unsigned char m = (data[i] & mask) ^ (v.data[i] & mask);
    CountWordSize p;
    if (m) { // True if they differ by at least one bit
      for (p = i*8; !(128 & m); p++) {
	m <<= 1;
      }
      return p;
    }
  }
  return lg;
}

template <typename CountWordSize>
CountWordSize LargeBinVect<CountWordSize>::GetLonguestCommonSuffix(const LargeBinVect<CountWordSize> &v) const {
  if (length > v.length) {
    return v.GetLonguestCommonSuffix(*this);
  } else { // length <= v.length
    const CountWordSize delta = v.length - length;
    for (CountWordSize i = length - 1; (i >= 0) && (i < length); i--) {
      if (GetBit(i) != v.GetBit(i+delta)) {
	return length - i - 1;
      }
    }
    return length;
  }
}

template <typename CountWordSize>
char * LargeBinVect<CountWordSize>::ToString() const {
  char* tmp;
  tmp = new char [length + 1];
  const CountWordSize w = length / 8;
  unsigned char m;
  for (CountWordSize i = 0; i < w; i++) {
    m = data[i];
    for (unsigned char b = 0; b < 8; b++) {
      tmp[i*8+b] = ((128 & m) ? '1' : '0');
      m<<=1;
    }
  }
  if (length % 8) {
    m = data[w];
    for (unsigned char b = 0; b < length % 8; b++) {
      tmp[w*8+b] = ((128 & m) ? '1' : '0');
      m<<=1;
    }
  }
  tmp[length]=0;
  return tmp;
}

template <typename CountWordSize>
CountWordSize LargeBinVect<CountWordSize>::GetLength() const {
  return length;
}

// Cast operators
template <typename CountWordSize, typename target_type>
target_type cast(const LargeBinVect<CountWordSize> &v) {
  target_type r = 0;
  CountWordSize i = sizeof(target_type)*8;
  if (i > v.GetLength()) {
    i = v.GetLength();
  }
  for (; i > 0; i--) {
    r <<= 1;
    r |= v.GetBit(v.GetLength() - i);
  }
  return r;
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator char() const {
  return cast<CountWordSize, char>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator unsigned char() const {
  return cast<CountWordSize, unsigned char>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator short int() const {
  return cast<CountWordSize, short int>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator unsigned short int() const {
  return cast<CountWordSize, unsigned short int>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator int() const {
  return cast<CountWordSize, int>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator unsigned int() const {
  return cast<CountWordSize, unsigned int>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator long int() const {
  return cast<CountWordSize, long int>(*this);
}

template <typename CountWordSize>
LargeBinVect<CountWordSize>::operator unsigned long int() const {
  return cast<CountWordSize, unsigned long int>(*this);
}

template <typename CountWordSize>
std::ostream &operator<<(std::ostream &os, const LargeBinVect<CountWordSize> &v) {
  char *tmp =  v.ToString();
  os << tmp;
  delete [] tmp;
  return os;
}

#ifdef TEST_LARGEBINVECT
// Compile this test program with the following command line:
// g++ largebinvect.cxx -DTEST_LARGEBINVECT -o test_largebinvect -Wall -ansi -pedantic -O4

using std::cout;
using std::endl;
#include <bitset>

template <typename CountWordSize>
void TestLargeBinVectors(char *title, CountWordSize max_test_length) {
  cout << "Testing class LargeBinVect<" << title << ">" << endl;
  cout << "SizeOf(" << title << ") = "
       << sizeof(CountWordSize) << "Bytes (= "
       << 8*sizeof(CountWordSize) << "bits)" << endl << endl;
  const char *binary_string1 = "1001010111010001011110100010101";
  const char *binary_string2 = "10010101110100010101";
  cout << "Building a short binvect from string '" << binary_string1 << "'" << endl;
  LargeBinVect<CountWordSize> * v1 = new LargeBinVect<CountWordSize>(binary_string1);
  cout << "=> '" << *v1 << "' of length " << (unsigned long int) v1->GetLength() << endl;
  LargeBinVect<CountWordSize> * v2 = new LargeBinVect<CountWordSize>(binary_string2);
  cout << "=> '" << *v2 << "' of length " << (unsigned long int) v2->GetLength() << endl;
  cout << "Compute the longuest common prefix (size 18) between" << endl;
  cout << "bv1='" << *v1 << "'" << endl;
  cout << "bv2='" << *v2 << "'" << endl;      
  cout << "It returns " << (unsigned long int)v1->GetLonguestCommonPrefix(*v2);
  cout << " which should also be " << (unsigned long int)v2->GetLonguestCommonPrefix(*v1) << endl;

  cout << "Compute the longuest common suffix (size 13) between" << endl;
  cout << "bv1='" << *v1 << "'" << endl;
  cout << "bv2='" << *v2 << "'" << endl;      
  cout << "It returns " << (unsigned long int)v1->GetLonguestCommonSuffix(*v2);
  cout << " which should also be " << (unsigned long int)v2->GetLonguestCommonSuffix(*v1) << endl;
  cout << endl;

  cout << "Resizing bv1 to 78 bits (filling with 0)" << endl;
  cout << "bv1 was ='" << *v1 << "'" << endl;
  v1->Resize(78);
  cout << "bv1 is  ='" << *v1 << "'" << endl;
  cout << "Resizing bv2 to 35 bits (filling with 1)" << endl;
  cout << "bv2 was ='" << *v2 << "'" << endl;
  v2->Resize(35, 1);
  cout << "bv2 is  ='" << *v2 << "'" << endl;
  cout << "Resizing bv1 to 35 bits ([not] filling with 1)" << endl;
  cout << "bv1 was ='" << *v1 << "'" << endl;
  v1->Resize(35, 1);
  cout << "bv1 is  ='" << *v1 << "'" << endl;
  cout << "Resizing bv2 to 13 bits ([not] filling with 0)" << endl;
  cout << "bv2 was ='" << *v2 << "'" << endl;
  v2->Resize(13);
  cout << "bv2 is  ='" << *v2 << "'" << endl;
  cout << "Resizing bv2 to 20 bits (filling with 0)" << endl;
  cout << "bv2 was ='" << *v2 << "'" << endl;
  v2->Resize(20);
  cout << "bv2 is  ='" << *v2 << "'" << endl;
  cout << endl;

  *v2 = LargeBinVect<CountWordSize>("1000000010000000100000001000000010000000100000001000000010000000");

  cout << "Testing cast operators" << endl;
  cout << "casting " << *v2 << " to char:" << (int) ((char) *v2)
       << ", which should be -128" << endl;
  cout << "casting " << *v2 << " to unsigned char:" << (unsigned int) (unsigned char) *v2
       << ", which should be 128" << endl;
  cout << "casting " << *v2 << " to short:" << (short int) *v2
       << ", which should be -32640" << endl;
  cout << "casting " << *v2 << " to unsigned short:" << (unsigned short int) *v2
       << ", which should be 32896" << endl;
  cout << "casting " << *v2 << " to int:" << (int) *v2
       << ", which should be -2139062144" << endl;
  cout << "casting " << *v2 << " to unsigned int:" << (unsigned int) *v2
       << ", which should be 2155905152" << endl;
  cout << "casting " << *v2 << " to long:" << (long int) *v2
       << ", which should be -9187201950435737472" << endl;
  cout << "casting " << *v2 << " to unsigned long:" << (unsigned long int) *v2
       << ", which should be 9259542123273814144" << endl;
  cout << endl;

  cout << "Deleting the binvectors" << endl;
  delete v1;
  delete v2;

  cout << endl;
  cout << "Creating empty binvector" << endl;
  v1 = new LargeBinVect<CountWordSize>();
  cout << "LargeBinVect is '" << *v1 << "'" << endl;
  cout << "Preincrementing '" << *v1 << "'" << endl;
  cout << "=> '" << ++(*v1) << "'" << endl;
  cout << "Postincrementing '" << (*v1)++ << "'" << endl;
  cout << "=> '" << *v1 << "'" << "'" << endl;
  cout << "Incrementing '" << *v1 << "'" << endl;
  cout << "=> '" << ++(*v1) << "'" << endl;
  cout << "Incrementing '" << *v1 << "'" << endl;
  cout << "=> '" << ++(*v1) << "'" << endl;
  cout << "Appending '111111111111111' to '" << *v1 << "'" << endl;
  *v1 += LargeBinVect<CountWordSize>("111111111111111");
  cout << "=> '" << *v1 << "'" << endl;
  cout << "Incrementing '" << *v1 << "'" << endl;
  cout << "=> '" << ++(*v1) << "'" << endl;

  cout << "Predecrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Postdecrementing '" << (*v1)-- << "'" << endl;
  cout << "=> '" << *v1 << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Resizing '" << *v1 << "' to 3 bits" << endl;
  v1->Resize(3);
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << "Decrementing '" << *v1 << "'" << endl;
  cout << "=> '" << --(*v1) << "'" << endl;
  cout << endl;
  cout << "Deleting the binvector" << endl;
  delete v1;

  cout << "===" << endl;
  for (CountWordSize i = 1; i < max_test_length; i*=2) {
    cout << endl;
    cout << "Make 0-filled binvect of length " << (unsigned long int) i << endl;
    LargeBinVect<CountWordSize> bv(i);
    cout << "=> '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Overwrite (operator=) the binvect with a 1-filled one of same length" << endl;
    bv = LargeBinVect<CountWordSize>(i, 1);
    cout << "=> '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Clearing the binvect with 0" << endl;
    bv.Clear();
    cout << "=> '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Clearing the binvect with 1" << endl;
    bv.Clear(1);
    cout << "=> '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    long unsigned int r = rand();
    cout << "Setting the bit range from " << i / 3
	 << " to " << 3 * i / 4
	 << " with the value " << std::bitset<20>(r) << endl;
    bv.SetBitRange(i / 3, 3 * i / 4, r);
    cout << "=> '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Setting even bits to 0" << endl;
    cout << "was '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;
    for (CountWordSize b = 0; b < i; b+=2) {
      bv.SetBit(b, 0);
    }
    cout << "=>  '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Left Shifting (Self assignment) binvect by 3" << endl;
    cout << "was '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;
    bv <<= 3;
    cout << "=>  '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Right Shifting (Self assignment) binvect by 4" << endl;
    cout << "was '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;
    bv >>= 4;
    cout << "=>  '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << "Assigning its inversed bit values" << endl;
    cout << "was '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;
    bv = ~bv;
    cout << "=>  '" << bv << "' of length " << (unsigned long int) bv.GetLength() << endl;

    cout << endl;
    cout << "Generating a random bitvector having same length" << endl;
    LargeBinVect<CountWordSize> rbv(i);
    rbv.RandomInit();
    cout << "=> '" << rbv << "' of length " << (unsigned long int) rbv.GetLength() << endl;

    cout << "Performing XOR between bv and rbv" << endl;
    LargeBinVect<CountWordSize> res = (bv ^ rbv);
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Performing OR between bv and rbv" << endl;
    res = (bv | rbv);
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Performing AND between bv and rbv" << endl;
    res = (bv & rbv);
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;
   
    cout << "Using the copy constructor" << endl;
    v1 = new LargeBinVect<CountWordSize>(res);
    cout << "res='" << res << "'" << endl;      
    cout << "=>  '" << *v1 << "' of length " << (unsigned long int) v1->GetLength() << endl;
    delete v1;

    cout << "Performing Self Assignment XOR between copy of bv and rbv" << endl;
    res = bv;
    res ^= rbv;
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Performing Self Assignment OR between copy of bv and rbv" << endl;
    res = bv;
    res |= rbv;
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Performing Self Assignment AND between copy of bv and rbv" << endl;
    res = bv;
    res &= rbv;
    cout << "bv ='" << bv << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Build prefix of rbv of length " << (unsigned long int) rbv.GetLength()/3 << endl;
    res = rbv(0, rbv.GetLength()/3);
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "Testing between rbv < res" << endl;
    cout << "=> " << (rbv < res ? "true" : "false") << endl;
    cout << "Testing between rbv <= res" << endl;
    cout << "=> " << (rbv <= res ? "true" : "false") << endl;
    cout << "Testing between rbv > res" << endl;
    cout << "=> " << (rbv > res ? "true" : "false") << endl;
    cout << "Testing between rbv >= res" << endl;
    cout << "=> " << (rbv >= res ? "true" : "false") << endl;
    cout << "Testing between rbv == res" << endl;
    cout << "=> " << (rbv == res ? "true" : "false") << endl;
    cout << "Testing between rbv != res" << endl;
    cout << "=> " << (rbv != res ? "true" : "false") << endl;

    CountWordSize p = 0;
    if (res.GetLength()) {
      p = rand() % res.GetLength();
      cout << "Spin bit at position " << (unsigned long int) p << " in res (using GetBit, SetBit and operator[])" << endl;
      cout << "res["<<(unsigned long int) p<<"] was " << res[p];
      res.SetBit(p, !res.GetBit(p));
      cout << " and now  is " << res[p]<<endl;
      cout << "=> '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;
    }

    cout << endl;
    cout << "Computing the longuest common prefix between res and rbv should give " << (res.GetLength() ? (unsigned long int) p : 0) << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "It returns " << (unsigned long int)res.GetLonguestCommonPrefix(rbv);
    cout << " which should also be " << (unsigned long int)rbv.GetLonguestCommonPrefix(res) << endl;

    cout << endl;
    cout << "Build suffix of rbv of length " << (unsigned long int) rbv.GetLength()/3 << endl;
    res = rbv(rbv.GetLength() - rbv.GetLength()/3 - 1, rbv.GetLength()/3);
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "=>  '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    if (res.GetLength()) {
      p = rand() % res.GetLength();    
      cout << "Spin bit at position " << (unsigned long int) p << " in res (using GetBit, SetBit and operator[])" << endl;
      cout << "res["<<(unsigned long int) p<<"] was " << res[p];
      res.SetBit(p, !res.GetBit(p));
      cout << " and now  is " << res[p]<<endl;
      cout << "=> '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;
    }

    cout << endl;
    cout << "Computing the longuest common suffix between res and rbv" << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "It returns " << (unsigned long int)res.GetLonguestCommonSuffix(rbv);
    cout << " which should also be " << (unsigned long int)rbv.GetLonguestCommonSuffix(res) << endl;

    cout << endl;
    cout << "Build centered factor of rbv of length " << (unsigned long int) rbv.GetLength()/3 << endl;
    res = rbv(rbv.GetLength()/3, rbv.GetLength()/3);
    cout << "=> '" << res << "' of length " << (unsigned long int) res.GetLength() << endl;

    cout << "Compute the longuest common prefix between res and rbv." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "It returns " << (unsigned long int)res.GetLonguestCommonPrefix(rbv);
    cout << " which should also be " << (unsigned long int)rbv.GetLonguestCommonPrefix(res) << endl;

    cout << "Compute the longuest common suffix between res and rbv." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;      
    cout << "It returns " << (unsigned long int)res.GetLonguestCommonSuffix(rbv);
    cout << " which should also be " << (unsigned long int)rbv.GetLonguestCommonSuffix(res) << endl;
    cout << "===" << endl;

    cout << "Append longuest common prefix between res and rbv to the longuest common suffix between res and rbv." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "rbv='" << rbv << "'" << endl;
    cout << "* res.GetLength() = " << (unsigned long int)res.GetLength() << endl;
    cout << "* res.GetLonguestCommonSuffix(rbv) = " << (unsigned long int) res.GetLonguestCommonSuffix(rbv) << endl;
    cout << "* res.GetLength() - res.GetLonguestCommonSuffix(rbv) = " << (unsigned long int)(res.GetLength() - res.GetLonguestCommonSuffix(rbv)) << endl;
    bv = res(res.GetLength() - res.GetLonguestCommonSuffix(rbv));
    res = res(0, res.GetLonguestCommonPrefix(rbv));
    cout << "longuest common prefix is '" << res << "'" << endl;
    cout << "longuest common suffix is '" << bv << "'" << endl;
    res += bv;
    cout << "res is ='" << res << "'" << endl;    
    cout << "===" << endl;

    res = rbv;
    cout << "Append a 0 to res." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "now='" << res.Append(0) << "'" << endl;

    cout << "Append a 1 to res." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "now='" << res.Append(1) << "'" << endl;
    
    cout << "Prepend a 0 to res." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "now='" << res.Prepend(0) << "'" << endl;

    cout << "Prepend a 1 to res." << endl;
    cout << "res='" << res << "'" << endl;
    cout << "now='" << res.Prepend(1) << "'" << endl;

    cout << "Split res at position " << (unsigned long int) i/3 << "." << endl;
    cout << "res='" << res << "'" << endl;
    rbv = res.Split(i/3);
    cout << "now='" << res << "'" << endl;
    cout << "suf='" << rbv << "'" << endl;
  }
}
int main(int argc, char** argv) {

  TestLargeBinVectors<unsigned char>((char *)"unsigned char", 128);
  TestLargeBinVectors<unsigned int>((char *)"unsigned int", 1000000);
  TestLargeBinVectors<long unsigned int>((char *)"long unsigned int", 1000000);
  
  return 0;
}
#endif
#endif
