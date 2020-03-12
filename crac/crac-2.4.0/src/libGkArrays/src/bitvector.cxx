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

#ifndef __BITVECTOR_CXX__
#define __BITVECTOR_CXX__

#include "bitvector.hxx"

#define BV BitVector<CountWordSize, InnerWordSize>
#define NO_TYPE
#define tpl(ret_type, fct) \
template <typename CountWordSize, typename InnerWordSize> \
ret_type BV::fct

tpl(NO_TYPE, BitVector(const CountWordSize length, const InnerWordSize nb_bits)):nb_bits(nb_bits), data(length*nb_bits, 0) {
  assert((((unsigned long int)length)*nb_bits) <= (unsigned long int)(InnerWordSize(~0)));
}

tpl(NO_TYPE, BitVector(const BitVector &v)): nb_bits(v.nb_bits), data(v.data) {
}

tpl(NO_TYPE, ~BitVector()) {
}

tpl(BV &,operator=(const BitVector &v)) {
  if (this != &v) {
    nb_bits = v.nb_bits;
    data = v.data;
  }
  return *this;
}

tpl(inline LargeBinVect<InnerWordSize>,  operator[](const CountWordSize i) const) {
  return data(i*nb_bits, nb_bits);
}

tpl(unsigned long int, GetValue(const CountWordSize i) const) {
  return data.GetValue(i*nb_bits, nb_bits);
}


// Set value at position i (>= 0, < length) to value
tpl(BV &, OldSetValue(const CountWordSize i, const InnerWordSize value)) {
  InnerWordSize v = value;
  for (InnerWordSize j = nb_bits; j > 0; j--) {
    data.SetBit(i * nb_bits + j - 1, 1&v);
    v >>= 1;
  }
  return *this;
}

// Set value at position i (>= 0, < length) to value
tpl(BV &, SetValue(const CountWordSize i, const InnerWordSize value)) {
  data.SetBitRange(i * nb_bits, (i + 1) * nb_bits - 1, value);
  return *this;
}

tpl(BV &, Clear()) {
  data.Clear(0);
  return *this;
}

tpl(BV &, Resize(const CountWordSize length)) {
  data.Resize(length*nb_bits);
  return *this;
}

tpl(BV &, CompleteResize(const CountWordSize length, const InnerWordSize nb_bits)) {
  this->nb_bits = nb_bits;
  data.Resize(length*nb_bits);
  return *this;
}

tpl(CountWordSize, GetLength() const) {
  return data.GetLength()/nb_bits;
}

tpl(InnerWordSize, GetBitLength() const) {
  return data.GetLength();
}

tpl(InnerWordSize,  GetMinimalBit(InnerWordSize w)) {
  for(InnerWordSize i = (sizeof(InnerWordSize) * 8)-1 ; i > 0 ; i--) {
    if(w>>i) {
      return i+1;
    }
  }
  return 1;
}


#ifdef TEST_BITVECTOR
// Compile this test program with the following command line:
// g++ bitvector.cxx -DTEST_BITVECTOR -o test_bitvector -Wall -ansi -pedantic -O4

#include <sys/time.h>
#include <sys/resource.h>
#include <vector>


using std::cout;
using std::endl;

template <typename CountWordSize, typename InnerWordSize>
void TestBitVectors(char *cws, char *iws, CountWordSize default_size) {
  cout << "Testing class BitVector<" << cws << ", " << iws << ">" << endl;
  cout << "SizeOf(" << cws << ") = "
       << sizeof(CountWordSize) << "Bytes (= "
       << 8*sizeof(CountWordSize) << "bits)" << endl;
  cout << "SizeOf(" << iws << ") = "
       << sizeof(InnerWordSize) << "Bytes (= "
       << 8*sizeof(InnerWordSize) << "bits)" << endl << endl;
  cout << endl;

  cout << "Testing constructor without arguments" << endl;
  BV *bv = new BV();
  cout << "It gives a vector of size " << (long unsigned int) bv->GetLength() << endl;
  cout << "having " << (long unsigned int) bv->GetBitLength() << " bits" << endl;
  cout << "Destroying the bitvector" << endl;
  delete bv;

  cout << "Testing constructor with first argument = " << (unsigned long int) default_size << endl;
  bv = new BV(default_size);
  cout << "It gives a vector of size " << (long unsigned int) bv->GetLength() << endl;
  cout << "having " << (long unsigned int) bv->GetBitLength() << " bits" << endl;
  cout << "Destroying the bitvector" << endl;
  delete bv;

  cout << "Testing constructor with first argument = " << (unsigned long int) default_size << " and second argument 3" << endl;
  bv = new BV(default_size, 3);
  cout << "It gives a vector of size " << (long unsigned int) bv->GetLength() << endl;
  cout << "having " << (long unsigned int) bv->GetBitLength() << " bits" << endl;

  cout << "Testing copy constructor" << endl;
  BV v(*bv);
  cout << "It gives a vector of size " << (long unsigned int) v.GetLength() << endl;
  cout << "having " << (long unsigned int) v.GetBitLength() << " bits" << endl;

  cout << "Testing assignment operator" << endl;
  *bv = v;
  cout << "It gives a vector of size " << (long unsigned int) bv->GetLength() << endl;
  cout << "having " << (long unsigned int) bv->GetBitLength() << " bits" << endl;
  cout << "Destroying the bitvector" << endl;
  delete bv;

  cout << "Accessing value at position " << (unsigned long int) (default_size - 3) << endl;
  cout << "v[" << (unsigned long int) (default_size - 3) << "] = "
       << (long unsigned int) v[default_size - 3] << endl;
  cout << "v.GetValue(" << (unsigned long int) (default_size - 3) << ") = "
       << (long unsigned int) v.GetValue(default_size - 3) << endl;

  cout << "Setting value 5 at position " << (unsigned long int) (default_size - 3) << endl;
  v.SetValue(default_size - 3, 5);

  cout << "v.GetValue(" << (unsigned long int) (default_size - 3) << ") = "
       << (long unsigned int) v.GetValue(default_size - 3) << endl;

  cout << "Setting value 5 from position 0 to 3" << endl;
  for (CountWordSize i = 0; i < 4; i++) {
    v.SetValue(i, 5);
  }
  for (CountWordSize i = 0; i < 10; i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }

  cout << "Clearing the vector" << endl;
  v.Clear();
  for (CountWordSize i = 0; i < 5; i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }
  cout << "..." << endl;
  for (CountWordSize i = v.GetLength() - 3; i < v.GetLength(); i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }

  cout << "Setting value 5 from position 0 to 3" << endl;
  for (CountWordSize i = 0; i < 4; i++) {
    v.SetValue(i, 5);
  }
  for (CountWordSize i = 0; i < 5; i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }
  cout << "..." << endl;
  for (CountWordSize i = v.GetLength() - 3; i < v.GetLength(); i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }

  cout << "Resize the vector to 3" << endl;
  v.Resize(3);
  cout << "It gives a vector of size " << (long unsigned int) v.GetLength() << endl;
  cout << "having " << (long unsigned int) v.GetBitLength() << " bits" << endl;
  for (CountWordSize i = 0; i < v.GetLength(); i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }

  cout << "Resize the vector to 10" << endl;
  v.Resize(10);
  cout << "It gives a vector of size " << (long unsigned int) v.GetLength();
  cout << " having " << (long unsigned int) v.GetBitLength() << " bits" << endl;
  for (CountWordSize i = 0; i < v.GetLength(); i++) {
    cout << "v.GetValue(" << (long unsigned int) i << ") = " << (long unsigned int) v.GetValue(i) << endl;
  }
  cout << "===" << endl << endl;
  
}

#define TestBV(type1, type2, v) TestBitVectors<type1, type2>((char*) #type1, (char*) #type2, v)

#define GetRUsage(fct, name)						\
  getrusage(RUSAGE_SELF, &usage);					\
  elapsed_utime = usage.ru_utime;					\
  elapsed_stime = usage.ru_stime;					\
  fct;									\
  getrusage(RUSAGE_SELF, &usage);					\
  elapsed_utime.tv_sec = usage.ru_utime.tv_sec - elapsed_utime.tv_sec;	\
  if (usage.ru_utime.tv_usec > elapsed_utime.tv_usec) {			\
    elapsed_utime.tv_usec = usage.ru_utime.tv_usec - elapsed_utime.tv_usec; \
  } else {								\
    elapsed_utime.tv_sec--;						\
    elapsed_utime.tv_usec = 1000000 - elapsed_utime.tv_usec + usage.ru_utime.tv_usec; \
  }									\
  elapsed_stime.tv_sec = usage.ru_stime.tv_sec - elapsed_stime.tv_sec;	\
  if (usage.ru_stime.tv_usec > elapsed_stime.tv_usec) {			\
    elapsed_stime.tv_usec = usage.ru_stime.tv_usec - elapsed_stime.tv_usec; \
  } else {								\
    elapsed_stime.tv_sec--;						\
    elapsed_stime.tv_usec = 1000000 - elapsed_stime.tv_usec + usage.ru_stime.tv_usec; \
  }									\
  cout << "- " #name ": " << endl;					\
  cout << "  - User time: ";						\
  cout << "~ " << (elapsed_utime.tv_sec*1e9 + elapsed_utime.tv_usec*1e3)/nb_tests << " nanosecond by function call" << endl; \
  cout << "  - System time: ";						\
  cout << "~ " << (elapsed_stime.tv_sec*1e9 + elapsed_stime.tv_usec*1e3)/nb_tests << " nanosecond by function call" << endl

int main(int argc, char** argv) {  
  // TestBV(<any>, unsigned char, 2^8/8-1);
  TestBV(unsigned char, unsigned char, 31);
  // TestBV(<unsigned short int or more>, unsigned short int, 2^16/16-1);
  TestBV(unsigned short int, unsigned short int, 4095);
  // TestBV(<unsigned int or more>, unsigned int, 2^32/32-1);
  // TestBV(unsigned long int, unsigned long int, 2^64/64-1);

  // TestBV(unsigned char, <unsigned short int or more>, 2^8-1);
  TestBV(unsigned char, unsigned short int, 255);
  // TestBV(unsigned short int, <unsigned int or more>, 2^16-1);
  TestBV(unsigned short int, unsigned long int, 65535);
  // TestBV(unsigned int, unsigned long int, 2^32-1);

  cout << "Creating a vector of size 1000000000 x 20 bits" << endl;
  BitVector<unsigned int, unsigned long int> bv(1000000000, 20);
  cout << "It gives a vector of size " << (long unsigned int) bv.GetLength() << endl;
  double l = (unsigned long int)bv.GetBitLength(); 
  cout << "having " << (unsigned long int) l << " bits";
  char u = 0;
  l/=8;
  while (l/1024 > 1) {
    l /= 1024;
    u++;
  }
  char unit[] = "?B";
  switch (u) {
  case 0:
    cout << endl;
    return 0;
  case 1:
    unit[0] = 'K';
    break;
  case 2:
    unit[0] = 'M';
    break;
  case 3:
    unit[0] = 'G';
    break;
  case 4:
    unit[0] = 'T';
    break;
  default:
    break;
  }
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(1);
  cout << " (=" << l << unit << ")" << endl;


  cout << "Finally comparing performances of OldSetValue and SetValue " << endl;

  const unsigned int nb_tests = 10000000;
  struct rusage usage;
  struct timeval elapsed_utime, elapsed_stime;
  cout << "Statistics (computed on " << nb_tests << " tests / this may take some time):" << endl;
  unsigned int x;
  cout.precision(3);

  x = 0;
  GetRUsage(
	    while (x++ < nb_tests) {
	      unsigned int i = lrand48() % 1000000000;
	      long unsigned int v = lrand48();
	      bv.OldSetValue(i, v);
	    }, OldSetValue);

  x = 0;
  GetRUsage(
	    while (x++ < nb_tests) {
	      unsigned int i = lrand48() % 1000000000;
	      long unsigned int v = lrand48();
	      bv.SetValue(i, v);
	    }, SetValue);

  x = 0;
  GetRUsage(
	    while (x++ < nb_tests) {
	      unsigned int i = lrand48() % 1000000000;
	      bv.GetValue(i);
	    }, GetValue);

  cout << "Resize the vector to size 20000 x 13 bits" << endl;
  bv.CompleteResize(20000, 13);
  cout << "It gives a vector of size " << (long unsigned int) bv.GetLength() << endl;
  l = (unsigned long int)bv.GetBitLength(); 
  cout << "having " << (unsigned long int) l << " bits";
  u = 0;
  l/=8;
  while (l/1024 > 1) {
    l /= 1024;
    u++;
  }
  unit[0] = '?';
  switch (u) {
  case 0:
    cout << endl;
    return 0;
  case 1:
    unit[0] = 'K';
    break;
  case 2:
    unit[0] = 'M';
    break;
  case 3:
    unit[0] = 'G';
    break;
  case 4:
    unit[0] = 'T';
    break;
  default:
    break;
  }
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(1);
  cout << " (=" << l << unit << ")" << endl;

  //  bv.Clear();
  for (unsigned int i = 0; i < bv.GetLength(); i++) {
    unsigned long int n = rand()%8000;
    bv.OldSetValue(i, n);
    if (bv.GetValue(i) != n) {
      cout << "ERROR at bv.GetValue(" << i << ") = " << bv.GetValue(i) << " <> n = " << n<< endl;
    }
    bv.SetValue(i, n);
    if (bv.GetValue(i) != n) {
      cout << "ERROR at bv.SetValue(" << i << ", " << n <<") = " << bv.GetValue(i) << " <> n = " << n<< endl;
    }
  }

  return 0;
}
#endif
#endif
