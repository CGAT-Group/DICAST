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
*                            Technologie).	                              *
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

#ifndef bitarray_h
#define bitarray_h

#include "basic.h"

class BitRankF {
private:
  ulong *data; //here is the bit-array
  bool owner;
  ulong n,integers;
  ulong factor,b,s;
  ulong *Rs; //superblock array
  uchar *Rb; //block array
  ulong BuildRankSub(ulong ini,ulong fin); //internal use of BuildRank
  void BuildRank(); //crea indice para rank
public:
  BitRankF(ulong *bitarray, ulong n, bool owner);
  ~BitRankF(); //destructor
  bool IsBitSet(ulong i);
  ulong rank(ulong i); //Rank from 0 to n-1
  ulong rank2(ulong i); //Rank from 0 to n-1
  ulong prev(ulong start); // gives the largest index i<=start such that IsBitSet(i)=true
  ulong select(ulong x); // gives the position of the x:th 1.
  ulong SpaceRequirementInBits();
  /*load-save functions*/
  int save(FILE *f);
  int load(FILE *f);
  BitRankF(FILE *f, int *error); 
};


class BitSelectNext {
private:
	ulong *datos; //arreglo de bits
	bool owner;
	ulong n;
	ulong integers;
public:
	// Crea arreglo segun numero de bits, semilla aleatoria y probabilidad
	BitSelectNext(ulong *bit, ulong n, bool owner); 
	~BitSelectNext(); //destructor
	ulong select_next(ulong i); // select_next
};

class BitRankFSparse {
private:
  BitRankF *block;
  BitRankF *sblock;
  ulong L;
public:
  BitRankFSparse(ulong *bitarray, ulong n);
  ~BitRankFSparse(); //destructor
  bool IsBitSet(ulong i);
  ulong rank(ulong i); //Rank from 0 to n-1
  ulong prev(ulong start); // gives the largest index i<=start such that IsBitSet(i)=true
  ulong select(ulong x); // gives the position of the x:th 1.
  ulong SpaceRequirementInBits();
};

#endif
