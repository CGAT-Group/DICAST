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

#ifndef BASICSINCLUDED
#define BASICSINCLUDED

#include <sys/types.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define mask31 0x0000001F

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))




/*numero de bits del entero de la maquina*/
#define W 32
/* W-1 */
#define Wminusone 31
/*numero de bits del entero de la maquina*/
#define WW 64
/*bits para hacer la mascara para contar mas rapido*/
#define bitsM 8
/*bytes que hacen una palabra */
#define BW 4
#define size_uchar 256


/* bits needed to represent a number between 0 and n */
inline uint bits (uint n){
  uint b = 0;
  while (n) { b++; n >>= 1; }
  return b;
}

/* reads bit p from e */
#define bitget(e,p) ((((e)[(p)/W] >> ((p)%W))) & 1)
/* sets bit p in e */
#define bitset(e,p) ((e)[(p)/W] |= (1<<((p)%W)))
/* cleans bit p in e */
#define bitclean(e,p) ((e)[(p)/W] &= ~(1<<((p)%W)))

/* numero de enteros necesarios para representar e elementos de largo n */
#define enteros(e,n) ((e)*(n))/W+(((e)*(n))%W > 0)

inline ulong GetField(ulong *A, register  ulong len, register ulong index) {
  register ulong i=index*len/W, j=index*len-W*i, result;
  if (j+len <= W)
    result = (A[i] << (W-j-len)) >> (W-len);
  else {
    result = A[i] >> j;
    result = result | (A[i+1] << (WW-j-len)) >> (W-len);
  }
  return result;
}

inline void SetField(ulong *A,register ulong len, register ulong index,register  ulong x) {
   ulong i=index*len/W, j=index*len-i*W;
   ulong mask = ((j+len) < W ? ~0u << (j+len) : 0) | ((W-j) < W ? ~0u >> (W-j) : 0);
   A[i] = (A[i] & mask) | x << j;
   if (j+len>W) {
      mask = ((~0u) << (len+j-W));
      A[i+1] = (A[i+1] & mask)| x >> (W-j);
   }
}

inline unsigned GetFieldW32(unsigned *A,register unsigned index) {
  return A[index];
}

inline void SetField32(unsigned *A, register unsigned index,register unsigned x) {
  A[index]=x;
}

inline unsigned GetFieldW16(unsigned *A,register unsigned index) {
  register unsigned i=index/2, j=(index&1)<<4, result;
  result = (A[i] << (16-j)) >> (16);
  return result;
}

inline unsigned GetFieldW4(unsigned *A,register unsigned index) {
  register unsigned i=index/8, j=(index&0x7)<<2;
  /*register unsigned i=index/8, j=index*4-32*i; */
  return (A[i] << (28-j)) >> (28);
}



#endif

