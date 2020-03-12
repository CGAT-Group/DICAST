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

#include "bitrankw32int.h"
#include "assert.h"
#include "math.h"
#include <sys/types.h>
#include <cstdlib>

/////////////
//Rank(B,i)// 
/////////////
//_factor = 0  => s=W*lgn
//_factor = P  => s=W*P
//Is interesting to notice
//factor=2 => overhead 50%
//factor=3 => overhead 33%
//factor=4 => overhead 25%
//factor=20=> overhead 5%

BitRankW32Int::BitRankW32Int( ulong *bitarray, ulong _n, bool owner, ulong _factor){
  data=bitarray;
  //cdata = (ulong *)calloc(_n/W+1, sizeof(ulong));  
  //ulong tmp = 0;
  //for(ulong i=0; i<(_n/W+1); i++)  {
  //  tmp += popcount32(data[i]);
  //  cdata[i] = tmp;
  //}
  this->owner = owner;
  this->as_cache = false;
  this->n=_n;
	ulong lgn=bits(n-1);
	this->factor=_factor;
	if (_factor==0) this->factor=lgn;
	else this->factor=_factor;
	b=32;
	s=b*this->factor;
  integers = n/W+1;
 	BuildRank();
}

BitRankW32Int::~BitRankW32Int() {
	delete [] Rs;
	if (owner) free(data);
	if (as_cache) free(cache_rank);
}

//Metodo que realiza la busqueda d
void BitRankW32Int::BuildRank(){
	ulong num_sblock = n/s;
	Rs = new ulong[num_sblock+1];// +1 pues sumo la pos cero
	for(ulong i=0;i<num_sblock+1;i++)
	  Rs[i]=0;
	ulong j;
	Rs[0]=0;
	for (j=1;j<=num_sblock;j++) {
		Rs[j]=Rs[j-1];
		Rs[j]+=BuildRankSub((j-1)*factor,factor);
	}
	BuildCacheRank();
}

ulong BitRankW32Int::BuildRankSub(ulong ini,ulong bloques){
	ulong rank=0,aux;
	for(ulong i=ini;i<ini+bloques;i++) {
		if (i < integers) {
			aux=data[i];
			rank+=popcount32(aux);
		}
	}
	return rank; //retorna el numero de 1's del intervalo

}

/*
 * Building cache rank array
 * Block of 8x32 set of data
 * [ RS ] [ sum/sum/sum sum/sum/sum ] [ data0 data1 data2 data3 data4 data5 data6 data7 ]
 * RS => contain total sum  ( 32bit )
 * sum => contain local sum ( 8bit  )
 * dataX => contain data n  ( 32bit )
 */
void BitRankW32Int::BuildCacheRank() {
  if ( this->as_cache == true ) return;
  ulong num_sblock = n/(8 * 32);
  num_sblock++;
  as_cache = true;
  cache_rank = (s_cache_rank *)calloc(num_sblock, sizeof(s_cache_rank));
  ulong nb_data = n/32;
  ulong pc, i, j, k;
  unsigned char lcount, tpc;
  pc = 0;
  k = 0;
  for (j=0; j<=nb_data; j+=8) {
    cache_rank[k].total_sum = pc;                          // Sum of all popcount32 seen before
    lcount = 0;
    for( i=j; i<(j+8); i++ ) {
      if ( i <= nb_data ) {
        cache_rank[k].data[i-j] = data[i];
        tpc = popcount32( data[i] );
        pc += tpc;
        cache_rank[k].local_sum[i-j] = lcount;
        lcount += tpc;
      }
    }
    k++;
   }
}


ulong BitRankW32Int::rank(ulong i) {
  if ( i > n ) return 0;
  ulong rang, cas, somm;
  unsigned short rest_i, rest;
  cas = (i / 32);
  rest_i = 31 - (i % 32);
  rang = cas / 8;
  rest = cas % 8;
  somm = cache_rank[rang].total_sum + cache_rank[rang].local_sum[rest];
  somm += popcount32(cache_rank[rang].data[rest] << rest_i );
  return somm;
}

bool BitRankW32Int::IsBitSet(ulong i) {
  return (1u << (i % W)) & data[i/W];
}

int BitRankW32Int::save(FILE *f) {
  if (f == NULL) return 20;
  if (fwrite (&n,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&factor,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (data,sizeof(ulong),n/W+1,f) != n/W+1) return 21;
  if (fwrite (Rs,sizeof(ulong),n/s+1,f) != n/s+1) return 21;
  return 0;
}

int BitRankW32Int::load(FILE *f) {
  if (f == NULL) return 23;
  if (fread (&n,sizeof(ulong),1,f) != 1) return 25;
  b=32; // b is a word
  if (fread (&factor,sizeof(ulong),1,f) != 1) return 25;
  s=b*factor;
  ulong aux=(n+1)%W;
  if (aux != 0)
    integers = (n+1)/W+1;
  else
    integers = (n+1)/W;
  data= (ulong *)calloc(n/W+1, sizeof(ulong));
  if (!data) return 1;
  if (fread (data,sizeof(ulong),n/W+1,f) != n/W+1) return 25;
  this->owner = true;
  this->as_cache = false;
  BuildCacheRank();
  Rs= new ulong[n/s+1];
  if (!Rs) return 1;
  if (fread (Rs,sizeof(ulong),n/s+1,f) != n/s+1) return 25;
  return 0;
}

BitRankW32Int::BitRankW32Int(FILE *f, int *error) {
  *error = BitRankW32Int::load(f);
}

ulong BitRankW32Int::SpaceRequirementInBits() {
  return (owner?n:0)+(n/s)*sizeof(ulong)*8 +sizeof(BitRankW32Int)*8; 
}

ulong BitRankW32Int::prev(ulong start) {
      // returns the position of the previous 1 bit before and including start.
      // tuned to 32 bit machine

      ulong i = start >> 5;
      int offset = (start % W);
      ulong answer = start;
      ulong val = data[i] << (Wminusone-offset);

      if (!val) { val = data[--i]; answer -= 1+offset; }

      while (!val) { val = data[--i]; answer -= W; }

      if (!(val & 0xFFFF0000)) { val <<= 16; answer -= 16; }
      if (!(val & 0xFF000000)) { val <<= 8; answer -= 8; }

      while (!(val & 0x80000000)) { val <<= 1; answer--; }
      return answer;
}

ulong BitRankW32Int::select(ulong x) {
  // returns i such that x=rank(i) && rank(i-1)<x or n if that i not exist
  // first binary search over first level rank structure
  // then sequential search using popcount over a int
  // then sequential search using popcount over a char
  // then sequential search bit a bit

  //binary search over first level rank structure
  ulong l=0, r=n/s;
  ulong mid=(l+r)/2;
  ulong rankmid = Rs[mid];
  while (l<=r) {
    if (rankmid<x)
      l = mid+1;
    else
      r = mid-1;
    mid = (l+r)/2;
    rankmid = Rs[mid];
  }
  //sequential search using popcount over a int
  ulong left;
  left=mid*factor;
  x-=rankmid;
        ulong j=data[left];
        ulong ones = popcount32(j);
        while (ones < x) {
    x-=ones;left++;
    if (left > integers) return n;
          j = data[left];
      ones = popcount32(j);
        }
  //sequential search using popcount over a char
  left=left*b;
  rankmid = popcount8(j);
  if (rankmid < x) {
    j=j>>8;
    x-=rankmid;
    left+=8;
    rankmid = popcount8(j);
    if (rankmid < x) {
      j=j>>8;
      x-=rankmid;
      left+=8;
      rankmid = popcount8(j);
      if (rankmid < x) {
        j=j>>8;
        x-=rankmid;
        left+=8;
      }
    }
  }

  // then sequential search bit a bit
        while (x>0) {
    if  (j&1) x--;
    j=j>>1;
    left++;
  }
  return left-1;
}
