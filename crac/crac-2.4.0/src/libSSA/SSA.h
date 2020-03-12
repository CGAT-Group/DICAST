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

#ifndef SSA_H
#define SSA_H
#include "HuffAlphabetRank.h"
#include "basic.h"
#include "bitarray.h"
#include "Huffman_Codes.h"

class TFMindex {
  private:
    ulong n;
    ulong samplerate;
    ulong C[size_uchar+1];
    TCodeEntry *codetable;
    THuffAlphabetRank *alphabetrank;
    BitRankW32Int *sampled;
    ulong *suffixes;
    ulong *positions;
    uchar map0;
    int remap0(uchar *text, ulong n) {
      int i;
      ulong j;
      ulong Freq_old[size_uchar];
      for(i=0;i<size_uchar;i++)
        Freq_old[i]=0;
      for(j=0;j<n;j++) {
        Freq_old[text[j]]++;
      }
      i=-1;
      // remap alphabet
      if (Freq_old[0]==0) return 0; //test if some character of T is zero
      for(j=1;j<size_uchar;j++)
        if(Freq_old[j]==0) {
          i=j;
          break;
        }
      if (i == -1 ) return i;
      // remap text
      for(j=0;j<n;j++)
        if (text[j] == 0) text[j] = i;
      return i;
    }
    uchar *pattern0(uchar *pattern, ulong m) {
      uchar *pat0;
      pat0 = (uchar *) malloc (sizeof(uchar)*m);
      for (ulong j=0;j<m;j++) {
        if (pattern[j] == 0) pat0[j]=map0;
        else pat0[j]=pattern[j];
      }
      return pat0;
    }
public:

    TFMindex(uchar *bwt, ulong length, ulong samplerate, bool free_text,ulong factor, int *error, uchar terminator=0);

    ulong count(uchar *pattern, ulong m);
    ulong countReverseComp(uchar *pattern, ulong m);

    void getRange(uchar *pattern, ulong m, ulong &sp, ulong &ep);
    void getReverseCompRange(uchar *pattern, ulong m, ulong &sp, ulong &ep);

    void iterCount(uchar c, ulong &sp, ulong &ep);
    void iterCountReverseComp(uchar c, ulong &sp, ulong &ep);

    ulong returnLocate(ulong sp, ulong ep, ulong nb_displayed, ulong **occ, long m=0);
  
    ulong locate(uchar *pattern, ulong m, ulong **occ, ulong nb_displayed=0);
    ulong locateReverseComp(uchar *pattern, ulong m, ulong **occ,ulong nb_displayed=0);

    int display(uchar *pattern, ulong m, ulong numc, ulong *numocc, uchar **snippet_text, ulong **snippet_lengths, ulong **occ);

    ulong get_length ();

    ~TFMindex() ;

    ulong SpaceRequirement();

    ulong SpaceRequirement_locate();

    ulong SpaceRequirement_count();

    int extract(ulong from, ulong to, uchar **snippet);

    int save(char *filename);

    int load(char *filename);

    TFMindex(char *filename, int *error);
};


int build_index(uchar *text, ulong length, char *build_options, void **index);
int save_index(void *index, char *filename);
int load_index(char *filename, void **index);
int free_index(void *index);
int index_size(void *index, ulong *size);
int index_size_count(void *index, ulong *size);
int index_size_locate(void *index, ulong *size);
////////////////////////
////Querying the Index//
////////////////////////
int count(void *index, uchar *pattern, ulong length, ulong *numocc);
int locate(void *index, uchar *pattern, ulong length, ulong **occ, ulong *numocc, ulong nb_displayed=0);
/////////////////////////
//Accessing the indexed//
/////////////////////////
int extract(void *index, ulong from, ulong to, uchar **snippet, ulong *snippet_length);
int display(void *index, uchar *pattern, ulong length, ulong numc, ulong *numocc, uchar **snippet_text, ulong **snippet_lengths);

int get_length(void *index, ulong *length);
////////////////////
////Error handling//
////////////////////
const char *error_index(int e);

#endif
