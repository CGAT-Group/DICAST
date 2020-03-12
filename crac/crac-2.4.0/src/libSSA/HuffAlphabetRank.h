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

#ifndef HuffAlphabetRank_h
#define HuffAlphabetRank_h

#include "Huffman_Codes.h"
#include "basic.h"
#include "bitrankw32int.h"

class THuffAlphabetRank {
// using fixed 0...255 alphabet
private:
   TCodeEntry *codetable;
   BitRankW32Int *bitrank;
   uchar ch;
   bool leaf;
   THuffAlphabetRank *left;
   THuffAlphabetRank *right;
/*    THuffAlphabetRank *parent; */

/*    THuffAlphabetRank *leaves; */
/*    void init_leaves(); */
public:
   THuffAlphabetRank(uchar *s, ulong n, TCodeEntry *codetable, ulong level, ulong factor);
   THuffAlphabetRank *getLeaf(uchar c);
   ulong rank(uchar c, ulong i); // returns the number of characters c before and including position i
   bool IsCharAtPos(uchar c, ulong i);
   uchar charAtPos(ulong i);
   uchar charAtPos2(ulong i, ulong *rank);
   ulong SpaceRequirementInBits();
   bool Test(uchar *s, ulong n);
   ~THuffAlphabetRank();
   int save(FILE *f);
   int load(FILE *f, TCodeEntry *_codetable);
   THuffAlphabetRank(FILE *f, TCodeEntry *_codetable, int *error);
};
#endif
