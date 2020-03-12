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

#include "HuffAlphabetRank.h"
#include <cstdlib>

   THuffAlphabetRank::THuffAlphabetRank(uchar *s, ulong n, TCodeEntry *codetable, ulong level, ulong factor) {
      left = NULL;
      right = NULL;
      bitrank = NULL;
//       parent = NULL;
      ch = s[0];
      this->codetable = codetable;
      ulong *Binbits = (ulong *)calloc(n/W+1, sizeof(ulong));;
      ulong sum=0,i;

      for (i=0;i<n;i++)
         if ((codetable[s[i]].code & (1u << level)) )  {
           SetField(Binbits,1,i,true);
           sum++;
         } else 
           SetField(Binbits,1,i,false);
      leaf=true;
      for (i=1;i<n;i++)
        if (s[i] != ch)
          leaf = false;
      if (leaf) {
        free(Binbits);
        delete [] s;
        return;
      }
      
      uchar *sfirst=NULL, *ssecond=NULL;
      if (n-sum > 0) {
        sfirst = new uchar[n-sum];
      }
      if (sum > 0) {
        ssecond = new uchar[sum];
      }
      ulong j=0,k=0;
      for (i=0;i<n;i++)
        if (GetField(Binbits,1,i)) ssecond[k++] = s[i];
        else sfirst[j++] = s[i];

      delete [] s;
      bitrank = new BitRankW32Int(Binbits,n,true,factor);
      if (j > 0) {
        left = new THuffAlphabetRank(sfirst,j,codetable,level+1,factor);
      } else if (sfirst) {
        delete [] sfirst;
      }
      if (k>0) {
        right = new THuffAlphabetRank(ssecond,k,codetable,level+1,factor);
      } else if (ssecond) {
        delete [] ssecond;
      }
//       if (level == 0)
//         init_leaves();
   }

// void THuffAlphabetRank::init_leaves() {
//   leaves = new THuffAlphabetRank[size_uchar];
//   for (uchar c = 0; c < size_uchar; c++) {
//     leaves[c] = getLeaf(c);
//   }
// }

THuffAlphabetRank *THuffAlphabetRank::getLeaf(uchar c) {
  if (codetable[c].count == 0) return NULL;
  ulong level = 0;
  ulong code = codetable[c].code;
  THuffAlphabetRank *currentNode = this;
  while (! currentNode->leaf) {
    if ((code & (1u<<level)) == 0) {
      currentNode = currentNode->left;
    } else {
      currentNode = currentNode->right;
    }
    ++level;
  }
  return currentNode;
}

   ulong THuffAlphabetRank::rank(uchar c, ulong i) { // returns the number of characters c before and including position i
      THuffAlphabetRank *temp=this;
      if (codetable[c].count == 0) return 0;
      ulong level = 0;
      ulong code = codetable[c].code;
      while (!temp->leaf) {
        if ((code & (1u<<level)) == 0) {
          i = i-temp->bitrank->rank(i);
          temp = temp->left;
        }
        else {
          i = temp->bitrank->rank(i)-1;
          temp = temp->right;
        }
        ++level;
      }
      return i+1;
   }

// ulong THuffAlphabetRank::select(uchar c, ulong i) {
//   THuffAlphabetRank *node = leaves[c];

//   if (node->parent == NULL)
//     return i;
  
//   bool bit = node->parent->right == node;

//   while(node->parent != NULL) {
//     node = node->parent;
//     node->bitrank->select
//   }
// }


   bool THuffAlphabetRank::IsCharAtPos(uchar c, ulong i) {
      THuffAlphabetRank *temp=this;
      if (codetable[c].count == 0) return false;
      ulong level = 0;
      ulong code = codetable[c].code;
      while (!temp->leaf) {
         if ((code & (1u<<level))==0) {
            if (temp->bitrank->IsBitSet(i)) return false;
            i = i-temp->bitrank->rank(i);
            temp = temp->left;
         }
         else {
            if (!temp->bitrank->IsBitSet(i)) return false;
            i = temp->bitrank->rank(i)-1;
            temp = temp->right;
         }
         ++level;
      }
      return true;
   }
   uchar THuffAlphabetRank::charAtPos(ulong i) {
      THuffAlphabetRank *temp=this;
      while (!temp->leaf) {
        if (temp->bitrank->IsBitSet(i)) {
          i = temp->bitrank->rank(i)-1;
          temp = temp->right;
        } else {
          i = i-temp->bitrank->rank(i);
          temp = temp->left;
        }
      }
      return temp->ch;
   }
   uchar THuffAlphabetRank::charAtPos2(ulong i, ulong *rank) {
     THuffAlphabetRank *temp=this;
     while (!temp->leaf) {
       if (temp->bitrank->IsBitSet(i)) {
         i = temp->bitrank->rank(i)-1;
         temp = temp->right;
       } else {
         i = i-temp->bitrank->rank(i);
         temp = temp->left;
       }
     }
     (*rank)=i;
     return temp->ch;
   }

   ulong THuffAlphabetRank::SpaceRequirementInBits() {
      ulong bits=sizeof(THuffAlphabetRank)*8;
      if (left!=NULL) bits+=left->SpaceRequirementInBits();
      if (right!=NULL) bits+=right->SpaceRequirementInBits();
      if (bitrank != NULL)
         bits += bitrank->SpaceRequirementInBits();
      return bits;
   }
   bool THuffAlphabetRank::Test(uchar *s, ulong n) {
      // testing that the code works correctly
      uchar C[size_uchar];
      ulong i,j;
      bool correct=true;
      for (j=0;j<size_uchar;j++)
         C[j] = 0;
      for (i=0;i<n;i++) {
         C[s[i]]++;
         if (C[s[i]]!=rank(s[i],i)) {
            correct = false;
            printf("%u (%c): %d<>%u\n",i,s[i],C[s[i]],rank(s[i],i));
         }
      }
      return correct;
   }
   THuffAlphabetRank::~THuffAlphabetRank() {
      if (left!=NULL) delete left;
      if (right!=NULL) delete right;
      if (bitrank!=NULL)
         delete bitrank;
   }

   int THuffAlphabetRank::save(FILE *f) {
     uchar sons;
     if (f == NULL) return 20;
     if (fwrite (&leaf,sizeof(uchar),1,f) != 1) return 21;
     if (fwrite (&ch,sizeof(uchar),1,f) != 1) return 21;
     if (!leaf) {
       if (bitrank->save(f) != 0) return 21;
       sons=0;
       if (left != NULL) sons=sons | 1;
       if (right != NULL) sons=sons | 2;
       if (fwrite (&sons,sizeof(uchar),1,f) != 1) return 21;
       if (left != NULL) 
         if (left->save(f)!=0) return 21;
       if (right != NULL)
         if (right->save(f)!=0) return 21;
     }
     return 0;
   }

   int THuffAlphabetRank::load(FILE *f, TCodeEntry *_codetable) {
     int error;
     uchar sons;
     left = NULL;
     right = NULL;
//      parent = NULL;
     bitrank = NULL;

     if (f == NULL) return 23;
     if (fread (&leaf,sizeof(uchar),1,f) != 1) return 25;
     if (fread (&ch,sizeof(uchar),1,f) != 1) return 25;
     codetable=_codetable;
     if (!leaf){ 
       bitrank= new BitRankW32Int(f,&error); 
       if (error != 0) return error;
       if (fread (&sons,sizeof(uchar),1,f) != 1) return 25;
       if ((sons&1) > 0 ) {
         left = new THuffAlphabetRank(f, _codetable, &error);
//          left->parent = this;
         if (error != 0) return error;
       }
       if ((sons&2) > 0 ) {
         right = new THuffAlphabetRank(f, _codetable, &error);
//          right->parent = this;
         if (error != 0) return error;
       }
     }
     return 0;
   }
   THuffAlphabetRank::THuffAlphabetRank(FILE *f, TCodeEntry *_codetable, int *error) {
     *error = load(f, _codetable);
//      init_leaves();
   }
