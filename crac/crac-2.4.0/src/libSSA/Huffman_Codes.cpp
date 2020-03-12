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

#include "Huffman_Codes.h"
//
// This routine does a quick count of all the
// characters in the input text.  
//


void count_chars(uchar *text, ulong n, TCodeEntry *counts )
{
    ulong i;
    for (i = 0 ; i < size_uchar ; i++ )
        counts[ i ].count = 0;
    for (i=0; i<n; i++)
        counts[text[i]].count++; 
}

TCodeEntry *makecodetable(uchar *text, ulong n)
{
    TCodeEntry *result = new TCodeEntry[ size_uchar ];
    
    count_chars( text, n, result );
/*    ulong count=0;
    for (int i=0;i<256; i++) {
       count += result[i].count;
       if (result[i].count>0)
          printf("%c,%d: %d\n",i , i,result[i].count);
    }
    printf("count=%d\n",count);*/
    priority_queue< node, vector< node >, greater<node> > q;
//
// First I push all the leaf nodes into the queue
//
    for ( ulong i = 0 ; i < size_uchar ; i++ )
        if ( result[ i ].count )
            q.push(node( i, result[ i ].count ) );
//
// This loop removes the two smallest nodes from the
// queue.  It creates a new internal node that has
// those two nodes as children. The new internal node
// is then inserted into the priority queue.  When there
// is only one node in the priority queue, the tree
// is complete.
//

    while ( q.size() > 1 ) {
        node *child0 = new node( q.top() );
        q.pop();
        node *child1 = new node( q.top() );
        q.pop();
        q.push( node( child0, child1 ) );
    }
//
// Now I compute and return the codetable
//
//     cout << "Char  Symbol   Code" << endl;
//     q.top().traverse();
    q.top().maketable(0u,0u, result);
    q.top().deleteme();
    q.pop();
    return result;
}

int save_codetable(FILE *f,TCodeEntry *codetable) {
  if (f == NULL) return 20;
  for (int i=0; i<size_uchar;i++) {
    if (codetable[i].save(f) != 0) return 21;
  }
  return 0;
}


int load_codetable(FILE *f,TCodeEntry **codetable) {
  if (f == NULL) return 23;
   TCodeEntry *result = new TCodeEntry[ size_uchar ];
  for (int i=0; i<size_uchar;i++) {
    if (result[i].load(f) != 0) return 25;
  }
  (*codetable)=result;
  return 0;
}

