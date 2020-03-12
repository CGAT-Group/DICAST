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

#ifndef maketable_h
#define maketable_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stack>
#include <string>
#include <queue>
#include "basic.h"

using namespace std;


// To append codes

//ulong SetBit(ulong x, ulong pos, ulong bit) {
//      return x | (bit << pos);
//}

// codetable := new TCodeEntry[256]

class TCodeEntry {
  public:
    ulong count;
    ulong bits;
    ulong code;
    TCodeEntry() {count=0;bits=0;code=0u;};
    int save(FILE *f) {
      if (f == NULL) return 20;
      if (fwrite (&count,sizeof(ulong),1,f) != 1) return 21;
      if (fwrite (&bits,sizeof(ulong),1,f) != 1) return 21;
      if (fwrite (&code,sizeof(ulong),1,f) != 1) return 21;
      return 0;
    };
    int load(FILE *f) {
      if (f == NULL) return 23;
      if (fread (&count,sizeof(ulong),1,f) != 1) return 25;
      if (fread (&bits,sizeof(ulong),1,f) != 1) return 25;
      if (fread (&code,sizeof(ulong),1,f) != 1) return 25;
      return 0;
    };
};

//
// The node class is used to represent both leaf
// and internal nodes. leaf nodes have 0s in the
// child pointers, and their value member corresponds
// to the character they encode.  internal nodes
// don't have anything meaningful in their value
// member, but their child pointers point to
// other nodes.
//

class node {
private:
    int weight;
    uchar value;
    const node *child0;
    const node *child1;
//
// Construct a new leaf node for character c
//
public:
    node( uchar c = 0, int i = -1 ) {
        value = c;
        weight = i;
        child0 = 0;
        child1 = 0;
    }
//
// Construct a new internal node that has
// children c1 and c2.
//
    node( const node* c0,  const node *c1 ) {
        value = 0;
        weight = c0->weight + c1->weight;
        child0 = c0;
        child1 = c1;
    }

//
// The comparison operators used to order the
// priority queue.
//
    bool operator>( const node &a ) const {
        return weight > a.weight;
    }
//
// The traverse member function is used to
// print out the code for a given node.  It
// is designed to print the entire tree if
// called for the root node.
//

    void traverse( string code = "" )  const
    {
        if ( child0 ) {
            child0->traverse( code + "0" );
            child1->traverse( code + "1" );
        } else {
            cout << " " << value << " " << (int) value << "      ";
            cout << setw( 2 ) << weight;
            cout << "     " << code << endl;
        }
    }
    void maketable( ulong code, ulong bits, TCodeEntry *codetable ) const
    {
        if ( child0 ) {
            child0->maketable( code | (0 << bits) , bits+1, codetable );
            child1->maketable( code | (1 << bits) , bits+1, codetable );
        } 
        else {
            codetable[value].code = code;    
            codetable[value].bits = bits;
        }
    }
    void deleteme() const {
        if ( child0 ) {
            if (child0->child0)
                child0->deleteme();
            if (child1->child0)
                child1->deleteme();
            delete child0;
            delete child1;
	}
    }
};

TCodeEntry *makecodetable(uchar *text, ulong n);
int save_codetable(FILE *f,TCodeEntry *codetable);
int load_codetable(FILE *f,TCodeEntry **codetable);
#endif
