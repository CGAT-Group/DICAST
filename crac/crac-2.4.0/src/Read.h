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

#include <string.h> // for strlen function
#include "types.h" // for uint type
#include "utils.h" // for majNucNCpy
#include "libReadsReader/readsReader.h" // for readIterator

#ifndef READ_H
#define READ_H

class Read{
  public:
    uint id;
    char *seq;
    char *name;
    char *qual;
    uint length;

    /*
     * Create a Read object
     * @param i id of the read in the file
     * @param s sequence of the read
     * @param n name of the read
     * @param q quality of the read
     * @param l read (desired) length
     */
    Read(uint i, const char*s, const char*n, const char*q, uint l=0) {
      if(l == 0) {
        init(i,s,n,q,strlen(s));
      } else {
        init(i,s,n,q,l);
      }
    }

    /*
     * Destructor
     */
    ~Read() {
      if(seq)
        delete[] seq;
      if(name)
        delete[] name;
      if(qual)
        delete[] qual;
    }

    /**
     * @return the length of the read
     */
    size_t getLength() { return length; }

  private:

    /*
     * Generic method to init the read object.
     * This method is (and should be) called by each constructor
     */
    void init(uint i, const char*s, const char*n, const char*q, uint l) {
      length = l;
      if(s) {
        seq = new char[l+1];
        majNucNCpy(seq,s,l);
        seq[l] = '\0';
      } else {
        seq = NULL;
      }
      if(n) {
        name = new char[strlen(n)+1];
        strcpy(name,n);
      } else {
        name = NULL;
      }
      if(q) {
        qual = new char[l+1];
        strncpy(qual,q,l);
        qual[l] = '\0';
      } else {
        qual = NULL;
      }
      id = i;
    }
};

#endif //READ_H
