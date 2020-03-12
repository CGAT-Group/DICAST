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

#ifndef CIGAR_H
#define CIGAR_H
#include <iostream>
#include <vector>
#include "cigarTypes.h"

class Cigar {
 private:
  /* Storing the whole CIGAR string */
  vector<cigar_type> cigar;

  /* Use X in cigar */
  bool use_x_in_cigar;

  /* Counters for each operator */
  uint nb_M, nb_I, nb_D, nb_S, nb_P, nb_H, nb_EQ, nb_X, nb_N, nb_chimeras;

  void update_counters(const cigar_type &cigar);

 public:
  
  Cigar(bool use_x_in_cigar = false);

  /* Command */
  /**
   * Append a new element in the CIGAR string
   * @param nb: the number related to the CIGAR operator (eg. 10 in 10M)
   * @param operator: the operator to be used 
   * @pre op in (M, I, D, S, P, H, =, X, N)
   * @post count() == old count() + 1
   *  AND get(count()-1) == {nb, op}
   */
  void append(uint nb, char op);

  /**
   * Append a new element in the CIGAR string
   * @param cigar: the CIGAR operator (with the number)
   * @post count() == old count() + 1
   *  AND get(count()-1) == cigar
   */
  void append(cigar_type cigar);

  /**
   * Filter consecutive operators with the same type, remove
   * operators whose length is 0 (apart from N).
   * @post count() <= old count()
   *   AND get(i).type != get(i+1).type (when get(i).type != 'N') for i = 0 ... count()-2 
   *   AND get(i).nb != 0 (when get(i).type != 'N') for i = 0 ... count()-1
   */
  void filter();

  /**
   * Prepend a new element in the CIGAR string
   * @param nb: the number related to the CIGAR operator (eg. 10 in 10M)
   * @param operator: the operator to be used 
   * @pre op in (M, I, D, S, P, H, =, X, N)
   * @post count() == old count() + 1
   *  AND get(0) == {nb, op}
   */
  void prepend(uint nb, char op);

  /**
   * Prepend a new element in the CIGAR string
   * @param cigar: the CIGAR operator (with the number)
   * @post count() == old count() - 1
   *  AND get(0) == cigar
   */
  void prepend(cigar_type cigar);

  /* Queries */

  /**
   * @return Number of operators in the string
   */
  uint count() const;

  /**
   * @param i: Element from the CIGAR string to retrieve
   * @pre 0 <= i < count()
   * @return the i-th element in CIGAR.
   */
  const cigar_type& get(uint i) const;

  /**
   * @return the length of the alignement on the reference 
   *        (ie sume of the numbers associated with a S, M, I, = or X operator).
   */
  uint getQueryAlignementLength() const;

  /**
   * @return the length of the alignement on the reference 
   *        (ie sume of the numbers associated with a M, D, N, = or X operator).
   */
  uint getReferenceAlignementLength() const;

  /**
   * @return true is the operator assigned to the ith element correspond to 
   * the reference sequence (ie. M, =, X, N, D)
   */
  bool isReferenceBasedElement(uint i) const;

  /**
   * @return true is the operator assigned to the ith element correspond to 
   * the read sequence (ie. M, =, X, I, S, H)
   */
  bool isQueryBasedElement(uint i) const;

  /*
   * @retrun the number of deleted nucleotids from the reference.
   *        (ie sum of the numbers associated with a D or N operator).
   */
  uint getNbDeletions() const;

  /**
   * @return the number of positions that are in the read
   *         (ie sum of the numbers associated with a M, I, S, = or X operator).
   */
  uint getNbPositionsRead() const;

  /**
   * @return the number of positions that are in the read
   *         (ie sum of the numbers associated with a X, I, D, or N operator).
   */
  uint getEditDistance() const;

  /**
   * @return the number of chimeras in the Cigar chain
   *        (ie sum of '0N' elements)
   */
  uint getNbChimeras() const;

  /**
   * Reverse the cigar elements.
   * Example : 8S5M1D3M => 3M1D5M8S
   */
  void reverse();

  cigar_type &operator[](uint i);
};

ostream &operator<<(ostream &, const Cigar &);


#endif
