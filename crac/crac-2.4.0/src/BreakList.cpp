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
#include "BreakList.h"
#include <cassert>
// #include <iostream>
// using namespace std;

BreakList::BreakList(SupportBreak **breaks, size_t nb_breaks):vector<SupportBreak *>(nb_breaks, NULL) {
  for (size_t i = 0; i < nb_breaks; i++) {
    // cout << breaks[i] << endl;
    (*this)[i] = breaks[i];
  }
  if (nb_breaks)
    parameters = breaks[0]->getParameters();
}

/**
 * @return the length between the breaks of index i and j
 * @pre 0 <= i < size() && 0 <= j < size() && i != j
 */
size_t BreakList::getGapSize(size_t i, size_t j) const{
  assert (0 <= i && i < size() && 0 <= j && j < size() && i != j);
  if (j < i)
    return getGapSize(j, i);

  return (*this)[j]->getPositionStartBreak() - (*this)[i]->getPositionEndBreak() - 1;
}

size_t BreakList::getGapSize(size_t i) const{
  return getGapSize(i, i+1);
}

size_t BreakList::getTotalBreakSize(size_t i, size_t j) const{
  assert (0 <= i && i < size() && 0 <= j && j < size() && i != j);
  if (j < i)
    return getTotalBreakSize(j, i);

  return (*this)[j]->getPositionEndBreak() - (*this)[i]->getPositionStartBreak() + 1;
}

size_t BreakList::getTotalBreakSize(size_t i) const{
  return getTotalBreakSize(i, i+1);
}

bool BreakList::hasShortBreakSpan(size_t i, size_t j) const {
  return getTotalBreakSize(i, j) <= parameters->max_nb_overlapping_breaks * ((*this)[i])->getThreshold();
}

bool BreakList::hasShortGapBetweenBreaks(size_t i, size_t j) const {
  return getGapSize(i, j) <= parameters->max_bases_randomly_matched;
}

void BreakList::merge(size_t i, size_t j) {
  assert (0 <= i && i < size() && 0 <= j && j < size() && i != j);
  if (j < i)
    merge(j, i);
  else {
    SupportBreak *first_break = (*this)[i];
    SupportBreak *second_break = (*this)[j];
    SupportBreak *final_break = new SupportBreak(*first_break, *second_break);

    this->erase(this->begin() + j);
    (*this)[i] = final_break;

    delete first_break;
    delete second_break;
  }
}

void BreakList::merge(size_t i) {
  assert (0 <= i && i < size()-1);
  merge(i, i+1);
}

size_t BreakList::conditional_neighbor_merge(size_t i, std::function <bool(size_t, size_t)> f) {
  if (i > 0 && f(i-1, i)) {
    merge(i-1);
    return i-1;
  } else if (i < size() - 1 && f(i, i+1)) {
    merge(i);
    return i;
  }
  return ~0; 
}

