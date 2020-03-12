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
#include "SupportBreak.h"
#include <vector>
#include <functional>
#include "Parameters.h"

using namespace std;

class BreakList: public vector<SupportBreak *> {
private:
  Parameters *parameters;
public:

  BreakList(SupportBreak **breaks, size_t nb_breaks);

  /**
   * @return the length between the breaks of index i and j
   * @pre 0 <= i < size() && 0 <= j < size() && i != j
   */
  size_t getGapSize(size_t i, size_t j) const;

  /**
   * @return the length between the breaks of index i and i+1
   * @pre 0 <= i < size()-1 
   */
  size_t getGapSize(size_t i) const;

  /**
   * @return the length from the start of the break of index i to the end of the break of index j
   * @pre 0 <= i < size() && 0 <= j < size() && i != j
   */
  size_t getTotalBreakSize(size_t i, size_t j) const;

  /**
   * @return the length from the start of break i to the end of break i+1
   * @pre 0 <= i < size()-1 
   */
  size_t getTotalBreakSize(size_t i) const;

  /**
   * @return true iff the distance from the start of the leftmost break
   * to the end of the rightmost break can be explained "by chance"
   * (the chance, here, depends on the threshold and the admissible
   * number of overlapping breaks).
   */
  bool hasShortBreakSpan(size_t i, size_t j) const;

  /**
   * @return true iff the distance between the breaks can be explained by
   * ``chance''.
   * The chance here, corresponds to the number of bases randomly matched consecutively.
   */
  bool hasShortGapBetweenBreaks(size_t i, size_t j) const;

  /**
   * Merge breaks i and j in break of minimal index.
   * @pre 0 <= i < size() && 0 <= j < size() && i != j
   * @post size() == old size() - 1
   */
  void merge(size_t i, size_t j);

  /**
   * Merge breaks i and i+1 in break i
   */
  void merge(size_t i);

  size_t conditional_neighbor_merge(size_t i, std::function <bool(size_t, size_t)> f);
};
