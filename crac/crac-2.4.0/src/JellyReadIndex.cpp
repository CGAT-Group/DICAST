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

#include <config.h>

#ifdef HAVE_LIBJELLYFISH
#include "JellyReadIndex.h"

uint JellyReadIndex::getFactorLength() const {
  return mer_dna::k();
}

uint JellyReadIndex::getNbTags() const {
  // TODO this method is not implemented right now!!!
  return 0;
}

uint *JellyReadIndex::getSupport(Read *r) const {
  size_t support_length = r->getLength() - getFactorLength() + 1;
  uint* support = new uint[support_length]();
  
  mer_dna          m, rm;
  unsigned int length = 0;
  uint64_t val = 0;
  uint j = 0;
  int last_n = -1;

  for(uint i = 0; i < r->getLength(); i++) {
    if(m.shift_left(r->seq[i]) == 'N') {
      last_n = i;
      //continue;
    }
    rm.shift_right(mer_dna::complement(r->seq[i]));
    if(++length < mer_dna::k()) {
      //support[j] = 0;
      //j++;
      continue;
    }

    // Check if current k-mer has a N
    if(last_n != -1 && i - last_n < mer_dna::k()) {
      support[j] = 0;
    } else {
      val = 0;
      if(!ary->ary()->get_val_for_key(m < rm ? m : rm, &val))
        val = 0;
      support[j] = val;
    }
    j++;
  }

  return support;
}
#endif //JELLYREADINDEX_H
