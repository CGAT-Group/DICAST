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

#include "Read.h"

#ifndef READINDEX_H
#define READINDEX_H

class ReadIndex {

  public:

    /*
     * Virtual desconstructor => mandatory
     */
    virtual ~ReadIndex() {};

    /**************************************
     *          VIRTUAL METHODS
     *************************************/

    /**
     * @return the factor length used to index the reads collection
     */
    virtual uint getFactorLength() const = 0;

    /**
     * @return the number of tags (or reads) indexed in the reads index structure
     */
    virtual uint getNbTags() const = 0;

    /**
     * @param a Read object instance
     * @return an array whose length is r->getLength() - getFactorLength()
     *         and where the value at position i is the number of occurrences
     *         of the given k-factor that starting at position i 
     *         among all indexed k-factors in the reads collection
     */
    virtual uint *getSupport(Read *r) const = 0;

};

#endif //READINDEX_H

