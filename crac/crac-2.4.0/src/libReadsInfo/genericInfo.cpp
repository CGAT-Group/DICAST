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

#include "genericInfo.h"

GenericInfo::GenericInfo() {
  stranded_pos = ~0;
}

GenericInfo::~GenericInfo() {}

void GenericInfo::cigarString(ostream &os) {
  os << getCigarInfo();
}

cigar_type GenericInfo::getCigarInfo() {
  cigar_type cigar;
  cigar.nb = 0;
  cigar.type = 0;
  return cigar;
}

uint GenericInfo::getPosition(int strand) {
  return (strand == 1) ? position 
    : (((uint)stranded_pos != (uint)~0) 
       ? stranded_pos 
       : getReadLength() - 1 - position);
}

uint GenericInfo::getReadLength() {
  return read_length;
}

void GenericInfo::setPosition(uint p) {
  position = p;
}

void GenericInfo::setReadLength(uint l) {
  read_length = l;
}

void GenericInfo::setStrandedPosition(uint pos) {
  stranded_pos = pos;
}

uint GenericInfo::getBreakId() {
  return breakId;
}

void GenericInfo::setBreakId(uint id) {
  breakId = id;
}

uint GenericInfo::getCauseId() {
  return causeId;
}

void GenericInfo::setCauseId(uint id) {
  causeId = id;
}
