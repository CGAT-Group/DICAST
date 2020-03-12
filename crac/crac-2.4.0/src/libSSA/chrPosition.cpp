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

#include "chrPosition.h"

ChrPosition::ChrPosition():rel_position(0), chr(NULL), strand(0) {}

ChrPosition::ChrPosition(ulong rel_position, char *chr, int chr_id, int strand):
rel_position(rel_position), chr(chr), chr_id(chr_id), strand(strand) {}

ChrPosition::ChrPosition(const ChrPosition &pos) {
  rel_position = pos.getRelativePosition();
  chr = pos.getChrPosition();
  chr_id = pos.getChrId();
  strand = pos.getStrand();
}

ulong ChrPosition::getRelativePosition() const {
  return rel_position;
}

char *ChrPosition::getChrPosition() const {
  return chr;
}

int ChrPosition::getChrId() const {
  return chr_id;
}

int ChrPosition::getStrand() const {
  return strand;
}

void ChrPosition::setStrand(int new_strand) {
  strand = new_strand;
}

ostream &operator<<(ostream &os, ChrPosition &p) {
  os << p.getChrPosition() << "|"<< p.getStrand() << "," << p.getRelativePosition();
  return os;
}

ChrPosition ChrPosition::operator+(int i) {
  return ChrPosition(*this) += i;
//   return pos;
}

ChrPosition &ChrPosition::operator+=(int i) {
  rel_position += i;
  return *this;
}

ChrPosition ChrPosition::operator-(int i) {
  return ChrPosition(*this) -= i;
}

ChrPosition &ChrPosition::operator-=(int i) {
  *this += -i;
  return *this;
}

ChrPosition ChrPosition::operator+(ulong i) {
  return ChrPosition(*this) += i;
}

ChrPosition &ChrPosition::operator+=(ulong i) {
  rel_position += i;
  return *this;
}

ChrPosition ChrPosition::operator-(ulong i) {
  return ChrPosition(*this) -= i;
}

ChrPosition &ChrPosition::operator-=(ulong i) {
  *this += -i;
  return *this;
}
