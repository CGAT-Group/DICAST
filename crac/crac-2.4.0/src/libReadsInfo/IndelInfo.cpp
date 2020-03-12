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

#include <cassert>
#include "IndelInfo.h"

IndelInfo::IndelInfo(ChrPosition chrPos, uint position, uint nb_ins, uint nb_del,
                     float score, bool is_duplicated):
  chrPos(chrPos),nb_ins(nb_ins),nb_del(nb_del),score(score),is_duplicated(is_duplicated) {
  setPosition(position);
}

IndelInfo::~IndelInfo() {
//   delete chrPos;
}

uint IndelInfo::getNbIns() {
  return nb_ins;
}

uint IndelInfo::getNbDel() {
  return nb_del;
}

ChrPosition &IndelInfo::getChrPosition() {
  return chrPos;
}

float IndelInfo::getScore() {
  return score;
}

bool IndelInfo::isDuplicated() {
  return is_duplicated;
}

void IndelInfo::output(ostream &os) {
  os << *this;
}

void IndelInfo::samOutput(ostream &os, int strand) {
  if (getNbIns() > 0 && getNbIns() < (uint)~0) {
    os << "Ins:" << getScore() << ":" << getPosition(strand) << ":" << getChrPosition()
       << ":" << getNbIns();       
  } else if (getNbDel() > 0 && getNbDel() < (uint) ~0) {
    os << "Del:" << getScore() << ":" << getPosition(strand) << ":" << getChrPosition()
       << ":" << getNbDel();       
  } else {
    cerr << "Nothing was expected here : " << __FILE__ << ":" << __LINE__ << endl;
  }
}

cigar_type IndelInfo::getCigarInfo() {
  cigar_type cigar;
  if (getNbIns() == 0 && getNbDel() == 0) {
    cigar.nb = 1;
    cigar.type = 'X';
  } else if (getNbIns() != (uint)~0 && getNbIns() > 0) {
    cigar.nb = getNbIns();
    cigar.type = 'I';
  } else {
    assert(getNbDel() != (uint)~0);
    cigar.nb = getNbDel();
    cigar.type = 'D';
  }
  return cigar;
};


ostream &operator<<(ostream &os, IndelInfo& i) {
  if (i.isDuplicated()){
    os <<"duplicate ";
  }else{
    os <<"single ";
  }
  os << i.getChrPosition() << " pos_indel="<<i.getPosition()<<" nb_ins=";
  uint nb = i.getNbIns();
  if (nb == (uint)~0)
    os << "unknown";
  else 
    os << nb;
  os <<" nb_del=";
  nb = i.getNbDel();
  if (nb == (uint)~0)
    os << "unknown";
  else 
    os << nb;
  return os;
}
