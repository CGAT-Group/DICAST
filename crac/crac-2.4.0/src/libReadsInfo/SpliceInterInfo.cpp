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

#include "SpliceInterInfo.h"

SpliceInterInfo::SpliceInterInfo(ChrPosition chrPos, uint position, ChrPosition chrPos2, 
				 bool is_duplicated, int chim_class, float chimera_score, const string chimera_score_info):
  chrPos(chrPos), chrPos2(chrPos2), is_duplicated(is_duplicated), chim_class(chim_class),
  chimera_score(chimera_score), chimera_score_info(chimera_score_info) {
  setPosition(position);
  is_valid = true;
}

SpliceInterInfo::~SpliceInterInfo() {
}

ChrPosition &SpliceInterInfo::getChromosomeDest() {
  return chrPos2;
}

ChrPosition &SpliceInterInfo::getChrPosition() {
  return chrPos;
}

int SpliceInterInfo::getChimeraClass() {
  return chim_class;
}

float SpliceInterInfo::getChimeraScore(){
  return chimera_score;
}

string SpliceInterInfo::getChimeraScoreInfo() const{
  return chimera_score_info;
}

cigar_type SpliceInterInfo::getCigarInfo() {
  cigar_type cigar;

  if(chim_class == 2) {
    if(getChrPosition().getStrand() == 1) {
    cigar.nb = getChromosomeDest().getRelativePosition() - getChrPosition().getRelativePosition() - 1;
    } else {
    cigar.nb = getChrPosition().getRelativePosition() - getChromosomeDest().getRelativePosition() - 1;
    }
  } else {
    cigar.nb = 0;
  }
  cigar.type = 'N';
  return cigar;
}

bool SpliceInterInfo::isDuplicated() {
  return is_duplicated;
}

void SpliceInterInfo::output(ostream &os) {
  os << *this;
}

void SpliceInterInfo::samOutput(ostream &os, int strand) {
  // if (strand == 1){
    os << "chimera:" << getPosition(strand) << ":" << getChrPosition() << ":"
       << getChromosomeDest() << ":" << chim_class << ":" << chimera_score << ":" << chimera_score_info;
  // }else{
  //   os << "chimera:" << getPosition(strand) << ":" <<  getChromosomeDest() << ":"
  //      << getChrPosition() << ":" << chim_class << ":" << chimera_score << ":" << chimera_score_info;
  // }
}

ostream &operator<<(ostream &os, SpliceInterInfo &i) {
  if (i.isDuplicated()){
    os << "duplicate ";  
  }else{
    os << "single ";
  }
  os  << i.getChimeraClass() << " " << i.getChimeraScore() << " " << i.getChimeraScoreInfo() << " " << i.getChrPosition() 
      << " " <<i.getChromosomeDest() << " pos_junction="<<i.getPosition();
  return os;
}
