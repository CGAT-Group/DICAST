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

#include "SNPInfo.h"
#include "libSSA/utils.h"

SNPInfo::SNPInfo(ChrPosition chrPos, uint position,
		 float score, bool is_duplicated, 
		 char actual_nuc, char expected_nuc):
  chrPos(chrPos), score(score),
  is_duplicated(is_duplicated),
  actual_nuc(actual_nuc),
  expected_nuc(expected_nuc) 
{
  setPosition(position);
  if (actual_nuc == 0) {
    if (expected_nuc == 0) {
      type = SNP_UNKNOWN;
    } else {
      type = SNP_DELETION;
    }
  } else {
    if (expected_nuc == 0) {
      type = SNP_INSERTION;
    } else {
      type = SNP_SUBSTITUTION;
    }
  }
}

SNPInfo::~SNPInfo() {
}

char SNPInfo::getActualNucleotide(int strand) {
  if (getSNPType() == SNP_UNKNOWN || getSNPType() == SNP_DELETION)
    return '?';
  return (strand == 1) ? actual_nuc : complementDNA(actual_nuc);;
}

char SNPInfo::getExpectedNucleotide(int strand) {
  if (getSNPType() == SNP_UNKNOWN || getSNPType() == SNP_INSERTION)
    return '?';
  return (strand == 1) ? expected_nuc : complementDNA(expected_nuc);
}

ChrPosition &SNPInfo::getChrPosition() {
  return chrPos;
}

float SNPInfo::getScore() {
  return score;
}

snp_type SNPInfo::getSNPType() {
  return type;
}

bool SNPInfo::isDuplicated() {
  return is_duplicated;
}

void SNPInfo::output(ostream &os) {
  os << *this;
}

void SNPInfo::samOutput(ostream &os, int strand) {
  os << "SNP:" << getScore() << ":" << getPosition(strand) << ":" << getChrPosition()
     << ":" << getExpectedNucleotide(strand) << ":" << getActualNucleotide(strand);
}

cigar_type SNPInfo::getCigarInfo() {
  cigar_type cigar;
  cigar.nb = 1;
  cigar.type = 0;               // Just initializing so that gcc is happy
  if (getSNPType() == SNP_SUBSTITUTION) {
    cigar.type = 'X';
  } else if (getSNPType() == SNP_DELETION) {
    cigar.type = 'D';
  } else if (getSNPType() == SNP_INSERTION) {
    cigar.type = 'I';
  }
  return cigar;
}

ostream &operator<<(ostream &os, SNPInfo& i) {
  if (i.isDuplicated()){
    os << "duplicate ";  
  }else{
    os << "single ";
  }
  os << i.getChrPosition() 
     << " pos_SNV="<<i.getPosition() 
     << " ";
  if (i.getSNPType() == SNP_UNKNOWN)
    os << "?->?";
  else {
    if (i.getSNPType() == SNP_INSERTION)
      os << "?";
    else 
      os << i.getExpectedNucleotide() ;
    
    os << "->";
    
    if (i.getSNPType() == SNP_DELETION)
      os << "?";
    else 
      os << i.getActualNucleotide();
    os << " score=" << i.getScore() ; 
	  
  }
  return os;
}
