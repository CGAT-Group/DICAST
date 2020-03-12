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

#ifndef SNP_INFO_H
#define SNP_INFO_H
#include "../types.h"
#include "genericInfo.h"
#include "../libSSA/chrPosition.h"
#include "cigarTypes.h"

class SNPInfo : public GenericInfo {
private:
  ChrPosition chrPos;
  float score;
  snp_type type;
  bool is_duplicated;
  char actual_nuc;
  char expected_nuc;

public:
  SNPInfo(ChrPosition chrPos, uint position, float score,
	  bool is_duplicated,
	  char actual_nuc=0, char expected_nuc=0);
  ~SNPInfo();

  /**
   * Return the nucleotide actually found in the sequence,
   * or 0 if getSNPType() == SNP_DELETION || getSNPType() == SNP_UNKNOWN
   * @param strand: if strand == -1 return the complementary of the nucleotide
   *                instead of the nucleotide itself
   */
  char getActualNucleotide(int strand=1);

  /**
   * Return the nucleotide as expected in the reference genome
   * or 0 if getSNPType() == SNP_INSERTION || getSNPType() == SNP_UNKNOWN
   * @param strand: if strand == -1 return the complementary of the nucleotide
   *                instead of the nucleotide itself
   */
  char getExpectedNucleotide(int strand=1);

  ChrPosition &getChrPosition();
  float getScore();
  snp_type getSNPType();
  bool isDuplicated();
  
  void output(ostream &);
  void samOutput(ostream &, int=1);
  cigar_type getCigarInfo();
};

ostream &operator<<(ostream &os, SNPInfo &);

#endif
