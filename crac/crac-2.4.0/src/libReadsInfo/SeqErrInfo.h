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

#ifndef SEQERR_INFO_H
#define SEQERR_INFO_H
#include "IndelInfo.h"
#include "../libSSA/chrPosition.h"

class SeqErrInfo : public IndelInfo {
 private:
  float score; 
  char *genome_sequence;
  char *err_sequence;
  char *genome_sequence_revcomp;
  char *err_sequence_revcomp;
  uint g_seq_length;
  uint t_seq_length;

 public:
  SeqErrInfo(ChrPosition chrPos, uint position, uint nb_ins, uint nb_del,
	     float score, bool is_duplicated,
	     char *genome=NULL, char *tag=NULL, uint g_length=~0, uint t_length=~0);

  ~SeqErrInfo();

  /**
   * @return err_sequence, the sequence, corresponding to the error, 
   * on the tag.
   */
  char *getErrorSequence(int strand = 1);

  /**
   * @return t_seq_length, the length of the string returned by
   * getErrorSequence().
   */
  uint getErrorSequenceLength();

  /**
   * @return err_sequence, the sequence, corresponding to the error, 
   * on the genome.
   */
  char *getGenomeSequence(int strand = 1);
  
  /**
   * @return g_seq_length, the length of the array returned by 
   * getGenomeSequence (\0 excluded).
   */
  uint getGenomeSequenceLength();

  
  /* bool isDuplicated(); */

  void samOutput(ostream &os, int=1);
};

ostream &operator<<(ostream &os, SeqErrInfo& i);

#endif
