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

#include "SeqErrInfo.h"
#include "../libSSA/utils.h"

SeqErrInfo::SeqErrInfo(ChrPosition chrPos, uint position, uint nb_ins, 
		       uint nb_del, float score, bool is_duplicated,
		       char *genome, char *tag,
		       uint g_length, uint t_length): 
  IndelInfo(chrPos, position, nb_ins, nb_del, score, is_duplicated),
  genome_sequence(genome),
  err_sequence(tag), g_seq_length(g_length), t_seq_length(t_length)
{
  genome_sequence_revcomp = NULL;
  err_sequence_revcomp = NULL;
  if (genome_sequence) {
    genome_sequence_revcomp = new char[g_seq_length+1];
    genome_sequence_revcomp[g_seq_length] = 0;
    for (uint i=0; i < g_seq_length; i++) 
      genome_sequence_revcomp[g_seq_length - i - 1] = complementDNA(genome_sequence[i]);
  }
  
  if (err_sequence) {
    err_sequence_revcomp = new char[t_seq_length+1];
    err_sequence_revcomp[t_seq_length] = 0;
    for (uint i = 0; i < t_seq_length; i++)
      err_sequence_revcomp[t_seq_length - i - 1] = complementDNA(err_sequence[i]);
  }
}

SeqErrInfo::~SeqErrInfo() {
  if (genome_sequence != NULL && genome_sequence[0]!=0) {
    free(genome_sequence);
    delete [] genome_sequence_revcomp;
  }
  if (err_sequence != NULL && err_sequence[0]!=0) {
    delete [] err_sequence;
    delete [] err_sequence_revcomp;
  }
}

char *SeqErrInfo::getErrorSequence(int strand) {
  return (strand == 1) ? err_sequence : err_sequence_revcomp;
}

uint SeqErrInfo::getErrorSequenceLength() {
  return t_seq_length;
}

char *SeqErrInfo::getGenomeSequence(int strand) {
  return (strand == 1) ? genome_sequence : genome_sequence_revcomp;
}

uint SeqErrInfo::getGenomeSequenceLength() {
  return g_seq_length;
}

// bool SeqErrInfo::isDuplicated() {
//   return is_duplicated;
// }

void SeqErrInfo::samOutput(ostream &os, int strand) {
  if (getErrorSequenceLength() == getGenomeSequenceLength() 
      && getErrorSequenceLength() < (uint)~0)
    // Substitution
    os << "Error:Sub:" << getPosition(strand) << ":" << getScore() << ":" 
       << getGenomeSequence(strand) << ":" 
       << getErrorSequence(strand);
  else if ((getNbIns() == (uint)~0) && (getNbDel() == (uint)~0))
    // Unknown
    os << "Error:Unknown:" << getPosition(strand) << ":" << getScore();
  else if (getNbDel() > 0)
    // Deletion
    os << "Error:Del:" << getPosition(strand) << ":" << getScore() << ":"
       << getNbDel() << ":" 
       << (getGenomeSequence() == NULL ? "<snip>" : getGenomeSequence(strand));
  else if (getNbIns() > 0)
    os << "Error:Ins:" << getPosition(strand) << ":" << getScore() << ":"
       << getNbIns() << ":" 
       << (getErrorSequence() == NULL ? "<snip>" : getErrorSequence(strand));
}

ostream &operator<<(ostream &os, SeqErrInfo& i) {
   if (i.isDuplicated()){
    os <<"duplicate ";
  }else{
    os <<"single ";
  }
  os << i.getChrPosition() 
     << " pos_error="<<i.getPosition()
     << " nb_ins=";
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

  if ((i.getGenomeSequence() == NULL 
       || i.getGenomeSequence()[0]==0)
      && (i.getErrorSequence() == NULL
	  || i.getErrorSequence()[0]==0))
    return os;

  os << " ";

  if (i.getGenomeSequence() == NULL) 
    os << "?";
  else
    os << i.getGenomeSequence();
  
  os << "->";

  if (i.getErrorSequence() == NULL)
    os << "?";
  else
    os << i.getErrorSequence();
  os << " score=" << i.getScore();
  return os;
}
