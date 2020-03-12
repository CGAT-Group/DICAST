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

#ifndef CANDIDAT_BREAK_H
#define CANDIDAT_BREAK_H

#include <utility>
#include "Parameters.h"
#include "libSSA/locateOnGenome.h"
#include "libSSA/chrPosition.h"
#include "types.h"

class SupportBreak;

class CandidatBreak {

  uint pos_start_break;		/* Position in the support of the start of
				   the break.*/
  uint loc_start_break;           /* the last absolute location before the 
				   break */

  int strand_start_break;       /* the strand of the start_break position */

  uint pos_end_break;		/* Position in the support of the end of 
				   the break. */
  
  uint loc_end_break;   /* the first absolute location after the 
			   break */
  

  int strand_end_break;         /* the strand  of the end_break position*/

  SupportBreak* supportBreak; 	    /* The supportBreak we are working on. */

  Parameters* parameters;	/* Parameters chosen by the user or chosen 
				   by default by the great developers. */

  bool is_duplicated;             /* true iff there are at least two possible causes for a same break */

  bool is_single_consistent;             /* true if the candidat locations follow the single location */

  long long int single_distance;           /* distance between (loc_start_break || loc_end_break) and loc_single in case of is_single_consistent */

 public:

  CandidatBreak(uint pos_start, uint loc_start, int strand_start, 
		uint pos_end, uint loc_end, int strand_end,
		SupportBreak *sb, Parameters *p);


  CandidatBreak(const CandidatBreak &c); 

  ~CandidatBreak();

  /**
   * Check the correspondence between Candidat and the single Location.
   * - same strand
   * - same chr
   * - loc_single <= loc_start_break if loc_single is before break
   *   or loc_single => loc_end_break if loc_single is after break  
   */
  void checkSingleCorrespondence(uint pos_single, uint loc_single, int strand_single);

  /**
   * @return chr for the pos_end_break
   */
  const uchar *getChrEndBreak();

  /**
   * @return a ChrPosition for the pos_end_break
   */
  ChrPosition *getChrPositionEndBreak();

  /**
   * @return a ChrPosition for the pos_start_break
   */
  ChrPosition *getChrPositionStartBreak();

  /**
   * @return chr id for the pos_end_break
   */
  ulong getChrIdEndBreak();

  /**
   * @return chr id for the pos_start_break
   */
  ulong getChrIdStartBreak();

  /**
   * @return chr for the pos_start_break
   */
  const uchar *getChrStartBreak();
  
  /**
   * @return the reference genome.
   */
  LocateOnGenome *getGenome();

  /**
   * @return the distance between loc_start_break and loc_end_break
   *         on the reference genome.
   */
  gap_size_t getGenomeGapLength();
  
  /**
   * @return the support length: 
   *         supportBreak->getLength()
   */
  uint getLength();

  /**
   * @return loc_end_break
   */
  uint getLocEndBreak();

  /**
   * @return loc_start_break
   */
  uint getLocStartBreak();

  /**
   * @return pos_end_break
   */
  uint getPosEndBreak();

  /**
   * @return pos_start_break
   */
  uint getPosStartBreak();

  /**
   * @return the difference between pos_start_break and pos_end_break
   *         on the read (the break length).
   */
  uint getReadBreakLength();

  /**
   * @return single_distance;
   */
  long long int getSingleDistance();

  /**
   * @return strand_end_break
   */
  int getStrandEndBreak();

  /**
   * @return strand_start_break
   */
  int getStrandStartBreak();

  /**
   * @return true if pos_end_break == (getReadLength() -1)
   */
  bool hasNoEndBreak();

  /**
   * @return true if pos_start_break == 0
   */
  bool hasNoStartBreak();
  
  /**
   * @return true if the candidat is duplicated (two similar candidats 
   *         for a same SupportBreak),
   *         false otherwise.
   */
  bool isDuplicated();

  /**
   * @return true is the candidat is a SNV, an indel or a splice,
   *         false is the candidat looks like a chimera. 
   */
  bool isNiceCandidat();

  /**
   * @return true is loc_start_break and loc_end_break are on the same chr,
   *         false otherwise.
   */
  bool isSameChr();
  
  /**
   * @return true is loc_start_break and loc_end_break are on the same strand,
   *         false otherwise.
   */
  bool isSameStrand();

  /**
   * @return true is the candidat is single consistent 
   *                         (see checkSingleCorrespondence procedure)
   *         false otherwise.
   */
  bool isSingleConsistent();
 

  /**
   * isDuplicated = flag.
   */
  void setDuplicated(bool flag);
  
  /**
   * loc_end_break = loc_end.
   */
  void setLocEndBreak(uint loc_end);
  
  /**
   * loc_start_break = loc_start.
   */
  void setLocStartBreak(uint loc_start);

  /**
   * pos_end_break = pos_end.
   */
  void setPosEndBreak(uint pos_end);
  
  /**
   * pos_start_break = pos_start.
   */
  void setPosStartBreak(uint pos_start);
  
  /**
   * isSingleConsistent = flag.
   */
  void setSingleConsistent(bool flag);

  /**
   * strand_end_break = flag.
   */
  void setStrandEndBreak(int flag);

  /**
   * strand_start_break = flag.
   */
  void setStrandStartBreak(int flag);

};
#endif
