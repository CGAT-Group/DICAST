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

#ifndef BIO_UNDETERMINED_H
#define BIO_UNDETERMINED_H
#include "genericInfo.h"
#include "../libSSA/chrPosition.h"

class BioUndeterminedInfo : public GenericInfo {
  public:
    enum Type {
    UNEXPECTED_CASE               = 1,
    BREAK_TOO_SMALL               = 2,
    BREAK_TOO_LARGE               = 3,
    NO_START_BREAK                = 4,
    NO_END_BREAK                  = 5,
    COMBINATION_ERROR_BIOLOGICAL  = 6,
    REPETITION_AROUND_THE_BREAK   = 7,
    AMBIGUOUS_SNP                 = 8,
    AMBIGUOUS_INDEL               = 9,
    AMBIGUOUS_SPLICE              = 10,
    AMBIGUOUS_CHIMERA             = 11,
    WEAK_CHIMERA                  = 12,
    DISCORDANT_PAIR               = 13,
    LOW_SUPPORT_INTRA_EVENT       = 14,
    };

 private:
  ChrPosition *pos;
  uint position;
  Type type;
  char *message;
  bool detailed;
 public:

  BioUndeterminedInfo(ChrPosition *pos, uint position, Type, char* message, bool detailed = false);
  ~BioUndeterminedInfo();

  ChrPosition *getChrPosition();
  char *getMessage();

  void output(ostream &);
  void samOutput(ostream &, int=1);
};

ostream &operator<<(ostream &os, BioUndeterminedInfo &);

#endif
