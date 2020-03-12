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
#include "bamUtils.h"

int bam_reg2bin(int beg, int end)
{
  --end;
  if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
  if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
  if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
  if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
  if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
  return 0;
}

int cigarop2int(char op) {
  int cigar_op_mask = 0;
  switch ( op ) {
    case 'M':
      cigar_op_mask = BAM_CMATCH;
      break;
    case 'I':
      cigar_op_mask = BAM_CINS;
      break;
    case 'D':
      cigar_op_mask = BAM_CDEL;
      break;
    case 'N':
      cigar_op_mask = BAM_CREF_SKIP;
      break;
    case 'S':
      cigar_op_mask = BAM_CSOFT_CLIP;
      break;
    case 'H':
      cigar_op_mask = BAM_CHARD_CLIP;
      break;
    case 'P':
      cigar_op_mask = BAM_CPAD;
      break;
    case '=':
      cigar_op_mask = BAM_CEQUAL;
      break;
    case 'X':
      cigar_op_mask = BAM_CDIFF;
      break;
    case 'B':
      cigar_op_mask = BAM_CBACK;
      break;
    default:
      std::cerr << "Wrong cigar operator" << std::endl;
  }
  return cigar_op_mask;
}
