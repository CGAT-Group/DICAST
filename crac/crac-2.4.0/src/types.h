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

#ifndef TYPES_H
#define TYPES_H
#include <sys/types.h>

typedef enum enum_snp_type {SNP_UNKNOWN, 
			    SNP_SUBSTITUTION,
			    SNP_INSERTION,
			    SNP_DELETION} 
  snp_type;

typedef enum enum_error_context {FIRST_SUBSTITUTION, /* First sub in a break 
                                                      * meaning there may
                                                      * have more
                                                      * to come */
                                 SECOND_SUBSTITUTION,
                                 INSERTION,
                                 DELETION,
                                 UNKNOWN_START_MISSING,
                                 UNKNOWN_END_MISSING}
  error_context;

typedef enum enum_pos_location_break {START_BREAK, END_BREAK}
pos_location_break;
#define NB_POS_BREAK 2

typedef enum eenum_strand {FORWARD_STRAND, REVERSE_STRAND}
strand_t;
#define NB_ENUM_STRAND  2

typedef enum enum_paired_end_orientation {FORWARD_REVERSE, REVERSE_FORWARD, FORWARD_FORWARD}
paired_end_orientation_t;
#define NB_ENUM_PAIRED_END_ORIENTATION 3

/* #ifndef BUILD_CRAC64 */
/* #define uintSA u_int32_t */
/* #else  */
/* #define uintSA u_int64_t */
/* #endif */

#ifndef __USE_MISC
typedef u_int32_t uint; 
#endif

#define ulong uint

typedef u_int64_t uint64;
typedef int64_t longint;
typedef u_int8_t byte;

typedef unsigned char uchar;

#ifndef gap_size_t
#define gap_size_t longint
#define GAP_SIZE_MAX  (~((gap_size_t)1 << 63))

#endif

#endif
