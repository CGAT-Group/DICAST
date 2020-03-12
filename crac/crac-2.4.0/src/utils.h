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

#ifndef UTILS_H
#define UTILS_H
#include "types.h"
#include <gkArraysTypes.h>
#include <sys/time.h>
#include <string>
#include <fstream>

#define PRINT_VAR(s) cerr << #s << " = " << s << endl;

uint convNuc(char nuc);
char intToNuc(uint c);

char majNuc(char nuc);

char *majNucNCpy(char *dest, const char *source, size_t n);

uintSA DNAtoInt(char *dna, uintSA dna_length);

void intToDNA(uint64 code, uint dna_length, char *dna);


int comparUint(const void *a1, const void* a2);

bool pairCompare(const std::pair<uint, uint>& firstElem, const std::pair<uint, uint>& secondElem);

uint64 factorsToInt(uintSA pos, uchar *dna, uint length);

uint posBitInMask(uint mask);

float getChrono();

std::string intToString(uint i);

/**
 * Create an output stream, that is gzipped if required.
 */
std::ostream *create_stream(bool gzipped, char *name);

/* Used to name the semaphores for each thread */
std::string concatStringAndInt(const char*str, int i);

/* Used to name the semaphores for each thread */
std::string concatStringAndFloat(const char*str, float i);

#endif
