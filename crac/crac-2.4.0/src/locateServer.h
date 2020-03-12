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

#ifndef LOCATE_SERVER_H
#define LOCATE_SERVER_H
#include <iostream>
#include "types.h"
using namespace std;

/**
 * @return the number of digits of a in base 10
 */
uint nb_digits(uint a);

/**
 * Write in os the number writeIt with whitespaces preceding it.
 * The number of whitespaces is chosen such that writeIt, a and b
 * are written on the same number of columns (including the whitespaces).
 * @param os: the output stream
 * @param writeIt: the number to be output
 * @param a: another number that will be written (or has already been written)
 *           on the same column.
 * @param b: another number that will be written (or has already been written)
 *           on the same column.
 */
ostream &writeWithSpaces(ostream &os, uint writeIt, uint a, uint b);


/**
 * Create a subserver that will be dedicated to a specific client.
 * @params args: an array that contains NB_ARGS_SUBSERVER parameters.
 */ 
void *launchSubServer(void *args);

/**
 * Create the main server.
 * @param ct: ClassifyTags object for classifying reads
 * @params fifo_name: fifo input name for the main server
 * @params fifo_output_name: fifo output name for the main server.
 */
void locateServer(ClassifyTags *ct,
                 const char *fifo_name, const char *fifo_output_name);


#endif
