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

#ifndef GENOME_INFO_H
#define GENOME_INFO_H

#include <utility>
#include <iostream>
#include "basic.h"

#include "chrPosition.h"

using namespace std;

/**
 * Maximal length of a chromosome (or sequence) name.
 */
#define MAX_LENGTH_CHR_NAME 60

/**
 * A class that gives access to the information stored in .conf file
 */
class GenomeInfo {
private:
  /**
   * Array containing the sequence names
   */
  uchar **chr_names;
  /**
   * Number of sequences
   */
  ulong nb_chr;
  /**
   * Array containing the length (in nt) of each sequence
   */
  ulong *length_chr;
  /**
   * Array containing the cumulative length of each sequence.
   * First cell's value is 0
   */
  ulong *total_length;


public:

  /**
   * Read the file whose name is stored in <basename>.conf.
   * Retrieves the information from the file and store them 
   * in the class.
   * @param basename: name of the file without the .conf extension
   */
  GenomeInfo(char *basename);

  /**
   * Construct the class given the sequences names, their number
   * and their lengths.
   * Beware: the arrays are not copied.
   * @param chr_names: Names of the sequences
   * @param nb_chr: Number of sequences
   * @param length_chr: Length of the sequences
   */
  GenomeInfo(uchar **chr_names, ulong nb_chr, ulong *length_chr);

  ~GenomeInfo();

  /**
   * @pre i >= 0 && i < getNbChr()
   * @return The name of the sequence number i
   */
  uchar *getChrName(ulong i);

  /**
   * @pre i >= 0 && i < getNbChr()
   * @return The length of the sequence number i (whose name is getChrName(i))
   */
  ulong getChrLength(ulong i);

  /**
   * @return the total length of the chromosome ie. the cumulated length of 
   *         the chromosomes.
   */
  ulong getGenomeLength();

  /**
   * @param chr: Name of the chromosome
   * @param position: Position in the chromosome
   * @return The position in the genome of <position> in chromosome <chr>
   *         or ~0 if no such chromosome or position in the chromosome exist.
   * @post p = getChrNameFromPosition(pos) 
   *       <==> pos == getGenomePosition(p.first, p.second)
   */
  ulong getGenomePosition(uchar *chr, ulong position);

  /**
   * @param id: Id of the chromosome
   * @param position: Position in the chromosome
   * @return The position in the genome of <position> in chromosome number <id>
   *         or ~0 if no such chromosome or position in the chromosome exist.
   * @post p = getChrNameFromPosition(pos) 
   *       <==> pos == getGenomePosition(p.first, p.second)
   */
  ulong getGenomePosition(ulong id, ulong position);

  /**
   * @return The number of sequences stored
   */ 
  ulong getNbChr();

  /**
   * @param name: Name of a sequence
   * @return the number associated to the sequence or getNbChr() if no such
   *         name corresponds to a known sequence.
   * @post  getChrName(getNumFromName(name)) and name are equals
   */ 
  ulong getNumFromName(uchar *name);


  /**
   * DEPRECATED
   *
   * please use getChrPositionFromPosition instead.
   *
   * @param pos: a position in the genome
   * @return the name of the chromosome associated to that position and
   *         the relative position in the chromosome.
   *         If the position isn't in the genome, the chromosome name is set
   *         to NULL
   */
  pair<const uchar *,const ulong> getChrNameFromPosition(ulong pos);

  ChrPosition *getChrPositionFromPosition(ulong pos);

  /**
   * @param pos: a position in the genome
   * @return the number of the chromosome associated to that position, and
   *         the relative position in the chromosome.
   *         If the position isn't in the genome, the chromosome number is
   *         set to getNbChr().
   * @post getNumFromName(getChrNameFromPosition(pos).first)
   *       == getNumFromPosition(pos).first
   */
  pair<const ulong, const ulong> getNumFromPosition(ulong pos);

  /**
   * @return true iff pos < getGenomeLength()
   */
  bool isValidPosition(ulong pos);

  /**
   * Save to the file whose name given in parameter
   * @param filename: the filename you must save to
   * @param ext: extension put after the filename (default .conf)
   */
  void save(char *name, const char *ext=".conf");
};

ostream &operator<<(ostream &os, GenomeInfo &g);

#endif
