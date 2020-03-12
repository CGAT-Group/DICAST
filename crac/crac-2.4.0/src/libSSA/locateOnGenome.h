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

#ifndef LOCATE_ON_GENOME
#define LOCATE_ON_GENOME

extern "C" {
#include <sys/time.h>
}

#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>
#include "SSA.h"
#include "chrPosition.h"
#include "GenomeInfo.h"
using namespace std;


class LocateOnGenome {
 private:
  GenomeInfo *genomeInfo;
  ostream *out_uniq, *out_multiple, *out_none, *out_sam;
  TFMindex *index;
  ulong nb_tags_treated;
  ulong nb_displayed;
  ulong nbOccs_lastTag;
  time_t chrono_sec;
  suseconds_t chrono_usec;
  bool display_nb;

  void display_occ(ostream *stream, ulong *occ, ulong n, ulong tag_length, bool sens);


 public:
  /**
   * @param file_base: Basename of the index file (ie. the common prefix
   *                   to the .ssa index and to the .conf file).
   */ 
  LocateOnGenome(uchar *file_base);

  ~LocateOnGenome();

  void displaySingleLocation(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens);
  void displayMultipleLocations(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens);
  void displayNoLocation(uchar *tag, ulong tag_length);

  float getElapsedTime();


  void getLocations(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens, ulong **poccs_sens, ulong **poccs_antisens);
  
  /**
   * @return all locations of a kmer at position_of_locations and its
   * number of occurrences.  In addition, position_of_locations is set
   * after pair<ChrPosition **,uint> which determines position of the
   * k-mer chosen for locations
   */
  pair<ChrPosition **,uint>getLocations(uchar *tag, uint klength, uint read_length, uint &position_of_locations);
  
 
  uint getPositionOfLocations();

  uint getNbLocations();

  /* 
   * @param chrPos : a ChrPosition object
   * @return the absolute position of chrPos
   */
  ulong getAbsolutePosition(ChrPosition *chrPos);

  /**
   * @param pos: absolute position
   * @return the chromosome which corresponds to that position
   */
  const uchar *getChromosome(ulong position);

  /**
   * @param chr_name: name of chr
   * @return an id of the chromosome which corresponds to that name.
   */
  ulong getIdChromosomeFromName(char *chr_name);

  /**
   * @param pos: absolute position
   * @return an id of the chromosome which corresponds to that position.
   * This id is mainly useful when one wants to compare two positions
   * and wants to know if they're on the same chromsome.
   */
  ulong getIdChromosome(ulong position);

  /**
   * @param pos: absolute position
   * @return the relative position on the genome.
   */
  ChrPosition *getChrPos(ulong pos, int strand=0);

  pair<uint, uint> getFMIndexRange(uchar *tag, ulong tag_length);

  pair<uint, uint> getFMIndexReverseRange(uchar *tag, ulong tag_length);

  void getOccurrencesFMIndexRange(ulong sp, ulong ep, ulong **occ, ulong m=0, ulong nb_disp=0);
  
  /**
   * @return A pointer to the GenomeInfo object for hvaing all the information
   *         on sequence length, name, IDs, ...
   */
  GenomeInfo * getGenomeInfo();

  uchar *getGenomeSubstring(ulong pos, ulong length, int strand=1,
                            char *chr_name = NULL);


  ulong getNbTagsTreated();
  void getNbOccurrences(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens);
  ulong getNbOccsLastLocatedTag();

  void locateTags(fstream &file, ulong tag_length, ulong nb_threads = 1);

  void locateTag(uchar *tag, ulong tag_length);


  /**
   * @param disp: true iff read number precedes each entry in the output files
   */
  void setDisplayNumbers(bool disp);

  void setOutputStreams (ostream *uniq, ostream *multiple, ostream *none);
  void setOutputUniq(ostream *s);
  void setOutputMultiple(ostream *s);
  void setOutputNone(ostream *s);
  void setOutputSAM(ostream *s);

  void setNbLocations(ulong nb_displayed);

  bool startChrono();
  
};

#endif
