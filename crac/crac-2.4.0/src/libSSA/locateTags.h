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

#ifndef LOCATETAGS_H
#define LOCATETAGS_H

#include "locateOnGenome.h"
#include "../Parameters.h"
#include "../utils.h"
#include <ostream>
#include <string>
#include <readsReader.h>
#include <stdlib.h>
#include <fcntl.h> /* For O_* constants */
#include <semaphore.h>
#include "../libReadsInfo/samLine.h"
#include "htslib/sam.h"

using namespace gkarrays;

typedef struct read_s {
  string name; /* Name of the read in the input file */
  string seq;  /* DNA sequence of the read */
  string qual; /* Quality of the DNA sequence */
} read_t;

class LocateTags{
 private:
  LocateOnGenome *genome; /* Genome object to locate reads on reference genome */
  Parameters *parameters; /* Parameters object */
  char *tags_file;
  uint nb_single;         /* Number of single locations */
  uint nb_multiple;       /* Number of multiple locations */
  uint nb_none;           /* Number of None locations */
  uint nb_tags;           /* Number of tags in the file */
  uint threshold;         /* Threshold to use */
  //sem_t stats_mutex;      /* Mutex when thread updates stats values (nb_single, nb_multiple...) */
  sem_t *stats_mutex;      /* Mutex when thread updates stats values (nb_single, nb_multiple...) */
  
 public: 
 
  /**
   * Constructor of locateTags
   *
   * @param genome the genome object used for getting locations
   * @param param the parameters object
   * @param tags_file the name of the FAST/FASTQ file containing reads to process
   * @param sam output stream to print the sam file
   * @param threshold the size of the prefix of the read we are going to use for mapping
   *                  if the size equal 0 locateTags will map the whole sequence of the read.
   */
  LocateTags(LocateOnGenome *genome, Parameters *param, char *tags_file, uint threshold = 0);
  
  /**
   * Destructor of locateTags
   */
  ~LocateTags(){};

  /**
   * Launch location process, and multiple threads if any.
   */
  void locate(samFile *sam,bam_hdr_t* sam_h);
  
  /**
   * locate the tag in the genome index and generate the coresponding sam line.
   *
   * @param read the read to process
   * @return a vector containing the sam lines of the mapping of this read
   */
  vector<SamLine*> *locateTag(read_t *read);
  
  /**
   * @return the number of tags that have been processed
   */
  uint getNbTags();

  /**
   * @return the number of single locations
   */
  uint getNbSingle();
  
  /**
   * @return the number of multiple locations
   */
  uint getNbMultiple();
  
  /**
   * @return the number of none locations
   */
  uint getNbNone();

 private:
  /**
   * Set a read_t object with info from iterato
   *
   * @readIt the iterator used to fill data of the read object
   * @read the read to fill
   */
  void setReadInfo(readIterator *readIt, read_t *read);

  /**
   * The function that is call for each thread
   */
  friend void *fillLocations(void *);

};  

void *fillLocations(void *);

typedef struct locate_thread_s {
  LocateTags *locate;
  sem_t *remaining_seats;  /* Number of remaining cells containing reads to process */
  sem_t *processed_reads;  /* Number of reads that have been processed and ready to write in the ouput */
  vector<SamLine*> **seats; /* Seats where the samLines of the processed reads are stored */
  read_t **reads;          /* Seats where the reads to process are stored */
  uint total_seats;        /* Total number of seats */
  bool keepGoing;          /* Boolean to stop the tread */
} locate_thread_t;

#endif
