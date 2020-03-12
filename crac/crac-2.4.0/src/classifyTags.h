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

#include <cstdlib>
#include "ReadIndex.h"
#include "Read.h"
#include "const.h"
#include "Parameters.h"
#include "libSSA/locateOnGenome.h"
#include "Classifier.h"
#include <fcntl.h> /* For O_* constants */
#include <semaphore.h>
#include "libReadsReader/readsReader.h"
#include "htslib/sam.h"

#ifndef CLASSIFY_TAGS
#define CLASSIFY_TAGS
using namespace std;

class ClassifyTags {

 private:
  LocateOnGenome *genome;
  ReadIndex *tags;
  readsReader *reads_reader;
  Parameters *parameters;
  uint nb_reads, nb_single, nb_multiple, nb_none, nb_duplication;

  /* Statistics */
  uint nb_classes[NB_MASKS+1];
  uint nb_explainable;

 public:
  ClassifyTags(LocateOnGenome *genome, 
	       ReadIndex *tags,
	       readsReader *reads_reader,
	       Parameters *param);

  /* Member access */
  LocateOnGenome *getGenomeIndex();
  Parameters *getParameters();
  ReadIndex *getReadIndex();
  readsReader *getReadsReader();

  /* Statistics */
  uint getNbReads();
  uint getNbNone();
  uint getNbMultiple();
  uint getNbExplainable();
  uint getNbAlmostNormal();
  uint getNbNormal();
  uint getNbDuplication();
  uint getNbSingle();
  uint getNbSNP();
  uint getNbBioTagIndel();
  uint getNbSeqErr();
  uint getNbRepetition();
  uint getNbSplice();
  uint getNbSpliceIntra();
  uint getNbSpliceInter();
  uint getNbSpliceNoCover();
  uint getNbNothing();
  uint getNbUndetermined();
  uint getNbBioUndetermined();
  uint getNbPairedEndChimera();
  

  /**  
   * variables used for duplication:
   * - max_localisation_duplication
   * - min_localisation_duplication
   * - percent_min_unique_repetition
   */

  void classify(ostream *snp=NULL
		, ostream *bioTagIndel=NULL
		, ostream *seqErr=NULL
		, ostream *splice=NULL
		, ostream *spliceIntra=NULL
		, ostream *spliceInter=NULL
		, ostream *undetermined=NULL
		, ostream *repetition=NULL
		, ostream *duplication=NULL
		, ostream *nothing=NULL
		, ostream *normal=NULL
		, ostream *almostNormal=NULL
		, ostream *multiple = NULL
		, ostream *none = NULL
		, ostream *bioUndermined = NULL
		, ostream *single = NULL
    , ostream *pairedEndChimera = NULL
    , samFile *sam_f = NULL
    , bam_hdr_t* sam_h = NULL
   );


   /** 
   * get headers for all ostreams   
   */
  void getHeaders(ostream *snp=NULL
		, ostream *bioTagIndel=NULL
		, ostream *seqErr=NULL
		, ostream *splice=NULL
		, ostream *spliceIntra=NULL
		, ostream *spliceInter=NULL
		, ostream *undetermined=NULL
		, ostream *repetition=NULL
		, ostream *duplication=NULL
		, ostream *nothing=NULL
		, ostream *normal=NULL
		, ostream *almostNormal=NULL
		, ostream *multiple = NULL
		, ostream *none = NULL
		, ostream *bioUndermined = NULL
		, ostream *single = NULL
    , ostream *pairedEndChimera = NULL);



  /**
   * Create a classifier object from a readIterator and increment it
   *
   * @readIt the iterator used to create the classifier
   * @isPairedEnd a boolean
   */
  Classifier *newClassifier(readIterator *readIt, bool isPairedEnd);

 private:

  /**
   * Set a read_t object with info from an iterator
   *
   * @readIt the iterator used to fill data of the read object
   * @read the read to fill
   */
  Read *newRead(readIterator *readIt);

  friend void *fillInfos(void *);

};

/**
 * Function called when creating a thread.
 * It fills a part of _store_tagsInfo array.
 */
void *fillInfos(void *);


typedef struct classify_thread_s {
  sem_t *remaining_seats;          /* Number of remaining cells for storing tagInfo */
  sem_t *processed_tags;
  Classifier **seats;               /* Seats where the Classifier are stored waiting to be classified*/
  uint total_seats;             /* Total number of seats */
  bool keepGoing;          /* Boolean to stop the tread */
  uint thread_num;
} classify_thread_t;


#endif
