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

#include "classifyTags.h"
#include "utils.h"
#include <config.h>
#include <pthread.h>
#include <cstdlib>
#include <string>
#include <sstream>

#ifdef HAVE_LIBPROGRESSBAR
#include <sys/ioctl.h>
#include <libProgressBar/progressBar.h>
#endif

using namespace std;

ClassifyTags::ClassifyTags(LocateOnGenome *genome,
			   ReadIndex *tags,
         readsReader *reads_reader,
			   Parameters *param):
  genome(genome),tags(tags),reads_reader(reads_reader),parameters(param),nb_reads(0)
{
  nb_single = 0;
  nb_none = 0;
  nb_duplication = 0;
  nb_multiple = 0;
}

void ClassifyTags::classify(ostream *snp
			    , ostream *bioTagIndel
			    , ostream *seqErr 
			    , ostream *splice
			    , ostream *spliceNoCover
			    , ostream *spliceInter
			    , ostream *undetermined
			    , ostream *repetition
			    , ostream *duplication
			    , ostream *nothing
			    , ostream *normal
			    , ostream *almostNormal
			    , ostream *multiple
			    , ostream *none
			    , ostream *bioUndermined
			    , ostream *single
          , ostream *pairedEndChimera
          , samFile *sam_f
          , bam_hdr_t *sam_h
          ) {

  pthread_t threads[parameters->nb_threads]; // Threads objects
  classify_thread_t params[parameters->nb_threads]; // Threads parameters
  readIterator *readItForBuffer = reads_reader->begin(); // ReadIterator to iterate through the reads and feed the threads
  Classifier **store_classifiers[parameters->nb_threads-1]; // An array of classifier for each thread
  sem_t remaining_seats[parameters->nb_threads-1]; // Semaphore to control classifier waiting to be classified
  sem_t processed_tags[parameters->nb_threads-1]; // Semaphore to control classifier that has been classified
  uint *positions = NULL; // Array to track possition in the store_classifiers arrays
  uint nb_classifiers; // Number of element in each store_classifiers arrays
  bool isPairedEnd = reads_reader->isPairedEnd(); // Can tell if the readIterators are pairedEnd or not

  // Variable to avoid re-declaration
  Classifier *classifier;


#ifdef HAVE_LIBPROGRESSBAR
  DoccY::ProgressBar *PB = NULL;
  unsigned int maxval = 0;
  unsigned int cpt = 0;
  unsigned int reset_cpt = 0;
  if (parameters->show_progressbar) {
    unsigned int cols = 80;
#ifdef TIOCGSIZE
      struct ttysize ts;
      cols = ioctl(STDIN_FILENO, TIOCGSIZE, &ts) ? 80 : ts.ts_cols;    
#elif defined(TIOCGWINSZ)
      struct winsize ts;
      cols = ioctl(STDIN_FILENO, TIOCGWINSZ, &ts) ? 80 : ts.ws_col;
#endif /* TIOCGSIZE */

    maxval = cols ? cols : 100;
    reset_cpt = tags->getNbTags() / maxval;
    PB = new DoccY::ProgressBar("Processing Reads", maxval, cols, cerr, false);
    PB->ShowPercent();
    PB->ShowTime();
    PB->update();
  }
#endif

  // Init statistics
  for (uint i=0; i <= NB_MASKS; i++)
    nb_classes[i] = 0;

  nb_explainable=0;

  // Init threads
  if (parameters->nb_threads >  1) {
    nb_classifiers = parameters->nb_tags_info_stored;

    positions = new uint[parameters->nb_threads - 1];

    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      positions[j] = 0;
      sem_init(&processed_tags[j], 0, 0);
      params[j].processed_tags = &processed_tags[j];
      params[j].total_seats = nb_classifiers;
      if(isPairedEnd) {
        store_classifiers[j] = (Classifier**)calloc(nb_classifiers, sizeof(PairedEndReadClassifier*));
      } else {
        store_classifiers[j] = (Classifier**)calloc(nb_classifiers, sizeof(SingleReadClassifier*));
      }
      params[j].seats = store_classifiers[j];
      params[j].thread_num = j;
      params[j].keepGoing = true;
    }
  } else {
    nb_classifiers = 1;
  }

  // We give some food (reads to treat) to the threads
  if (parameters->nb_threads > 1) {
    // add reads data in threads
    uint i = 0;
    uint j = 0;
    while (!readItForBuffer->isFinished() && i < nb_classifiers) {
      j = 0;
      while (j < (parameters->nb_threads -1) && !readItForBuffer->isFinished()) {
        store_classifiers[j][i] = newClassifier(readItForBuffer,isPairedEnd);
	nb_reads++;
	
        j++;
      }
      i++;
    }
    // threads 0..j-1 have buffered i Classifiers
    // threads j..nb_threads-1 have buffered i-1 Classifiers
    
    // Init the semaphore according to the number of reads stored
    // in each thread
    for (uint k = 0; k < parameters->nb_threads -1; k++) {
      if (k >= j) {
        sem_init(&remaining_seats[k], 0, i-1);
        params[k].remaining_seats = &remaining_seats[k];
      } else {
        sem_init(&remaining_seats[k], 0, i);
        params[k].remaining_seats = &remaining_seats[k];
      }
    }

    // Launching threads
    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      pthread_create(&threads[j], NULL, fillInfos, &params[j]);
    }
  } 

  uint i = 0;
  while (i < nb_reads || !readItForBuffer->isFinished()) {
    if (parameters->nb_threads > 1) {
      // get the number of the thread we are waiting for
      uint j = i % (parameters->nb_threads - 1);
      // wait until the read(s) has been processed
      sem_wait(&processed_tags[j]);
      // get the Classifier object resulting of the read location
      classifier = store_classifiers[j][positions[j]];
      // Add a new read to the thread buffer
      if (!readItForBuffer->isFinished()) {
        params[j].seats[positions[j]] = newClassifier(readItForBuffer,isPairedEnd);
	nb_reads++;
        // Tell the tread that it can locate this new read
        sem_post(&remaining_seats[j]);
      }
      // Update the position of the next read to get from
      // this thread
      positions[j] = (positions[j] + 1) % nb_classifiers;
    } else {
      nb_reads++;
      classifier = newClassifier(readItForBuffer,isPairedEnd);
      classifier->init();
      classifier->classify();
    }
      
    classifier->writeOutputs(snp,
      bioTagIndel,
      seqErr,
      splice,
      spliceNoCover,
      spliceInter,
      undetermined,
      repetition,
      duplication,
      nothing,
      normal,
      almostNormal,
      multiple,
      none,
      bioUndermined,
      single);
    //classifier->samOutput(*sam);
    classifier->samOutput(sam_f,sam_h);
    if(isPairedEnd && pairedEndChimera != NULL)
      ((PairedEndReadClassifier*)classifier)->writePairedEndChimera(pairedEndChimera);
    classifier->updateStatistics(&nb_classes,&nb_explainable,&nb_single,&nb_duplication,&nb_multiple,&nb_none);
    delete classifier;

    i++;

#ifdef HAVE_LIBPROGRESSBAR
    // If we are using paired-end tags, increments loop ind by 2
    if(isPairedEnd) {
      cpt += 2;
    } else {
      cpt++;
    }
    if ((cpt >= reset_cpt) && parameters->show_progressbar) {
      PB->Step();
      cpt = 0;
    }
#endif
  }
#ifdef HAVE_LIBPROGRESSBAR
  if (parameters->show_progressbar) {
    PB->SetVal(maxval);
    PB->update(false);
    cerr << endl;
    delete PB;
  }
#endif
  
  // unlock threads
  for (uint j = 0; j < parameters->nb_threads -1; j++) {
    params[j].keepGoing = false;
    sem_post(&remaining_seats[j]);
  }

  // deleting the iterators, buffers and semaphores
  if (readItForBuffer != NULL) 
    delete readItForBuffer;
  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    free(store_classifiers[j]);
    sem_destroy(&remaining_seats[j]);
    sem_destroy(&processed_tags[j]);
  }
  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    pthread_join(threads[j], NULL);
  }

  if (positions)
    delete [] positions;

}

LocateOnGenome *ClassifyTags::getGenomeIndex() {
  return genome;
}

Parameters *ClassifyTags::getParameters() {
  return parameters;
}

ReadIndex *ClassifyTags::getReadIndex() {
  return tags;
}

readsReader *ClassifyTags::getReadsReader() {
  return reads_reader;
}

uint ClassifyTags::getNbReads() {
  return nb_reads;
}

uint ClassifyTags::getNbNone() {
  return nb_none;
  // return nb_classes[posBitInMask(MASK_NONE)];
}

uint ClassifyTags::getNbMultiple() {
  return nb_multiple;
  // return nb_classes[posBitInMask(MASK_MULTIPLE)];
}

uint ClassifyTags::getNbSingle() {
  return nb_single;
  // return nb_classes[posBitInMask(MASK_SINGLE)];
}

uint ClassifyTags::getNbAlmostNormal() {
  return nb_classes[posBitInMask(MASK_ALMOST_NORMAL)];
}

uint ClassifyTags::getNbExplainable() {
  return nb_explainable;
}

uint ClassifyTags::getNbNormal() {
  return nb_classes[posBitInMask(MASK_NORMAL)];
}

uint ClassifyTags::getNbDuplication() {
  return nb_duplication;
  // return nb_classes[posBitInMask(MASK_DUPLICATION)];
}

uint ClassifyTags::getNbSNP() {
  return nb_classes[posBitInMask(MASK_SNP)];
}

uint ClassifyTags::getNbBioTagIndel() {
  return nb_classes[posBitInMask(MASK_BIOLOGICAL_TAG_INDEL)];
}

uint ClassifyTags::getNbBioUndetermined() {
  return nb_classes[posBitInMask(MASK_BIOLOGICAL_UNDETERMINED)];
}

uint ClassifyTags::getNbSeqErr() {
  return nb_classes[posBitInMask(MASK_SEQ_ERR)];
}


uint ClassifyTags::getNbRepetition() {
  return nb_classes[posBitInMask(MASK_REPETITION)];
}

uint ClassifyTags::getNbSplice() {
  return nb_classes[posBitInMask(MASK_SPLICE)];
}

uint ClassifyTags::getNbSpliceIntra() {
  return nb_classes[posBitInMask(MASK_INTRA_TRANSPLICING)];
}

uint ClassifyTags::getNbSpliceInter() {
  return nb_classes[posBitInMask(MASK_INTER_TRANSPLICING)];
}

uint ClassifyTags::getNbSpliceNoCover() {
  return nb_classes[posBitInMask(MASK_SPLICE_NO_COVER)];
}

uint ClassifyTags::getNbNothing() {
  return nb_classes[NB_MASKS];
}

uint ClassifyTags::getNbUndetermined() {
  return nb_classes[posBitInMask(MASK_UNDETERMINED_ERROR)];
}

uint ClassifyTags::getNbPairedEndChimera() {
  return nb_classes[posBitInMask(MASK_PAIRED_END_CHIMERA)];
}

void ClassifyTags::getHeaders(ostream *snp
			    , ostream *bioTagIndel
			    , ostream *seqErr 
			    , ostream *splice
			    , ostream *spliceNoCover
			    , ostream *spliceInter
			    , ostream *undetermined
			    , ostream *repetition
			    , ostream *duplication
			    , ostream *nothing
			    , ostream *normal
			    , ostream *almostNormal
			    , ostream *multiple
			    , ostream *none
			    , ostream *bioUndermined
			    , ostream *single
          , ostream *pairedEndChimera) {

  stringstream headerTemplate;
  headerTemplate << "#CRAC version " << PACKAGE_VERSION << endl;
  headerTemplate << "#Created on " <<__DATE__ << endl << "#" << endl;

  string chrPosPatternInfo = "#A pattern [ chr|strand,pos ] is an occurrence of location on the index (genome or reference)";
  
  // For SNV
  if (snp){
    *snp << headerTemplate.str() << chrPosPatternInfo << endl;
    *snp << "#A tag \"single\" is noted if there is one possibility of SNV, a tag \"duplicate\" is noted otherwise" << endl;
    *snp << "#read_id tag_snv loc_snv_on_genome pos_snv_on_read snv crac_score single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For indels
  if (bioTagIndel){
    *bioTagIndel << headerTemplate.str() << chrPosPatternInfo << endl;
    *bioTagIndel << "#A tag \"single\" is noted if there is one possibility of indel, a tag \"duplicate\" is noted otherwise" << endl;
    *bioTagIndel << "#read_id tag_indel loc_indel_on_genome pos_start_indel_on_read nb_of_insertion nb_of_deletion single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For error
  if (seqErr){
    *seqErr << headerTemplate.str() << chrPosPatternInfo << endl;
    *seqErr << "#A tag \"single\" is noted if there is one possibility of error, a tag \"duplicate\" is noted otherwise" << endl;
    *seqErr << "#Two different formats for error: "<< endl;
    *seqErr << "#read_id tag_error loc_error_on_genome pos_ponctual_error_on_read error crac_score single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
    *seqErr << "#or" << endl;
    *seqErr << "#read_id tag_error loc_error_on_genome pos_start_indel_error_on_read nb_of_insertion nb_of_deletion single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For splice and splice_no_cover
  if (splice){
    *splice << headerTemplate.str() << chrPosPatternInfo << endl;
    *splice << "#A tag \"single\" is noted if there is one possibility of splice, a tag \"duplicate\" is noted otherwise" << endl;
    *splice << "#read_id tag_splice loc_end_first_exon_on_genome pos_junction_on_read splice_length single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }

  if (spliceNoCover) {
    *spliceNoCover << headerTemplate.str() << chrPosPatternInfo << endl;
    *spliceNoCover << "#A tag \"single\" is noted if there is one possibility of splice, a tag \"duplicate\" is noted otherwise" << endl;
    *spliceNoCover << "#read_id tag_splice loc_end_first_exon_on_genome pos_junction_on_read splice_length single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For chimera
  if (spliceInter){
    *spliceInter << headerTemplate.str() << chrPosPatternInfo << endl;
    *spliceInter << "#A tag \"single\" is noted if there is one possibility of chimera, a tag \"duplicate\" is noted otherwise" << endl;
    *spliceInter << "#read_id tag_chimera chimera_class chimera_score loc_end_first_exon_on_genome loc_start_second_exon_on_genome pos_junction_on_read single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }

  if (pairedEndChimera){
    *pairedEndChimera << headerTemplate.str() << chrPosPatternInfo << endl;
    *pairedEndChimera << "#A paire of paired-end tags classified as \"single\" are noted if there is one possibility of chimera between them" << endl;
    *pairedEndChimera << "#chimera_class read_id1 single_loc_on_genome1(private) pos_single_loc_on_read1(private) read1 p_support1 p_loc1 read_id2 single_loc_on_genome2(private) pos_single_loc_on_read2(private) read2 p_support2 p_loc2" << endl;
  }

  // For undetermined cause and bioUndermined
  if (undetermined){
    *undetermined << headerTemplate.str() << chrPosPatternInfo << endl;
    *undetermined << "#read_id undetermined_cause_features single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (bioUndermined) {
    *bioUndermined << headerTemplate.str() << chrPosPatternInfo << endl;
    *bioUndermined << "#read_id bioUndetermined_cause_features single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For repetition
  if (repetition){
    *repetition << headerTemplate.str() << chrPosPatternInfo << endl;
    *repetition << "#read_id pos_start_repeat_on_read pos_end_repeat_on_read single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For duplication and multiple
  if (duplication){
    *duplication << headerTemplate.str() << chrPosPatternInfo << endl;
    *duplication << "#read_id occurrence_loc_on_genome(private) pos_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (multiple) {
    *multiple << headerTemplate.str() << chrPosPatternInfo << endl;
    *multiple << "#read_id one_occurrence_loc_on_genome(private) pos_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For nothing
  if (nothing){
    *nothing << headerTemplate.str();
    *nothing << "#read_id read p_support p_loc" << endl;
  }
  // For normal and almostNormal
  if (normal){
    *normal << headerTemplate.str() << chrPosPatternInfo << endl;
    *normal << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (almostNormal) {
    *almostNormal << headerTemplate.str() << chrPosPatternInfo << endl;
    *almostNormal << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For none
  if (none){
    *none << headerTemplate.str() << chrPosPatternInfo << endl;
    *none << "#read_id [loc_on_genome(private) pos_loc_on_read(private)] read p_support p_loc" << endl;
  }
  // For single
  if (single){
    *single << headerTemplate.str() << chrPosPatternInfo << endl;
    *single << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
}

Classifier *ClassifyTags::newClassifier(readIterator *readIt, bool isPairedEnd) {
  Classifier *classifier;
  if (isPairedEnd) {
    Read *r1, *r2;
    r1 = newRead(readIt);
    ++(*readIt);
    r2 = newRead(readIt);
    ++(*readIt);
    classifier = new PairedEndReadClassifier(r1,r2,getGenomeIndex(),getReadIndex(),getParameters());
  } else {
    Read *r = newRead(readIt);
    ++(*readIt);
    classifier = new SingleReadClassifier(r,getGenomeIndex(),getReadIndex(),getParameters());
  }
  return classifier;
}

/// PRIVATE ///

Read *ClassifyTags::newRead(readIterator *it) {
  Read *r;
  if(parameters->reads_length > 0) {
    r = new Read(it->getReadNumber(), it->getSequence(), it->getName(), it->getQuality(), parameters->reads_length);
  } else {
    r = new Read(it->getReadNumber(), it->getSequence(), it->getName(), it->getQuality());
  }
  return r;
}

void *fillInfos(void *args) {
  classify_thread_t params = ((classify_thread_t *)args)[0];
  uint current_pos = 0;
  while(((classify_thread_t *)args)[0].keepGoing) {
    // We need one seat for the current read
    sem_wait(params.remaining_seats);

    // This supplementary test is a security for the main thread to unlock
    // this thread
    if(((classify_thread_t *)args)[0].keepGoing) {
      // init the classifier attributes (support profile, loc profile)
      params.seats[current_pos]->init();
      // Run the analysis of the read's alignement
      params.seats[current_pos]->classify();
      // Tell the main thread that we have classified a read
      sem_post(params.processed_tags);
      // Update the position for the next classifier to process
      current_pos = (current_pos + 1) % params.total_seats;
    }
  }
  pthread_exit(NULL);
  return NULL;
}
