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
#include "locateTags.h"
#include "utils.h"
#include <config.h>
#include <sstream>
#include <algorithm>

using namespace std;

LocateTags::LocateTags(LocateOnGenome *genome, Parameters *param, char *tags_file, uint threshold):
  genome(genome),parameters(param),tags_file(tags_file),threshold(threshold)
{}


void LocateTags::locate(samFile *sam,bam_hdr_t* sam_h) {
  // Init stat variables
  nb_single = 0;
  nb_multiple = 0;
  nb_none = 0;
  nb_tags = 0;
  stats_mutex = sem_open("/stats_mutex", O_CREAT, 0644, 1);
  sem_unlink("/stats_mutex");

  pthread_t threads[parameters->nb_threads];
  locate_thread_t params[parameters->nb_threads];
  vector<SamLine*> **store_samLines[parameters->nb_threads-1];
  read_t **store_reads[parameters->nb_threads-1];
  readsReader reads(tags_file,threshold);
  readIterator *readItForBuffer = reads.begin();
  readIterator *readIt = reads.begin();
  sem_t *remaining_seats[parameters->nb_threads-1];
  sem_t *processed_reads[parameters->nb_threads-1];
  uint *positions = NULL;
  uint nb_samLines;
  vector<SamLine*> *sam_lines;
  string sem_name;

  if (parameters->nb_threads >  1) {
    nb_samLines = parameters->nb_tags_info_stored;

    positions = new uint[parameters->nb_threads - 1];

    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      positions[j] = 0;
      params[j].locate = this;
      sem_name = concatStringAndInt("/processed_reads",j);
      processed_reads[j] = sem_open(sem_name.c_str(), O_CREAT, 0644, 0);
      if(processed_reads[j] == SEM_FAILED) 
        cerr << "Semaphore named : " << sem_name << " failed to open." << endl;
      sem_unlink(sem_name.c_str());
      params[j].processed_reads = processed_reads[j];
      params[j].total_seats = nb_samLines;
      store_reads[j] = (read_t **)malloc(nb_samLines*sizeof(read_t *));
      params[j].reads = store_reads[j];
      store_samLines[j] = (vector<SamLine*> **)malloc(nb_samLines*sizeof(vector<SamLine*> *));
      params[j].seats = store_samLines[j];
      params[j].keepGoing = true;
    }
  } else {
    nb_samLines = 1;
  }

  if (parameters->nb_threads > 1) {
    // add reads data in threads
    uint i = 0;
    uint j = 0;
    uint readCpt = 0;
    while (!readItForBuffer->isFinished() && i < nb_samLines) {
      j = 0;
      while (j < (parameters->nb_threads -1) && !readItForBuffer->isFinished()) {
        read_t *r = new read_t();
        setReadInfo(readItForBuffer,r);
        store_reads[j][i] = r;
        ++(*readItForBuffer);
        readCpt++;
        j++;
      }
      i++;
    }
    // Init the semaphore according to the number of reads stored
    // in each thread
    for (uint k = 0; k < parameters->nb_threads -1; k++) {
      sem_name = concatStringAndInt("/remaining_seats",k);
      if (k >= j) {
        remaining_seats[k] = sem_open(sem_name.c_str(), O_CREAT, 0644, i-1);
        if(remaining_seats[k] == SEM_FAILED) 
          cerr << "Semaphore named : " << sem_name << " failed to open." << endl;
        sem_unlink(sem_name.c_str());
        params[k].remaining_seats = remaining_seats[k];
      } else {
        remaining_seats[k] = sem_open(sem_name.c_str(), O_CREAT, 0644, i);
        if(remaining_seats[k] == SEM_FAILED) 
          cerr << "Semaphore named : " << sem_name << " failed to open." << endl;
        sem_unlink(sem_name.c_str());
        params[k].remaining_seats = remaining_seats[k];
      }
    }

    // Launching threads
    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      pthread_create(&threads[j], NULL, fillLocations, &params[j]);
    }
  } 

  uint i = 0;
  while (!readIt->isFinished()) {
    if (parameters->nb_threads > 1) {
      // get the number of the thread we are waiting for
      uint j = i % (parameters->nb_threads - 1);
      // wait until the read has been processed
      sem_wait(processed_reads[j]);
      // get the sam line resulting of the read location
      sam_lines = store_samLines[j][positions[j]];
      // delete the read object that has been used to get the location
      delete store_reads[j][positions[j]];
      // Add a new read to the thread buffer
      if (!readItForBuffer->isFinished()) {
        read_t *r = new read_t();
        params[j].reads[positions[j]] = r;
        setReadInfo(readItForBuffer,params[j].reads[positions[j]]);
        ++(*readItForBuffer);
        // Tell the tread that it can locate this new read
        sem_post(remaining_seats[j]);
      }
      // Update the position of the next read to get from
      // this thread
      positions[j] = (positions[j] + 1) % nb_samLines;
    } else {
      // If we dont use multi-threading we directly compute
      // sam line for the current read
      read_t *read = new read_t();
      setReadInfo(readIt,read);
      sam_lines = locateTag(read);
      delete read;
    }
    // Write the sam line in the file
    //*sam << *sam_line;
    // FIXME
    for(vector<SamLine*>::iterator it = sam_lines->begin(); it != sam_lines->end(); ++it) {
      (*it)->writeBamRecord(sam,sam_h);
      delete *it;
    }
    delete sam_lines;
    ++(*readIt);
    i++;
    nb_tags++;
  }

  // unlock threads
  for (uint j = 0; j < parameters->nb_threads -1; j++) {
    params[j].keepGoing = false;
    sem_post(remaining_seats[j]);
  }
  
  // deleting iterators
  delete readIt;
  delete readItForBuffer;

  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    free(store_samLines[j]);
    free(params[j].reads);
    sem_close(processed_reads[j]);
    sem_close(remaining_seats[j]);
  }

  // Closing stats_mutex
  sem_close(stats_mutex);

  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    pthread_join(threads[j], NULL);
  }

  if (positions)
    delete [] positions;

}

vector<SamLine*> *LocateTags::locateTag(read_t *read) {
  ostringstream os;
  pair<ChrPosition **,uint>locs;
  uint pos_locs;
  vector<SamLine*> *sam_lines = new vector<SamLine*>();
  SamLine *base_line = new SamLine();

  // minuscule nucleotide to capital
  string tag_s = read->seq; 
  ulong tag_length = (ulong)tag_s.size();
  transform(tag_s.begin(), tag_s.end(), tag_s.begin(), static_cast<int(*)(int)>(toupper));
  uchar *tag = (uchar *)tag_s.c_str();
  //string *samLines = new string();

  if (threshold == 0){
    locs = genome->getLocations(&tag[0], tag_length, tag_length, pos_locs);
  }else{
    locs = genome->getLocations(&tag[0], threshold, tag_length, pos_locs);
  }

  // write generic SAM features
  base_line->setQname(read->name);
  base_line->setSeq(read->seq);
  base_line->setQual(read->qual);

  // write sam_line, if no loc
  if (locs.second == 0) {
    base_line->setSegmentUnmapped();
    base_line->setMapQ(255);
    base_line->addOptionalField("NH",locs.second);
    sam_lines->push_back(base_line);
    //base_line->writeLine(os);
    // set mutex for none stats
    sem_wait(stats_mutex);
    nb_none++;
    sem_post(stats_mutex);
    // end mutex zone
  }
  // write sam_line for each loc (or a single loc)
  else{    
    for(ulong i = 0; i < locs.second; i++) {
      if (i==0 || parameters->treat_multiple){
        // We use the copy constructor to get a fresh sam_line that is based
        // on our baseline
        SamLine *current_line = new SamLine(*base_line);
        string rname = string(locs.first[i]->getChrPosition());
        current_line->setRname(rname,locs.first[i]->getChrId());
        // +1 because SAM is 1-based
        current_line->setPos(locs.first[i]->getRelativePosition() + 1);
        current_line->setMapQ(254);

        Cigar cigar;
        if (threshold > 0){
          uint soft_size = tag_length-threshold-pos_locs;
          if (locs.first[i]->getStrand() == 1){
            if (pos_locs > 0)
              cigar.append(pos_locs,'S');
            cigar.append(threshold,'M');
            if (soft_size  > 0)
              cigar.append(soft_size,'S');
          }else{
            if (soft_size > 0)
              cigar.append(soft_size,'S');
            cigar.append(threshold,'M');
            if (pos_locs > 0)
              cigar.append(pos_locs,'S');
          }
        }else
          cigar.append(tag_length,'M');
        current_line->setCigar(cigar);
        if (locs.first[i]->getStrand() == -1){
          current_line->setSeqReverseComplemented();
          current_line->reverseComplementeSeq();
          current_line->reverseQual();
        }
        //multiple_alignement.writeLine(cout);
        // If this is not the primary line
        if(i>0) {
          current_line->setSecondaryAlignement();
        }
        current_line->addOptionalField("NH",locs.second);
        // set single and multiple and update stats
        // use a mutex to protect globa stats update
        sem_wait(stats_mutex);
        if (locs.second == 1)
          nb_single++;
        else 
          nb_multiple++;
        sem_post(stats_mutex);
        // end of mutex zone
        //current_line->writeLine(os);
        sam_lines->push_back(current_line);
      }
      delete locs.first[i];
    }
    free(locs.first);
    delete base_line;
  }
  //*samLines = os.str();  
  //delete [] tag;
  return sam_lines;
}

uint LocateTags::getNbTags(){
  return nb_tags;
}

uint LocateTags::getNbSingle(){
  return nb_single;
}

uint LocateTags::getNbMultiple(){
  return nb_multiple;
}

uint LocateTags::getNbNone(){
  return nb_none;
}

// private methods

void LocateTags::setReadInfo(readIterator *readIt,read_t *read) {
  if(readIt->getName() != NULL) {
    read->name = readIt->getName();
  } else { read->name = "*"; }
  read->seq = readIt->getSequence();
  if(readIt->getQuality() != NULL) {
    read->qual = readIt->getQuality();
  } else { read->qual = "*";}
}

void *fillLocations(void *args) {
  locate_thread_t params = ((locate_thread_t *)args)[0];
  LocateTags *locate = params.locate;
  uint current_pos = 0;
  while(((locate_thread_t *)args)[0].keepGoing) {
    // We need one seat for the current read
    sem_wait(params.remaining_seats);

    if(((locate_thread_t *)args)[0].keepGoing) {
      params.seats[current_pos] = locate->locateTag(params.reads[current_pos]);
      sem_post(params.processed_reads);
      // delete the read
      //delete params.reads[current_pos];
      current_pos = (current_pos + 1) % params.total_seats;
    }
  }
  pthread_exit(NULL);
  return NULL;
}
