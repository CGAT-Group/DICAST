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

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include "classifyTags.h"
#include "types.h"
#include "utils.h"
//#include <gkArrays.h>
//#include <gkArraysTypes.h>
#include "ReadIndex.h"
#include "libReadsReader/readsReader.h"
#include "libSSA/locateOnGenome.h"
#include "Parameters.h"
#include "libReadsInfo/tagInfo.h"
#include "Support.h"
#include "SupportBreak.h"
#include <semaphore.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>
#include <math.h>
#include <pthread.h>

using namespace std;
//using namespace gkarrays;

#define LINE_MAX_LENGTH 256
#define NB_POS_LOCATED 3
#define NB_STEP_LOCATION 7
#define SLEEP_DELAY_MS 1000

#define ARG_VALUE 0
#define ARG_CLASSIFY_TAGS 1
#define ARG_GENOME_SEMAPHORE 2
#define NB_ARGS_SUBSERVER 3

void dispLocateFactors(ostream & os, TagInfo *info, int pos,
		       LocateOnGenome *genome, 
		       uchar *tag, uint threshold, int increment
		       , uint max_done=NB_POS_LOCATED) {
  uint nb_done = 0;
  uchar *factor = new uchar[threshold+1];
  factor[threshold]=0;
  while (pos >= 0  && pos < (int)info->getSupportLength()
	 && nb_done < max_done) {
    strncpy((char *)factor, (char *)&tag[pos], threshold);
    if (info->getLocalisations()[pos] > 0) {
      os << "\t  pos. " << pos << ": ";
      genome->locateTag(factor, threshold);
    }
    pos += increment;
    nb_done++;
  }
  delete [] factor;
}

uint nb_digits(uint a) {
  return (uint)floor(log(a)/log(10))+1;
}

ostream &writeWithSpaces(ostream &os, uint writeIt, uint a, uint b) {
  uint nb_writeIt_digits = nb_digits(writeIt);
  uint max_nb_digits = max(nb_digits(a), nb_digits(b));
  max_nb_digits = max(max_nb_digits, nb_writeIt_digits);
  if (max_nb_digits > nb_writeIt_digits) {
    os << setw(max_nb_digits) << writeIt;
  } else {
    os << writeIt ;
  }
  return os;
}

void detailed_taginfo_output(ostream &os, uint tag_id, ClassifyTags *ct,
                             uint nb_pos_located, 
                             uint nb_step_location) {
  TagInfo *info;
  Classifier *classi;
  readIterator *readIt = ct->getReadsReader()->begin();
  bool is_paired = ct->getReadsReader()->isPairedEnd();

  // We incrment the iterator util we reach the read
  uint max = tag_id;
  // if we have paired-end reads and the query tag_id is the second member
  // of the paire, we want to stop at the first read of the pair in order to
  // build a proper PairedEndReadClassifier() !
  if(is_paired && tag_id % 2  != 0) {
    max--;
  }
  for(uint i = 0; i < max; i++) {
      ++(*readIt);
  }
    
  classi = ct->newClassifier(readIt,false);
  classi->init();
  classi->classify();
  if(!is_paired) {
    SingleReadClassifier *classifier = dynamic_cast<SingleReadClassifier*>(classi);
    info = classifier->getTagInfo();
  } else {
    PairedEndReadClassifier *classifier = dynamic_cast<PairedEndReadClassifier*>(classi);
    if (tag_id % 2 == 0) {
      info = classifier->getFirstTagInfo();
    } else {
      info = classifier->getSecondTagInfo();
    }
  }

  uchar *tag = (uchar *)info->getRead()->seq;
  LocateOnGenome *genome = ct->getGenomeIndex();
  ReadIndex *indexTags = ct->getReadIndex();
  uint tag_length = info->getRead()->getLength();
  uint threshold =  indexTags->getFactorLength();
  

  os << endl;
  os << "****** Tag " << tag_id << " ******" << endl;

  for (uint k = 1000; k >= 10; k/= 10) {
    if (k < tag_length) {
      for (uint i=0; i < tag_length; i++) {
        if (i/k > 0)
          os << (i/k) % 10 ;
        else
          os << " ";
      }
      os << endl;
    }
  }

  for (uint i=0; i < tag_length; i++) {
    os << i%10 ;
  }
  os << endl;
  os << tag << endl;

  os << endl;
  for (uint i= 0; i < info->getSupportLength(); i++) {
    writeWithSpaces(os,i,
                    info->getSupport()[i],
                    info->getLocalisations()[i]);
    os << " ";
  }
  os << endl;
  for (uint i= 0; i < info->getSupportLength(); i++) {
    writeWithSpaces(os,info->getSupport()[i],
                    i, info->getLocalisations()[i]);
    os << " ";
  }
  os << endl;
  for (uint i= 0; i < info->getSupportLength(); i++) {
    writeWithSpaces(os,info->getLocalisations()[i],
                    i, info->getSupport()[i]);
    os << " ";
  }
  os << endl;

  os << "Stats: " << info->getNbSingle() << " single located, "
              << info->getNbDuplicate() << " duplicated and " 
              << info->getNbMultiple() << " multiples" << endl;
  if (info->getStartPosRepeat() > -1) {
    os << "Repetition from position " << info->getStartPosRepeat();
    if (info->getEndPosRepeat() == -1)
      os << " until the end." << endl;
    else
      os << " to position " << info->getEndPosRepeat() << "." << endl;
  }
	  
  os << endl;
  if (info->getNbBreaks() > 0) {
    os << "Breaks: " << endl;
    for (uint i=0; i < info->getNbBreaks(); i++) {
      os << i << "\t" << info->getPositionStartBreak(i) 
	 << " -> "<< info->getPositionEndBreak(i);
      os << " (distance: " << info->getTagBreakLength(i);
	// 		    if (info->getIsExtended(i))
	// 		      os << " -- break extended";
      os  << ")";
      os << "\tscore (out/in): "
	 << info->getScoreOutsideBreak(i) << "/"
	 << info->getScoreInsideBreak(i);
      if (info->getScoreOutsideBreak(i) > 1) {
	os << ", score intra: " 
	   << info->getScoreComputedIntraExon(i)
	   << ", score inter: " 
	   << info->getScoreComputedInterExon(i)
	   << endl
	   << "\t\t\t\taverage (high/low): "
	   << info->getAverageHighInside(i)
	   << "/"
	   << info->getAverageLowInside(i)
	   << "\t score: "
	   << info->getScoreInsideAverages(i);
	if (info->getAverageLowInside(i) < 1 + EPSILON) {
	  os << endl << "\t\t\t\tfalling left: " 
	     << info->isSupportFallingLeft(i)
	     << "; right: " 
	     << info->isSupportFallingRight(i);
	}
      }
      os << endl;
      if (info->getNbCandidats(i) > 0 || !info->isRepeated(i)){
	os << "\tLocs: ";
	if (info->getPositionStartBreak(i) > 0)
	  os << *info->getChrPosStartBreak(i) << " (" 
	     << info->getLocationStartBreak(i) << ")";
	else
	  os << "??";
	os <<  " -> ";
	if (info->getPositionEndBreak(i) < info->getSupportLength()-1)
	  os << *info->getChrPosEndBreak(i) << " ("
	     << info->getLocationEndBreak(i) << ")";
	else
	  os << "??";
	if (info->getGenomeGapLength(i) != GAP_SIZE_MAX)
	  os << "\t gap length: " << info->getGenomeGapLength(i);
	os << endl;
	
	if (info->getPositionStartBreak(i) > 0) {
	  if (info->getPositionStartBreak(i) != info->getPositionStartBreak(i))
	    dispLocateFactors(os, info, info->getPositionStartBreak(i)-1,
			      genome,tag, threshold, 1, nb_pos_located);
	  else
	    dispLocateFactors(os, info, info->getPositionStartBreak(i)-1,
			      genome,tag, threshold, -1, nb_pos_located);
	  os << "\t-----" << endl;
	}
	if (info->getPositionEndBreak(i) < info->getSupportLength()-1) {
	  if (info->getPositionEndBreak(i) != info->getPositionEndBreak(i))
	    dispLocateFactors(os, info, info->getPositionEndBreak(i)+1,
			      genome, tag, threshold, -1, nb_pos_located);
	  else
	    dispLocateFactors(os, info, info->getPositionEndBreak(i)+1,
			      genome, tag, threshold, 1, nb_pos_located);
	}
      }else{
	os << "This break is not considered because there is a repeated kmer directly after it!" << endl;
      }
      os << "----------" << endl;
    }
  }

  if (info->hasNothing())
    os << "This tag has *not* been classified" << endl;
  else {
    os << "This tag is classified as:" << endl;
    if (info->isNormal()) {
      os << " Normal" << endl;
      dispLocateFactors(os, info, 0, genome, tag, 
                        threshold, info->getSupportLength() / nb_step_location,
                        nb_step_location+1);
    }
    if (info->isAlmostNormal()) {
      os << " Almost normal" << endl;
      dispLocateFactors(os, info, 0, genome, tag, 
                        threshold, info->getSupportLength() / nb_step_location,
                        nb_step_location+1);
    }
    if (info->isSingle())
      os << " Single" << endl;
    if (info->isMultiple())
      os << " Multiple" << endl;
    if (info->isNone()) 
      os << " None" << endl;
    if (info->hasRepetition())
      os << " Repetition" << endl;
    if (info->isDuplication()) {
      os << " Duplication" << endl;
      // 		    dispLocateFactors(os, info, 0, &genome, tag, threshold, 
      // 				      info->getSupportLength() / nb_step_location, 
      // 				      nb_step_location+1);
    } 
    if (info->hasSNP()) {
      os << " SNP (" << (uint) info->getNbSNP() << ")" << endl;
      for (uint j = 0; j < info->getNbSNP(); j++) {
	if (info->getInfosSNP()[j]->isDuplicated()){
	  os << "\tduplicate" ;
	}else{
	  os << "\tsingle" ;
	}
        os << " position " << info->getInfosSNP()[j]->getPosition()
                    << ", " << info->getInfosSNP()[j]->getExpectedNucleotide()
                    << " -> " << info->getInfosSNP()[j]->getActualNucleotide()
                    << "\t" << info->getInfosSNP()[j]->getChrPosition()
                    << "\t (score="<<info->getInfosSNP()[j]->getScore()<< ")"
                    << endl;
      }
    }
    if (info->hasSplice()) {
      os << " Splice (" << (uint) info->getNbSplice() << ")" << endl;
      for (uint j = 0; j < info->getNbSplice(); j++) {
        os << "\t" << *(info->getInfosSplice()[j]) << endl;
      }
    }
    if (info->getNbSpliceNoCover()) {
      os << " Coverless splice (" 
                  << (uint)info->getNbSpliceNoCover() << ")" << endl;
      for (uint j = 0; j < info->getNbSpliceNoCover(); j++) {
        os << "\t" << *(info->getInfosSpliceNoCover()[j])
                    << endl;		
      }
    }
    if (info->hasSpliceIntraChr()) {
      os << " Large splice (intra chromosome) (" 
                  << (uint)info->getNbSpliceIntra() << ")" << endl;
      for (uint j = 0; j < info->getNbSpliceIntra(); j++) {
	os << "\t" << *(info->getInfosSpliceIntra()[j])
                    << endl;		
      }
    }
    if (info-> hasSpliceInterChr()) {
      os << " Chimera (" << (uint) info->getNbSpliceInter() << ")" << endl;
      for (uint j = 0; j < info->getNbSpliceInter(); j++) {
	os << "\t" << *(info->getInfosSpliceInter()[j])
                    << endl;
      }
    }
    if (info->hasBioTagIndel()) {
      os << " Bio indel in the tag (" << (uint) info->getNbBioTagIndel()
                  << ")" << endl;
      for (uint j = 0; j < info->getNbBioTagIndel(); j++) {
	if (info->getInfosBioTagIndel()[j]->isDuplicated()){
	  os << "\tduplicate" ;
	}else{
	  os << "\tsingle" ;
	}
	os << " position " << info->getInfosBioTagIndel()[j]->getPosition()
                    << ": " << info->getInfosBioTagIndel()[j]->getNbIns()
                    << " insertion(s) and " 
                    << info->getInfosBioTagIndel()[j]->getNbDel() << " deletions"
                    << endl;
      }
    }
    if (info->hasBioUndetermined()) {
      os << " Undetermined biological reason (" 
                  << (uint)info->getNbBioUndetermined() << ")" << endl;
      for (uint j = 0; j < info->getNbBioUndetermined(); j++) {
        os << "\t position " << info->getInfosBioUndetermined()[j]->getPosition()
                    << ": " << info->getInfosBioUndetermined()[j]->getMessage()
                    << endl;
      }
    }
    if (info->hasSeqErr()) {
      os << " Sequence error (" << (uint) info->getNbSeqErr() << ")" 
                  << endl;
      for (uint j = 0; j < info->getNbSeqErr(); j++) {
	if (info->getInfosSeqErr()[j]->isDuplicated()){
	  os << "\tduplicate" ;
	}else{
	  os << "\tsingle" ;
	}
        os << " position " << info->getInfosSeqErr()[j]->getPosition();
        if (info->getInfosSeqErr()[j]->getNbDel() != 0) {
          if (info->getInfosSeqErr()[j]->getNbDel() == (uint)~0)
            os << " unknown number of deletions";
          else
            os << " " << info->getInfosSeqErr()[j]->getNbDel()
                        << " deletion(s)";
        }
        if (info->getInfosSeqErr()[j]->getNbIns() != 0) {
          if (info->getInfosSeqErr()[j]->getNbIns() == (uint)~0)
            os << " unknown number of insertions";
          else
            os << " " << info->getInfosSeqErr()[j]->getNbIns()
                        << " insertion(s)";
        }
        if (info->getInfosSeqErr()[j]->getNbIns() != (uint)~0
            && info->getInfosSeqErr()[j]->getNbDel() != (uint)~0) {
          os << "\t";
          if (info->getInfosSeqErr()[j]->getGenomeSequence()) 
            os << info->getInfosSeqErr()[j]->getGenomeSequence();
          else
            os << "?" ;
          os << " -> ";
          if (info->getInfosSeqErr()[j]->getErrorSequence()) 
            os << info->getInfosSeqErr()[j]->getErrorSequence();
          else 
            os << "?";
        }
        os << endl;
      }
    }
    if (info->hasUndeterminedError()) {
      os << " Undetermined error (" << (uint)info->getNbUndeterminedError() << ")" << endl;
      for (uint j = 0; j < info->getNbUndeterminedError(); j++) {
        os << "\t" << info->getInfosUndeterminedError()[j]->getMessage()
                    << endl;
      }
    }
    os << endl;
    info->samOutput(os);
    os << endl;


  } // end else hasNothing()

  delete classi;
  delete readIt;
}

void *launchSubServer(void *args) {
  void **theArgs = (void **)args;
  // Retrieving the arguments
  char *value = (char *)theArgs[ARG_VALUE];
  ClassifyTags *ct = (ClassifyTags*)theArgs[ARG_CLASSIFY_TAGS];
  sem_t genome_semaphore = * ((sem_t *)theArgs[ARG_GENOME_SEMAPHORE]);
  LocateOnGenome *genome = ct->getGenomeIndex();
  ReadIndex *indexTags = ct->getReadIndex();
  Parameters *params = ct->getParameters();
  delete [] theArgs;
  uint threshold = indexTags->getFactorLength();

  char *dontCare, *variable, *tag_length, *factor;
  char *end_ptr;
  uint tag_id;
  uint nb_pos_located = NB_POS_LOCATED, nb_step_location = NB_STEP_LOCATION;
  uint sleep_delay_us = SLEEP_DELAY_MS*1000;
  char line[LINE_MAX_LENGTH];
  bool quit = false;
  uint len_name = strlen(value)+9;
  char *fifo_inputname = new char[len_name];
  sprintf(fifo_inputname,"in.%s.fifo",value);
  char *fifo_outputname = new char[len_name+1];
  sprintf(fifo_outputname,"out.%s.fifo",value);

  if (mkfifo(fifo_inputname, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) {
    perror("Creation of fifo input file: ");
    exit(3);
  }
  if (mkfifo(fifo_outputname, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) {
    perror("Creation of fifo output file: ");
    exit(3);
  }

  FILE *fifo_input;
  fstream fifo_output;
  fifo_input = fopen(fifo_inputname, "r");

  while (! quit) {
    while (fgets(line, LINE_MAX_LENGTH, fifo_input) != NULL) {
      fifo_output.open((const char *)fifo_outputname,
		       fstream::out);
      fifo_output.precision(4);
      if (line[0] == 'q' && line[1] == 0) {
	fifo_output << "Quitting subserver" << endl;
	quit = true;
      } else {
	if (strncmp(line, "set|", 4) == 0) {
	  // Format set|<variable>|<value>
	  dontCare = strtok(line, "|");
	  variable = strtok(NULL, "|");
	  value = strtok(NULL, "|");
	  if (variable == NULL || value == NULL) {
	    fifo_output << "Error: bad format, should be set|<variable>|<value> "
			<< "(eg. set|nb_dispalyed|10)" << endl;
	  } else {
	    if (strncmp("nb_displayed", variable, 12) == 0) {
	      genome->setNbLocations(atol(value));
	      fifo_output << "nb_dispayed set to " << value << endl;
	    } else if (strncmp("nb_pos_located", variable, 14) == 0) {
	      nb_pos_located = atol(value);
	      fifo_output << "nb_pos_located set to " << value << endl;
	    } else if (strncmp("nb_step_location", variable, 16) == 0) {
	      nb_step_location = atol(value);
	      fifo_output << "nb_step_location set to " << value << endl;
	    } else if (strncmp("sleep_delay_ms", variable, 14) == 0) {
	      sleep_delay_us = atol(value)*1000;
	      fifo_output << "Now sleeping " << value << " ms." << endl;
            } else if (strncmp("deep_snp",variable, 8) == 0) {
              params->deep_snp_search ^= true;
              fifo_output << "SNP search ";
              if (! params->deep_snp_search) 
                fifo_output << "dis";
              fifo_output << "activated" << endl;
	    } else {
	      fifo_output << "Error: <variable> sould be one of these values: " 
			  << endl
			  << " nb_displayed, nb_pos_located, nb_step_location, sleep_delay_ms"
			  << endl
			  << " nb_displayed: number of locations displayed for a given position"
			  << endl
			  << " nb_pos_located: number of different positions located before and "
			  << endl
			  << "                 after a break" << endl
			  << " nb_step_location: incremental step when locating normal or duplicates"
			  << endl
			  << " sleep_delay_ms: Sleep the given amount of time (in ms) between two iterations"
			  << endl;
	    }
	  }
        } else if (strncmp("name|",line, 5) == 0) {
          // Strig in format name|<read name>
          uint length = strlen(line);
          uint name_length = length - 5;
          char *name = new char[name_length+1];
          strncpy(name, &line[5], name_length);
          name[name_length]=0;

          // Searching read whose name is name.
          readIterator *it = ct->getReadsReader()->begin();
          uint tag_id = 0;

          while (!it->isFinished()  && it->getName() 
                 && strcmp(name, it->getName())) {
            ++(*it);
            tag_id++;
          }
          if (!it->isFinished() && it->getName()) {
            genome->setOutputStreams (&fifo_output, &fifo_output, &fifo_output);
            // Name found
            detailed_taginfo_output(fifo_output, tag_id, ct,
                                    nb_pos_located, nb_step_location);
          } else {
            fifo_output << "Error: read name " << name << " was not found" 
                        << endl;
          }
          delete it;
          delete [] name;
	} else if (strncmp("tag|",line, 4) == 0) {
    fifo_output << "This command is not supported anymore.." << endl;
	  //// String in format tag|<factor length>|<factor>
	  //dontCare = strtok(line, "|");
	  //tag_length = strtok(NULL, "|");
	  //factor = strtok(NULL, "|");
	  //if (tag_length == NULL || factor == NULL) {
	  //  fifo_output << "Error: bad format, should be tag|<factor_length>|<factor> "
		//	<< "(eg. tag|13|ACAGTAGCATCAG" << endl;
	  //} else {
	  //  uint nb_result;
	  //  pair <uint, uint> *occs = indexTags->getTagsWithFactor(factor,
		//						   (uint) atol(tag_length),
		//						   nb_result);
	  //  fifo_output << nb_result << " Tags with " << factor << ": " << endl;
	  //  for (size_t i=0; i < nb_result; i++) {
	  //    fifo_output << occs[i].first << ":" << occs[i].second << "\t";
	  //    if ((i+1) % 5 == 0)
		//fifo_output << endl;
	  //  }	      
	  //  fifo_output << endl;
	  //}	    
	} else if (strncmp("GkSA|", line, 5) == 0) {
    fifo_output << "This command is not supported anymore.." << endl;
	  //dontCare = strtok(line, "|");
	  //value = strtok(NULL, "|");
	  //if (value == NULL) {
	  //  fifo_output << "Error: bad format, should be GkSA|<int>" << endl;
	  //} else  {
	  //  fifo_output << "GkSA[" << value << "]=" 
		//	<< indexTags->getGkSA((uintSA)atoll(value)) << endl;
	  //}
	} else if (strncmp("genome|", line, 7) == 0) {
	  // genome|<chr>|<pos>|<length>
	  dontCare = strtok(line, "|");
	  variable = strtok(NULL, "|");
	  value = strtok(NULL, "|");
	  tag_length = strtok(NULL, "|");
	  if ((variable == NULL) || (value == NULL) || (tag_length == NULL)) {
	    fifo_output << "Error: bad format, should be genome|<chr_name>|[-]<pos>|<length> (set chr_num to 0 for absolute position)" << endl;
	  } else {
	    long long int pos = atoll(value);
	    char *chr_name = variable;
	    int strand;
	    if (value[0] == '-') {
	      strand = -1;
	      pos = -pos;
	    } else 
	      strand = 1;		
	    uchar *dna = genome->getGenomeSubstring(pos, (ulong)atoll(tag_length), strand, chr_name);
	    fifo_output << dna << endl;
	    free(dna);
	  }
	} else {
    sem_wait(&genome_semaphore);
	  genome->setOutputStreams (&fifo_output, &fifo_output, &fifo_output);

	  if (strncmp(line, "loc|", 4) == 0) {
	    // String in format loc|<tag_length>|<factor>
	    dontCare = strtok(line, "|");
	    tag_length = strtok(NULL, "|");
	    factor = strtok(NULL, "|");
	    if (tag_length == NULL || factor == NULL) {
	      fifo_output << "Error: bad format, should be loc|<tag_length>|<factor> "
			  << "(eg. loc|13|ACAGTAGCATCAG" << endl;
	    } else {
	      genome->locateTag((uchar *)factor, atol(tag_length));
	    } // end else tag_length == NULL || factor == NULL
	  } else {
	    // String in format <tag id> or <factor>
	    tag_id = strtol(line, &end_ptr, 10);
            if (tag_id >= indexTags->getNbTags()) {
              fifo_output << "Error: " << tag_id << " is too large." << endl
                          << "There are only " << indexTags->getNbTags() 
                          << " reads indexed." << endl;
            } else if (*end_ptr != '\0') {
    fifo_output << "This command is not supported anymore.." << endl;
	  //    if (end_ptr == NULL || end_ptr[0] != ':') {
		//fifo_output << "Error: bad format, you must give the id of tag alone or "
		//	    << " <tag id>:<tag pos>" << endl;
	  //    } else {
		//uint tag_pos = strtol(&end_ptr[1], NULL, 10);
	  //    
		//uint nb_tags = indexTags->getNbTagsWithFactor(tag_id, tag_pos);
		//pair<uint, uint> *tags = indexTags->getTagsWithFactor(tag_id, tag_pos);
		//char *tag_factor = indexTags->getTagFactor(tag_id, tag_pos, threshold);
	  //    
		//if (nb_tags >= 2) {
		//  fifo_output << nb_tags << " tags contain the factor " << tag_factor << "." << endl;
		//  fifo_output << "Tags: ";
		//  for (uint i=0; i < nb_tags-1; i++) {
		//    fifo_output << tags[i].first << " pos. " << tags[i].second << "; ";
		//    if ((i+1)%10 == 0)
		//      fifo_output << endl;
		//  }
		//  fifo_output << tags[nb_tags-1].first << " pos. " << tags[nb_tags-1].second << "; ";
		//} else {
		//  fifo_output << "This is the only tag with the factor " << tag_factor << "." ;
		//}
		//fifo_output << endl;
		//delete [] tag_factor;
		//delete [] tags;	      
    //     }
	    } else {
	      //   for (tag_id = 32; tag_id < 33; tag_id++) {
              detailed_taginfo_output(fifo_output, tag_id, ct,
                                      nb_pos_located, nb_step_location);
	      //  }
	    } // end else end_ptr == line
	  } // else strcmp with loc|
    sem_post(&genome_semaphore);
	} // else -> locating on genome
      }
      fifo_output.close();
    }
    usleep(sleep_delay_us);
  }
  fclose(fifo_input);
  unlink(fifo_inputname);
  unlink(fifo_outputname);
  delete [] fifo_inputname;
  delete [] fifo_outputname;

  return NULL;
}

void locateServer(ClassifyTags *ct,
                  const char *fifo_name, const char *fifo_output_name) {
  char line[LINE_MAX_LENGTH];
  bool quit = false;

  if (mkfifo(fifo_name, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) {
    perror("Creation of fifo input file: ");
    exit(3);
  }
  if (mkfifo(fifo_output_name, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) {
    perror("Creation of fifo output file: ");
    exit(3);
  }

  cout << "Fifo created, let's go!" << endl;

  FILE *myfifo;
  myfifo = fopen(fifo_name, "r");

  char *dontCare;
  char *value=NULL;

  fstream out_fifo;
  
  uint sleep_delay_us = SLEEP_DELAY_MS*1000;
  sem_t genome_semaphore;
  sem_init(&genome_semaphore, 0, 1);

  while (! quit) {
    while (fgets(line, LINE_MAX_LENGTH, myfifo) != NULL) {
      out_fifo.open((const char *)fifo_output_name,
		    fstream::out);
      if (line[0] == 'q' && line[1] == 0 ) {
	out_fifo << "Quitting!" << endl;
	quit = true;
      } else if (strncmp(line, "connect@",8) == 0) {
	dontCare = strtok(line, "@");
	value = strtok(NULL, "@");
	
	pthread_t subServer;
	void **args = new void*[NB_ARGS_SUBSERVER];
	args[ARG_VALUE] = value;
	args[ARG_CLASSIFY_TAGS] = ct;
	args[ARG_GENOME_SEMAPHORE] = &genome_semaphore;
        
	pthread_create(&subServer, NULL, launchSubServer,args);
      } else {
	out_fifo << "Usage: connect@<id>" << endl;
      } // end else connect
      out_fifo.close();
    } // end while
    usleep(sleep_delay_us);
  }

  fclose(myfifo);
  sem_destroy(&genome_semaphore);

  // When server exits
  unlink(fifo_name);
  unlink(fifo_output_name);

}
