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

#include <stdlib.h>
#include <utility>
#include <vector>
#include <algorithm>
#include "Support.h"
#include "BreakList.h"
#include "utils.h"
#include "Bitset.h"

Support::Support(Parameters *p, uint *s, Read *r
		,LocateOnGenome *g, ReadIndex *ind, Support *psupport):
  parameters(p), support(s), read(r), length(r->getLength()-ind->getFactorLength()+1),
  genome(g), tags(ind), pair_support(psupport), nb_pos_located(0) {

  
  ulong nb_occs_sens, nb_occs_antisens;
//   ulong *occs_sens, *occs_antisens;
  ulong total_occ;
    bool located = true;

  ulong max_nb_breaks = NB_BREAKS_ALLOCATION_IN_SUPPORT;
  breaks = (SupportBreak**)malloc(sizeof(SupportBreak*)*max_nb_breaks);
  ulong pos_break=0;
  
//   pos_occ = new uint[tag_length-threshold+1];
//   strand_occ = new int[tag_length-threshold+1];

  nb_locs = new uint[length];

  start_pos_repeat = -1;
  end_pos_repeat = -1;
  nb_breaks = 0;
  nb_single = 0;
  nb_multiple = 0;
  nb_duplicate = 0;
  almostNormal = false;
  nb_locs_max = 0;
  position_of_location = read->getLength();
  best_location = NULL;
  locations = NULL;
  
  ranges_forward = new pair<uint, uint>[length];
  ranges_reverse = new pair<uint, uint>[length];

  average_support = 0;
  for (ulong j=0; j < length; j++) {
    ranges_forward[j] = genome->getFMIndexRange((uchar *)&read->seq[j],ind->getFactorLength());
    nb_occs_sens = (ranges_forward[j].first <= ranges_forward[j].second) ? 
      (ranges_forward[j].second - ranges_forward[j].first + 1):0;
    
    ranges_reverse[j] = genome->getFMIndexReverseRange((uchar *)&read->seq[j],ind->getFactorLength());
    nb_occs_antisens = (ranges_reverse[j].first <= ranges_reverse[j].second) ? 
      (ranges_reverse[j].second - ranges_reverse[j].first + 1):0;
    nb_locs[j] = nb_occs_antisens + nb_occs_sens;

    //if support contains 'N' associated k-mer is not in the index and the support is equal to 0
    //TODO: associate a bioUndetermined error when support[j] == 0
    if (support[j] == 0){
      support[j] = 1;
    }

    average_support += support[j];    
  }
  average_support /= length;
  
  // process to choose the corresponding k-mer for read location
  // cerr << "First compute single" << endl;
  computeBestLocation(true);
  if (position_of_location == read->getLength() && pair_support && pair_support->getNbLocsMax() == 1){
    // cerr << "First compute all" << endl;
    computeBestLocation(false);
  }

  if (position_of_location != read->getLength()){
    setLocationsAt(position_of_location);
    if (getNbLocsMax()>0){
      if (pair_support
	  && getPairSupport()->getNbLocsMax() > 0
	  && ((getPairSupport()->getNbLocsMax() <= parameters->max_mapping_comparison && getNbLocsMax() <= parameters->max_localisation_duplication) 
	      || (getNbLocsMax() <= parameters->max_mapping_comparison && getPairSupport()->getNbLocsMax() <= parameters->max_localisation_duplication)
	      || getNbLocsMax() == 1
	      || getPairSupport()->getNbLocsMax() == 1)
	  ){
	// Only extract intersection of locations between the two paired reads
	checkPairedConnection();
      }
      // if !paired or nb_loc == 1, so the best location is the first element of locations
      if (pair_support == NULL || (getNbLocsMax() == 1 && getBestLocation() == NULL)){
      	setBestLocation(locations[0]);
      }
    }
  }
  
  
  // temporary repeat start, end and length
  int minRepeat = 0;
  int start_pos_r = -1;
  int end_pos_r = -1;
  for (ulong j=0; j < length; j++) {    
    total_occ = nb_locs[j];
    
    // almost normal
    // TODO: check the condition (it seems strange)
    // if (!almostNormal && (j+1) < length) {
    //   almostNormal =  ((support[j] > support[j+1] 
    // 			  && support[j]*(1.0 - parameters->percent_support_variation_almost_normal) > support[j+1]) 
    // 			 || (support[j] < support[j+1] && support[j]*1/(1.0 - parameters->percent_support_variation_almost_normal) < support[j+1]));    
    // }    
    if (total_occ > 0) {
      nb_pos_located++;
      // end of repetition in case of repetition (we suppose one
      // repetition per tag)
      if (total_occ < parameters->min_occ_repetition 
	  && start_pos_repeat != -1 && end_pos_repeat == -1 ){
	end_pos_repeat = j-1;
      }
      // only the longest repetition is saved
      if (end_pos_repeat == (int)(j-1)) {
	if (getRepeatLength() > minRepeat){
	  minRepeat = getRepeatLength();
	  start_pos_r = start_pos_repeat;
	  end_pos_r = end_pos_repeat;
	}
	start_pos_repeat = -1;
	end_pos_repeat = -1;
      }
      
      // duplication
      if (total_occ <= parameters->max_localisation_duplication) {
	
	if (total_occ >= parameters->min_localisation_duplication
	    && total_occ <= parameters->max_localisation_duplication) {
	  nb_duplicate++;
	}
      } else {
     	nb_multiple++;
	if(start_pos_repeat == -1 
	   && total_occ >= parameters->min_occ_repetition){
	  start_pos_repeat = j;
	}
      }
          
      if (total_occ == 1){
	nb_single++;
      }
      
      // We didn't attribute the score of the previous break 
      // and we're yet quite far from the break, so we attribute it right now!
      if (! located) {
	// We have more breaks than expected.
	// We expand the array.
	if (nb_breaks > max_nb_breaks) {
	  max_nb_breaks *= 2;
	  breaks = (SupportBreak**)realloc(breaks, sizeof(SupportBreak*)*max_nb_breaks);
	}

	  breaks[nb_breaks-1] = new SupportBreak(pos_break, j-1, 
						 this, parameters,ranges_forward,
						 ranges_reverse);
      }

      located = true;

    }else{
      if (located) {
	nb_breaks++;
	pos_break = j;
      }
      located = false;
    } 
  }
  
  // save the biggest repetition
  if (getRepeatLength() <= minRepeat){
    start_pos_repeat = start_pos_r;
    end_pos_repeat = end_pos_r;
  }
    
  if (! located) {
    if (nb_pos_located) {
      // We have more breaks than expected.
      // We expand the array.
      if (nb_breaks > max_nb_breaks) {
	max_nb_breaks *= 2;
	breaks = (SupportBreak**)realloc(breaks, sizeof(SupportBreak*)*max_nb_breaks);
      }
      breaks[nb_breaks-1] = new SupportBreak(pos_break, length-1,
					     this, parameters,ranges_forward,
					     ranges_reverse);
    } else {
      nb_breaks--;
    }
  } 

  
  // Now we try to merge breaks
  if (nb_breaks > 1) {
    tryToMergeBreaks();
  } else if (nb_breaks== 1) {
    breaks[0]->adjustBreakWithSplicingSites();
  }


  if (getBestLocation() == NULL
      // && (!pair_support || pair_support->getNbLocsMax() != 1)
      ){
    // cerr << "Second compute all" << endl;
    // A last computation (more softly) if we did not find a good one before
    bool before_found = true;
    if (position_of_location == read->getLength()){
      computeBestLocation(false,true);
      before_found = false;
    }
    
    // There is at least a k-mer located
    if (position_of_location != read->getLength()){
      if (!before_found){
	setLocationsAt(position_of_location);
	if (getNbLocsMax()>0){
	  if (pair_support
	      && getPairSupport()->getNbLocsMax() > 0
	      && (
		  // (getPairSupport()->getNbLocsMax() <= parameters->max_mapping_comparison && getNbLocsMax() <= parameters->min_localisation_duplication) 
		  // || (getNbLocsMax() <= parameters->max_mapping_comparison && getPairSupport()->getNbLocsMax() <= parameters->min_localisation_duplication)
		  getNbLocsMax() == 1
		  || getPairSupport()->getNbLocsMax() == 1)	    
	      ){
	    // Only extract intersection of locations between the two paired reads
	    checkPairedConnection();
	  }
	}
      }
      // In the worst case, the best location is the first element of locations
      if (getBestLocation()==NULL && getNbLocsMax()>0){
        setBestLocation(locations[0]);
      }
    }
  }
  
  almostNormal = almostNormal && ! hasRepeat() && isContinuous();
}

Support::Support(const Support &s): parameters(s.parameters), 
                                    read(s.read),
                                    length(s.length),
                                    genome(s.genome),
                                    tags(s.tags), 
                                    pair_support(s.pair_support),
                                    nb_pos_located(s.nb_pos_located),
                                    nb_single(s.nb_single),
                                    nb_multiple(s.nb_multiple), nb_duplicate(s.nb_duplicate),
                                    almostNormal(s.almostNormal),
                                    start_pos_repeat(s.start_pos_repeat),
                                    end_pos_repeat(s.end_pos_repeat), nb_breaks(s.nb_breaks)
                           
{
  support = new uint[length];
  nb_locs = new uint[length];
  ranges_reverse = new pair<uint,uint>[length];
  ranges_forward = new pair<uint,uint>[length];
  for (uint i = 0; i < length; i++) {
    support[i] = s.support[i];
    nb_locs[i] = s.nb_locs[i];
    ranges_forward[i] = s.ranges_forward[i];
    ranges_reverse[i] = s.ranges_reverse[i];
  }

  if (s.best_location) {
    best_location = new ChrPosition(*s.best_location);
  }
  
  if (s.locations != NULL) {
    locations = (ChrPosition **)malloc(sizeof(ChrPosition *) * s.nb_locs_max);
    
    for (uint i=0; i < s.nb_locs_max ; i++){
      locations[i] = new ChrPosition(*(s.locations[i]));
    }
    
    nb_locs_max = s.nb_locs_max;
    position_of_location = s.position_of_location;
  }
  
  breaks = new SupportBreak*[nb_breaks];
  for (uint i=0; i < nb_breaks; i++) {
    breaks[i] = new SupportBreak(*(s.breaks[i])); 
    breaks[i]->setSupport(this);
    breaks[i]->setRangesForward(ranges_forward);
    breaks[i]->setRangesReverse(ranges_reverse);
  }
  
}

Support::~Support() {
  for (uint i=0; i < nb_breaks; i++) {
    delete breaks[i];
  }
  free(breaks);
  
  if (locations){
    for (uint i=0 ; i < nb_locs_max ; i++){
      delete locations[i];
    }
    free(locations);
  }
  if (best_location)
    delete best_location;


  delete [] ranges_forward;
  delete [] ranges_reverse;
  //   delete [] pos_occ;
  //   delete [] strand_occ;
  delete [] nb_locs;
  delete [] support;
}


float Support::getAverage() {
  return average_support;
}

float Support::getVariance() {
  float variance=0;
  for (uint i=0 ; i<getLength() ; i++){
    variance += pow((getSupport(i)-getAverage()),2);
  }
  variance /= getLength();
  return variance;
}

float Support::getStandardDeviation() {
  if (getVariance())
    return sqrt(getVariance());
  else 
    return 0;
}

float Support::getCoefficientVariation() {
  if (getAverage())
    return getStandardDeviation()*1.0 / getAverage();
  return 0;
}

SupportBreak *Support::getBreak(uint i, bool consider_strand) {
  if (! consider_strand || ! getBestLocation() 
      || getBestLocation()->getStrand() == 1)
    return breaks[i];
  else
    return breaks[getNbBreaks() - i - 1];
}

int Support::getEndPosRepeat() {
  return end_pos_repeat;
}

ReadIndex *Support::getIndexTags(){
  return tags;
}

LocateOnGenome *Support::getGenome() {
  return genome;
}

uint Support::getLength() {
  return length;
}


uint Support::getNbBreaks() {
  return nb_breaks;
}

uint Support::getNbDuplicate() {
  return nb_duplicate;
}

uint *Support::getNbLocs() {
  return nb_locs;
}

uint Support::getNbLocs(uint i) {
  return nb_locs[i];
}

uint Support::getNbLocsMax() {
  return nb_locs_max;
}

uint Support::getNbMultiple() {
  return nb_multiple;
}

uint Support::getNbPositionsLocated() {
  return nb_pos_located;
}

uint Support::getNbSingle() {
  return nb_single;
}

Parameters *Support::getParameters() {
  return parameters;
}

ChrPosition *Support::getBestLocation(){
  return best_location;
}

uint Support::getPositionOfLocation() {
  return position_of_location;
}

bool Support::hasPositionOfLocation() {
  return (position_of_location != read->getLength());
}

ChrPosition *Support::getLocation(uint i) {
  return locations[i];
}

ChrPosition **Support::getLocations() {
  return locations;
}

pair<uint, uint> Support::getKmerRange(uint i, uint strand) {
  if (i<length && i>=0) {
    if (strand == FORWARD_STRAND)
      return ranges_forward[i];
    else
      return ranges_reverse[i];
  } else {
    throw ILLEGAL_STATE_EXCEPTION;
  }
}

int Support::getRepeatLength() {
  if (getEndPosRepeat() == -1 && getStartPosRepeat() == -1)
    return 0;
  else if (getStartPosRepeat() == -1)
    return getEndPosRepeat() + 1;
  else if (getEndPosRepeat() == -1)
    return (int)(getLength() - getStartPosRepeat() + 1);
  else
    return (int)(getEndPosRepeat() -getStartPosRepeat() + 1);
}


int Support::getStartPosRepeat() {
  return start_pos_repeat;
}

uint *Support::getSupport() {
  return support;
}

Support *Support::getPairSupport(){
  return pair_support;
}

uint Support::getSupport(uint i) {
  return support[i];
}

Read *Support::getRead(){
  return read;
}

uint Support::getThreshold(){
  return tags->getFactorLength();
}

bool Support::hasRepeat() {
  return (getStartPosRepeat() != -1 || getEndPosRepeat() != -1) 
    && getRepeatLength() >= parameters->percent_min_unique_repetition*getNbPositionsLocated();
}

bool Support::isAlmostNormal(){
  return almostNormal;
}

bool Support::isContinuous() {
  return getLength() == (getNbMultiple() + getNbSingle() + getNbDuplicate());
}

bool Support::isDuplicate() {
  return !isNone() 
    && getNbDuplicate() > 0
    && getNbDuplicate() >= (uint) (parameters->percent_min_duplicate 
                                   * getNbPositionsLocated());
} 

bool Support::isMultiple() {
  return ! isNone() 
    && getNbMultiple() > 0
    && getNbMultiple() >= (uint) (parameters->percent_min_multiple 
				  * getNbPositionsLocated());
}
  
bool Support::isNone() {
  bool is_none = ((getNbSingle() == 0 && getNbMultiple() == 0 && getNbDuplicate() == 0)
		  || getNbLocsMax() == 0);
  if (nb_breaks == 1 && !is_none)
    is_none = (breaks[0]->isAtStart()
	       && breaks[0]->isAtEnd()
               );
  return is_none;
}

bool Support::isSingle() {
  return !isNone() 
    && getNbSingle() > 0
    && getNbSingle() >= (uint) (parameters->percent_min_single
                                * getNbPositionsLocated());
}

void Support::setBestLocation(ChrPosition *best_candidate){
  if (best_location)
    delete best_location;
  best_location = new ChrPosition(*best_candidate);
}

// PRIVATE
void Support::setLocationsAt(uint pos_loc){
  ulong *occs_fwd=NULL;
  uint nb_occs_fwd=0;
  ulong *occs_rev=NULL;
  uint nb_occs_rev=0;
  // multiple is initialized to nb_max_diplayed (--max-locs)
  ulong nb_locs_multiple = getGenome()->getNbLocations(); 

  // First, on forward strand
  if (ranges_forward[pos_loc].first <= ranges_forward[pos_loc].second) {
    genome->getOccurrencesFMIndexRange(ranges_forward[pos_loc].first,
				       ranges_forward[pos_loc].second,
				       &occs_fwd);
    nb_occs_fwd = ranges_forward[pos_loc].second - ranges_forward[pos_loc].first + 1;
  }
  // We check that the number of occurrences is smaller than nb_locs_multiple before the reverse strand treatment
  if (nb_occs_fwd < nb_locs_multiple) {
    if (ranges_reverse[pos_loc].first <= ranges_reverse[pos_loc].second) {
      genome->getOccurrencesFMIndexRange(ranges_reverse[pos_loc].first,
					 ranges_reverse[pos_loc].second,
					 &occs_rev);
      nb_occs_rev = ranges_reverse[pos_loc].second - ranges_reverse[pos_loc].first + 1;
    }
  }
  // Again, we check that the total number of occurrences is smaller than nb_locs_multiple
  // and we adjust the nb_locs_multiple if the limit is not reached
  if ((nb_occs_fwd+nb_occs_rev) < nb_locs_multiple)
    nb_locs_multiple=nb_occs_rev+nb_occs_fwd;
  
  // We save locations
  if (nb_locs_multiple>0){
    if (locations){
      for (uint i=0 ; i<getNbLocsMax() ; i++){
	delete locations[i];
      }
      free(locations);
    }
    locations = (ChrPosition **)malloc(sizeof(ChrPosition *) * nb_locs_multiple);
    uint nb_displayed = min(nb_occs_fwd,nb_locs_multiple);
    for (uint i=0 ; i<nb_displayed ; i++){
      locations[i] = getGenome()->getChrPos(occs_fwd[i],1);
    }
    for (uint i=nb_displayed ; i<nb_locs_multiple ; i++){
      locations[i] = getGenome()->getChrPos(occs_rev[i-nb_displayed],-1);
    }
  }
  
  // We delete tempory tables
  if (nb_occs_fwd)
    free(occs_fwd);
  if (nb_occs_rev)
    free(occs_rev);

  // set nb locs max
  nb_locs_max = nb_locs_multiple;
}

void Support::computeBestLocation(bool only_single, bool accept_single) {
  // Look for a not isolated 1 in nb_locs
  uint longest_run = 0, current_run = 0;
  uint start_longest_run = 0, start_current_run = 0;
  bool in_run = false;
  bool drop = false;
  bool find_start_break = false;
  uint start_last_break = ~0;
  uint offset_break = 1;
  uint nb_chuncks = 0;
  
  for (uint i=0; i < getLength(); i++) {
    // We prefer single_loc in single mode but a multipe_loc still good
    if (nb_locs[i] > 0
	&& support[i] > 0
	&& ((only_single && nb_locs[i] == 1) || !only_single)){
      
      // We could not consider a chunck in a drop
      if (!drop){
	// if !in_run, a new chunck is consider (probably the first one)
	if (!in_run){
	  in_run = true;
	  start_current_run = i;
	  current_run = 1;
	}
	// The chunck is increased outside a drop (the process seems good)
	else{
	  current_run++;
	}
      }

      // cerr << "in_run: " << in_run << " ,in_drop: " << drop <<" ,i: " << i << " ,start_current: " << start_current_run << " ,current_run: " << current_run << endl;
      
      // Single_loc inside a drop is not good (it could be a sequence error)
      // We consider a chunck before the drop
      if (in_run){
	if (!drop && i>0){
	  if (isFallenSupport(support[i-1],support[i])){
	    // cerr << "Fall and drop" << endl;
	    if ((!only_single || isSingleConsistent(start_current_run,current_run-1,parameters->min_bases_before_break)) 
		&& (!find_start_break || isGoodBreak(start_current_run, current_run-1, start_last_break, offset_break))
		&& !isOscillateLoc(start_current_run,current_run-1)){
	      if ((current_run-1) > longest_run){
		longest_run = current_run-1;
		start_longest_run = start_current_run;
		// cerr << "Good chunck1: " << start_current_run << ", " << longest_run << endl;
	      }
	      nb_chuncks++;
	    }
	    drop = true;
	    in_run = false;
	    // the drop length due to a seq error is at least k
	    i += (getThreshold()-1);   
	  }
	}
	if (i < (getLength()-1)){
	  if (isFallenSupport(support[i+1],support[i])){
	    // cerr << "Fall and no drop" << endl;
	    // we consider a new chunck after the drop and we do not consider the current chunck
	    if (drop
		|| (find_start_break && !isGoodBreak(start_current_run, current_run, start_last_break, offset_break))
		|| !isSingleConsistent(start_current_run,current_run,parameters->min_bases_before_break)){
	      in_run = false;
	      // cerr << "Not a good chunck2: " << start_current_run << ", " << current_run << endl;
	    }
	    drop = false;
	  }
	}
      }
    }
    // bad condition, we stop and save the chunck if it is good
    else {
      if (in_run) {
	in_run = false;
	if ((!only_single || isSingleConsistent(start_current_run,current_run,parameters->min_bases_before_break)) 
	    && !drop
	    && (!find_start_break || isGoodBreak(start_current_run, current_run, start_last_break, offset_break))
	    && !isOscillateLoc(start_current_run,current_run)){
	  // cerr << "Good chunck3: " << start_current_run << ", " << current_run << endl;
	  if (current_run > longest_run){
	    longest_run = current_run;
	    start_longest_run = start_current_run;
	  }
	  nb_chuncks++;
	}
      }
      // save the last start break_length
      if (nb_locs[i]==0){
	if (!find_start_break){
	  find_start_break = true;
	  start_last_break = i;
	  offset_break = 1;
	  // cerr << "First break: " << i << endl;
	}else if (start_last_break+offset_break == i){
	  offset_break++;
	}else {
	  // cerr << "other break: " << i << " ,start_last: " << start_last_break << " ,offset: " << offset_break << endl;
	  start_last_break = i;
	  offset_break = 1;
	}
      }
    }
  }
  
  // we check the last chunck
  if (in_run) {
    if ((!only_single || isSingleConsistent(start_current_run,current_run,parameters->min_bases_before_break)) 
	&& !drop
	&& (!find_start_break || isGoodBreak(start_current_run, current_run, start_last_break, offset_break))
	&& !isOscillateLoc(start_current_run,current_run)){
      // cerr << "Good chunck4: " << start_current_run << ", " << current_run << endl;
      if (current_run > longest_run){
	longest_run = current_run;
	start_longest_run = start_current_run;
      }
      nb_chuncks++;
    }
  }

  // Best position is i/ single_loc case: the position in the middle of the longest chunck
  //                 ii/ multiple_loc case: the position of the median locs 
  if (nb_chuncks*longest_run >= parameters->percent_min_chunck*read->getLength()){
    if (only_single){
      position_of_location = start_longest_run + (longest_run - 1)/2;
    }else{
      position_of_location = getMedianRun(start_longest_run,longest_run,accept_single);
    }
  }
  else if (!only_single){
    position_of_location = getMedianRun(0,getLength(),accept_single);
  }
  // cerr << "single? "<< only_single << ", start_run: " << start_longest_run << " ,longest: " << longest_run << " ,position of loc: " << position_of_location << " ,nb_loc: " << nb_locs[position_of_location] << endl;
}

bool Support::isFallenSupport(uint s1, uint s2){
  float ratio_drop = 1;
  if (s1 > 0)
    ratio_drop= s2*1.0/s1;
   
  return ((s2 == 1 && s1 > 2) || (s2 > 1 && s1 > 3 && ratio_drop <= parameters->percent_support_variation_almost_normal));
}

bool Support::isSingleConsistent(uint start, uint length, uint windows){
  // checking for the chunck length 
  bool consistent = true; //(length >= windows || getLength()/2 <= windows) ;    
  uint end = start+length;
  for (uint i = 0 ; i < windows && consistent ; i++){
    // if (start > 0)
    //   cerr << "start: " << start-1 << " ,loc_start: " << nb_locs[start-1] << endl;
    // if (end < getLength())
    //   cerr << "end: " << end << " ,loc_end: " << nb_locs[end] << endl;
    
    // checking for the locations around the chunck single 
    consistent = ((start == 0 || nb_locs[start-1] <= 1 || (!isFallenSupport(support[start],support[start-1]) && !isFallenSupport(support[start-1],support[start])))
		  && (end == getLength() || nb_locs[end] <= 1 || (!isFallenSupport(support[end],support[end-1]) && !isFallenSupport(support[end-1],support[end]))));
    
    if (start>0)
      start--;
    if (end<getLength())
      end++;
  }
  return consistent;
}

bool Support::isPairedEndOrientationConsistent(ChrPosition *pos1, ChrPosition *pos2) {
  switch (parameters->paired_end_orientation) {
    case FORWARD_REVERSE:
	    return (pos2->getStrand() == -1)? pos1->getRelativePosition() <= (pos2->getRelativePosition() + (long) getPairSupport()->getRead()->getLength()) : pos2->getRelativePosition() < (pos1->getRelativePosition() + (long) getPairSupport()->getRead()->getLength());
    case REVERSE_FORWARD:
	    return (pos2->getStrand() == 1)? pos1->getRelativePosition() <= (pos2->getRelativePosition() + (long) getPairSupport()->getRead()->getLength()) : pos2->getRelativePosition() < (pos1->getRelativePosition() + (long) getPairSupport()->getRead()->getLength());
    case FORWARD_FORWARD:
	    return pos1->getRelativePosition() <= (pos2->getRelativePosition() + (long) getPairSupport()->getRead()->getLength());
    default:
      cerr << "Wrong paired_end_orientation value" << endl;
      exit(1);
  }
}

bool Support::isOscillateLoc(uint start, uint length){
  if (length == 0)
    return false;
    
  uint min_loc = ~0;
  uint max_loc = 0;
  float average = 0;
  float variance = 0;
  float standard_deviation = 0;
  float standard_error = 0;
  for (uint i=0 ; i<length ; i++){
    average += nb_locs[start+i];
    if (min_loc > nb_locs[start+i])
      min_loc = nb_locs[start+i];
    if (max_loc < nb_locs[start+i])
      max_loc = nb_locs[start+i];
  }
  average /= length;

  for (uint i=0 ; i<length ; i++){
    variance += pow((nb_locs[start+i]-average),2);
  }
  variance /= length;
  
  standard_deviation = sqrt(variance);

  standard_error = standard_deviation/sqrt(length);
  
   // cerr << "average: " << average << " ,variance: " << variance << " ,standard_deviation: " << standard_deviation << " ,standard_error: " << standard_error << " ,min: " << min_loc << " ,max: " << max_loc << endl;
    if (min_loc <= parameters->max_localisation_duplication && max_loc > min_loc*100)
    return true;
  else
    return (max_loc > parameters->max_localisation_duplication
	    && (average * parameters->percent_max_loc_variation < standard_error));
}

bool Support::isGoodBreak(uint start_chunck, uint chunck_length, uint start_break, uint break_length){
  // cerr << "start_break: " << start_break <<" ,break_length: " << break_length << " ,start_chunck: " << start_chunck << " ,chunck_length: " << chunck_length << endl;
  
  // if (break_length < parameters->min_break_length)
  //   return false;
  if (start_chunck > start_break && (start_chunck-start_break+chunck_length) < (getThreshold()-1))
    return false;
  else if (start_chunck < start_break && (start_break-start_chunck+break_length) < (getThreshold()-1))
    return false;
  return true;
}


uint Support::getMedianRun(uint start, uint length, bool accept_single){
  vector <pair<uint,uint> > sort_locs;
  uint total = 0;
  for (uint i=0 ; i<length ; i++){
    if (nb_locs[start+i] > 1 || (nb_locs[start+i] > 0 && accept_single && support[start+i]>1)) {
      pair <uint,uint> ploc;
      ploc.first = nb_locs[start+i];
      ploc.second = start+i;
      sort_locs.push_back(ploc);
      total++;
    }
  }
  if (total > 0){
    sort(sort_locs.begin(), sort_locs.end(), pairCompare);
    uint i = (total-1)/2;
    pair <uint,uint> ploc = sort_locs[i];
    return ploc.second;
  }
  else if (length > 0 && !accept_single){
    return 1;
  }
  return read->getLength();
}

void Support::tryToMergeBreaks() {
  BreakList break_list(breaks, nb_breaks);

  // Is the first break located?
  // Should we do a loop?
  if (break_list[0]->hasNoStartBreak() && ! break_list[1]->isNiceBreak()
      && break_list.getGapSize(0) < parameters->min_bases_before_break) {
    break_list.merge(0);
  }

  size_t last_index = break_list.size() - 1;
  // Is the last break located?
  // Should we do a loop?
  if (last_index > 0 && break_list[last_index]->hasNoEndBreak() 
      && ! break_list[last_index-1]->isNiceBreak()
      && break_list.getGapSize(last_index-1) < parameters->min_bases_before_break) {
    break_list.merge(last_index-1);
  }

  // Merging short chimeric breaks

  size_t merged = ~0;
  for (size_t i = 0; i < break_list.size(); i++) {
    if (break_list[i]->isShortChimeric()) {
      merged = break_list
        .conditional_neighbor_merge(i,
                                     [=] (size_t i, size_t j) {
                                      return break_list.hasShortBreakSpan(i,j)
                                       && break_list[i]->isNiceMerge(break_list[j]);
                                     });
      if (merged == (size_t) ~0) {
        // size_t current_break = i;
        merged = break_list
          .conditional_neighbor_merge(i,
                                      [=] (size_t i, size_t j) {
                                        return break_list.hasShortBreakSpan(i,j)
                                        || break_list.hasShortGapBetweenBreaks(i, j);
                                      });

        if (merged == (size_t) ~0) {
          // Merge no (start/end) break with the neighbour if they are close
          merged = break_list
            .conditional_neighbor_merge(i,
                                        [=] (size_t i, size_t j) {
                                          return (break_list[i]->hasNoStartBreak()
                                                  || break_list[j]->hasNoEndBreak())
                                          && break_list.getGapSize(i, j)
                                          <= parameters->max_bases_randomly_matched;});
        }
      }

      if (merged != (size_t) ~0) {
        i = merged;
      }
    }
  }
    

  // Merging chimeric_breaks
  for (size_t i = 0; i < break_list.size(); i++) {
    if (break_list[i]->isChimeric()) {
      // Merge if we obtain a nice merge and the total size is not too large or the gap size is small
      merged = break_list
        .conditional_neighbor_merge(i,
                                    [=] (size_t i, size_t j) {
                                      return break_list[i]->isNiceMerge(break_list[j])
                                      && (break_list.hasShortBreakSpan(i,j) 
                                          || break_list.hasShortGapBetweenBreaks(i,j));});
      
      // If we have a chimeric break next to us, we merge it if it is not too
      // far.
      if (merged == (size_t)~0) {
        if (i < break_list.size() - 1 && break_list[i+1]->isChimeric()
            && (break_list.hasShortBreakSpan(i, i+1) 
                || break_list.hasShortGapBetweenBreaks(i ,i+1)))
        break_list.merge(i);
      
        else {
          merged = break_list
            .conditional_neighbor_merge(i,
                                        [=] (size_t i, size_t j) {
                                          return (break_list[i]->hasNoStartBreak()
                                                  || break_list[j]->hasNoEndBreak())
                                          && break_list.getGapSize(i,j) 
                                          <= parameters->max_bases_randomly_matched;});
        }
      }
      if (merged != (size_t)~0)
        i = merged;
    }
  }

  // Merging other breaks
  for (size_t i = 0; i < break_list.size() - 1; i++) {
    if (! break_list[i]->hasNoStartBreak() && ! break_list[i]->hasNoEndBreak()
        && break_list.getTotalBreakSize(i) <= getThreshold()
        && (! break_list[i]->isVeryNiceBreak() || ! break_list[i+1]->isVeryNiceBreak())
        && break_list[i]->isVeryNiceMerge(break_list[i+1])) {
      break_list.merge(i);
    } else if (break_list[i]->getPositionEndBreak() 
               >= break_list[i+1]->getPositionStartBreak()) {
      // If the breaks overlap, we merge them.
      break_list.merge(i);
    }
  }

  uint end_of_previous_break = 0;
  bool relaunch_merge = false;
  // Just keep the necessary breaks
  for (size_t i = 0; i < break_list.size(); i++) {
    breaks[i] = break_list[i];

    // adjust start_break and end_break positions of splicing events (splices, chimeras) 
    // for overlaping cases (break_length < k-1)
    // we use GT--AG splicing sites
    breaks[i]->adjustBreakWithSplicingSites();
    if (breaks[i]->getPositionStartBreak() < end_of_previous_break)
      relaunch_merge = true;
    end_of_previous_break = breaks[i]->getPositionEndBreak();
  }
    
  if (break_list.size() < nb_breaks) {
    // Realloc memory so that we don't take unnecessary memory
    breaks = (SupportBreak **)realloc(breaks, 
                                      break_list.size() * sizeof(SupportBreak *));
    nb_breaks = break_list.size();
  }

  if (relaunch_merge)
    tryToMergeBreaks();
}

void Support::checkPairedConnection(){
  int opposite_strand = parameters->paired_end_orientation == FORWARD_FORWARD? 1 : -1;
  Bitset *false_positives = new Bitset(getNbLocsMax());
  Bitset *false_positives_pair = new Bitset(getPairSupport()->getNbLocsMax());
  Bitset *false_positives_better = new Bitset(getNbLocsMax());
  Bitset *false_positives_pair_better = new Bitset(getPairSupport()->getNbLocsMax());
  uint min_distance = ~0;
  bool find_candidat = false;
  bool find_candidat_better = false;
  uint min_distance_better = ~0;
  //int best_strand1, best_strand2;
  //uint best_pos1, best_pos2;
  //char *best_chr1;
  //char *best_chr2;
  uint nb_candidats=0, nb_candidats_better=0;
  uint best_candidate1, best_candidate2;
   
  for (uint i=0 ; i < getNbLocsMax() ; i++){
    for (uint j=0 ; j < getPairSupport()->getNbLocsMax() ; j++){
      int strand2 = getLocation(i)->getStrand();
      char *chr2 = getLocation(i)->getChrPosition();
      long pos2 = (long) getLocation(i)->getRelativePosition();
      int strand1 = getPairSupport()->getLocation(j)->getStrand()*opposite_strand;
      char *chr1 = getPairSupport()->getLocation(j)->getChrPosition();
      long pos1 = (long) getPairSupport()->getLocation(j)->getRelativePosition();
      if (strcmp(chr1,chr2)==0 && strand1 == strand2){
	// cerr << endl << "[ chr1, pos1, strand1 ]:" << "[ " << chr1 <<", " << pos1 <<", " << strand1 << " ]" << endl;
	// cerr << "[ chr2, pos2, strand2 ]:" << "[ " << chr2 <<", " << pos2 <<", " << strand2 << " ]" << endl;    
	// cerr << "nb_candidats: " << nb_candidats << endl;
	if (labs(pos2-pos1) < min_distance && !find_candidat_better){
	  min_distance = labs(pos2-pos1);
	  //best_pos2 = pos2;
	  //best_chr2 = chr2;
	  //best_strand2 = strand2;
	  //best_pos1 = pos1;
	  //best_chr1 = chr1;
	  //best_strand1 = strand1;
    best_candidate1 = j;
    best_candidate2 = i;
	}
	false_positives->setBit(i);
	false_positives_pair->setBit(j);
	find_candidat = true;
	nb_candidats++;
	if (isPairedEndOrientationConsistent(getPairSupport()->getLocation(j),getLocation(i))
	    // && labs(pos1-pos2) <= parameters->max_splice_length
	    ){
	  find_candidat_better = true;
	  nb_candidats_better++;
	  // cerr << "nb_candidats_better: " << nb_candidats_better << endl;
	  if (labs(pos2-pos1) < min_distance_better){
	    min_distance_better = labs(pos2-pos1);
	    //best_pos2 = pos2;
	    //best_chr2 = chr2;
	    //best_strand2 = strand2;
	    //best_pos1 = pos1;
	    //best_chr1 = chr1;
	    //best_strand1 = strand1;
      best_candidate1 = j;
      best_candidate2 = i;
	  }
	  false_positives_better->setBit(i);
	  false_positives_pair_better->setBit(j);
	}
      }
    }
  }

  if (find_candidat){
    setBestLocation(getLocation(best_candidate2));
    getPairSupport()->setBestLocation(getPairSupport()->getLocation(best_candidate1));
    float canditats_ratio = nb_candidats_better*1.0/nb_candidats;
    if (find_candidat_better && (nb_candidats < parameters->max_good_paired_reads || canditats_ratio > parameters->percent_very_good_paired)){
      removeLocations(false_positives_better);
      getPairSupport()->removeLocations(false_positives_pair_better);
    }else{
      removeLocations(false_positives);
      getPairSupport()->removeLocations(false_positives_pair);
    }
  }
  delete false_positives;
  delete false_positives_pair;
  delete false_positives_better;
  delete false_positives_pair_better;
}

void Support::removeLocations(Bitset *false_positives){
  uint j=0;
  uint nb_false_positives=0;
  for (uint i=0 ; i<getNbLocsMax() ; i++){
    ChrPosition *tmp = new ChrPosition(*locations[i]);
    delete locations[i];
    if (false_positives->isBitSet(i)){
      locations[j] = tmp;
      j++;
    }else{
      delete tmp;
      nb_false_positives++;
    }
  }
  nb_locs_max -= nb_false_positives;
}
