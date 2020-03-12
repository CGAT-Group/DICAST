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

#include <string>
#include <limits.h>
#include "SupportBreak.h"
#include "Support.h"
#include "utils.h"
#include <math.h>

using namespace std;

SupportBreak::SupportBreak(uint start, uint end, 
			   Support *s, Parameters *p,
			   pair<uint, uint> *range_occs_forward,
			   pair<uint, uint> *range_occs_reverse,
         uint nb_merges):
  k(s->getThreshold()),
  pos_start_break(start),pos_end_break(end),
  support(s),parameters(p),
  range_occs_fwd(range_occs_forward), range_occs_rev(range_occs_reverse),
  nb_merges(nb_merges){
  init();
}

SupportBreak::SupportBreak(SupportBreak &b1, SupportBreak &b2) {
  SupportBreak *first, *last;
  if (b1.support != b2.support) {
    throw ILLEGAL_STATE_EXCEPTION;
  }
  if (b2.getPositionEndBreak() > b1.getPositionEndBreak()) {
    first = &b1;
    last = &b2;
  } else {
    first = &b2;
    last = &b1;
  }
  uint pos_break = first->getPositionStartBreak();
  uint pos_end = last->getPositionEndBreak();

  // TODO: why check the two conditions? Is it not the same?
  if (first->hasNoStartBreak() && pos_break == 0) {
    // If the k-mers are located, we have extended the breaks.  We come back
    // to the original positions so that the extension can be performed
    // (again) after the merge without problem
    while (pos_break < pos_end && first->support->getNbLocs()[pos_break] > 0)
      pos_break++;    
  }
  if (last->hasNoEndBreak() && pos_end == last->getLength() - 1) {
    while (pos_end > pos_break && last->support->getNbLocs()[pos_end] > 0)
      pos_end--;
  }

  k = b1.getThreshold();
  pos_start_break = pos_break;
  pos_end_break = pos_end;
  support = first->support;
  parameters = first->parameters;
  range_occs_fwd = first->range_occs_fwd;
  range_occs_rev = first->range_occs_rev;
  nb_merges = b1.getNbMerges() + b2.getNbMerges() + 1;
  init();
}

SupportBreak::~SupportBreak(){
  for (int i=0; i < nb_candidats; i++) {
    delete candidats[i];
  } 
  free(candidats);
  if (candidat_chosen != NULL)
    delete candidat_chosen;
}

void SupportBreak::computeScore(){
  computeScore(pos_start_break, pos_end_break);
}

void SupportBreak::computeScore(uint pos_start, uint pos_end){
  // we can not check a location after getLength()-1 in Support
  uint max_pos_end = min(pos_end, support->getLength()-1);
  uint max_pos_start = min(pos_start, support->getLength()-1);
  // Computing score outside break
  setInsideScore(max_pos_start, max_pos_end);

  // Computing score outside break
  setOutsideScore(max_pos_start, max_pos_end);
}

float SupportBreak::getMedianInside() {
  return median_inside;
}

float SupportBreak::getVarianceInside() {
  float variance_inside=0;
  uint nb_inside=0;
  for (uint i=pos_start_break ; i<=pos_end_break ; i++){
    if (support->getNbLocs()[i] == 0) {
      variance_inside += pow((support->getSupport(i)-getInsideScore()),2);
      nb_inside++;
    }
  }
  variance_inside /= nb_inside;
  return variance_inside;
}

float SupportBreak::getStandardDeviationInside() {
  if (getVarianceInside())
    return sqrt(getVarianceInside());
  else 
    return 0;
}

float SupportBreak::getCoefficientVariationInside() {
  if (getInsideScore())
    return getStandardDeviationInside()*1.0 / getInsideScore();
  return 0;
}

float SupportBreak::getAverageLowInside() {
  return average_low_inside;
}

float SupportBreak::getAverageHighInside() {
  return average_high_inside;
}

CandidatBreak *SupportBreak::getCandidatChosen(){
  return candidat_chosen;
}

int SupportBreak::getChimeraClass(){
  if (getTagBreakLength() < getThreshold() && hasLongEnoughTagBreak()){
    if  (getChrId(START_BREAK) != getChrId(END_BREAK))
      return 1;
    if (getStrandLocation(START_BREAK) != getStrandLocation(END_BREAK))
      return 4;
    if (getGenomeGapLength() < 0)
      return 3;
    if (labs(getGenomeGapLength()) > parameters->max_splice_length)
      return 2;
  }
  return 0;
}

const uchar *SupportBreak::getChr(pos_location_break p){
  if (((p == START_BREAK && hasNoStartBreak()) 
       || (p == END_BREAK && hasNoEndBreak()))
      ) 
    return NULL;
  if (p == START_BREAK)
    return getGenome()->getChromosome(loc_start_break);
  else
    return getGenome()->getChromosome(loc_end_break);
}


ulong SupportBreak::getChrId(pos_location_break p){
  if (((p == START_BREAK && hasNoStartBreak())
       || (p == END_BREAK && hasNoEndBreak()))
      )
    return ~0;
  if (p == START_BREAK)
    return getGenome()->getIdChromosome(loc_start_break);
  else
    return getGenome()->getIdChromosome(loc_end_break);
}
 
ChrPosition *SupportBreak::getChrPosition(pos_location_break p){
  if (((p == START_BREAK && hasNoStartBreak()) 
      || (p == END_BREAK && hasNoEndBreak()))
      )
    return NULL;
  if (p == START_BREAK) {
    return getGenome()->getChrPos(loc_start_break,
					   strand_start_break);
  }
  else 
    return getGenome()->getChrPos(loc_end_break,strand_end_break);
}


LocateOnGenome *SupportBreak::getGenome(){
  return support->getGenome();
}

gap_size_t SupportBreak::getGenomeGapLength(){
  if (hasNoStartBreak() || hasNoEndBreak() 
      )
    return GAP_SIZE_MAX;
  if (strand_start_break != strand_end_break
      || (getChrId(START_BREAK) != getChrId(END_BREAK)))
    return min(labs((long long int) loc_end_break - loc_start_break),
	       labs((long long int) loc_start_break - loc_end_break)) - 1;
  if (strand_end_break == -1)
    return (gap_size_t)loc_start_break - (gap_size_t)loc_end_break-1;
  return (gap_size_t)loc_end_break - (gap_size_t)loc_start_break-1;
}

uint SupportBreak::getInsideQuartile1() {
  return inside_quartile1;
}

uint SupportBreak::getInsideQuartile4() {
  return inside_quartile4;
}

float SupportBreak::getInsideScore(){
  return score_during_break;
}

uint SupportBreak::getLength(){
  return getTagLength() - getThreshold() + 1;
}

ChrPosition *SupportBreak::getLocationEndBreak() {
  return getChrPosition(END_BREAK);
}

ChrPosition *SupportBreak::getLocationStartBreak() {
  return getChrPosition(START_BREAK);
}

int SupportBreak::getNbCandidats(){
  return nb_candidats;
}

uint SupportBreak::getNbMerges(){
  return nb_merges;
}

int SupportBreak::getNbGenomeIndels(){
  if (this->getGenomeGapLength() == GAP_SIZE_MAX)
    return INT_MAX;
  return (int)(this->getGenomeGapLength() 
	       - (int)this->getTagBreakLength());
}

int SupportBreak::getNbTagIndels(){
  //return (int)((int)this->getTagBreakLength() 
  //	       - (int)getThreshold());
  return (int)((int)this->getTagBreakLength() 
	       -this->getGenomeGapLength());
}

float SupportBreak::getOutsideScore(){
  return score_outside_break;
}

Parameters *SupportBreak::getParameters() {
  return parameters;
}

uint SupportBreak::getPositionEndBreak(int strand){
  if (strand == 1)
    return pos_end_break;
  return getLength() - 1 - getPositionStartBreak();
}

uint SupportBreak::getPositionEndBreakOnTheRead() {
  if (! hasPositionEndBreakOnTheRead()) {
    return getPositionEndBreak() + getThreshold() - 1;
  }
  return override_pos_end_break;
}

ulong SupportBreak::getPositionLocation(pos_location_break p){
  if (((p == START_BREAK && hasNoStartBreak()) 
    || (p == END_BREAK && hasNoEndBreak()))
      // || isRepeated()
      )
    return ~0;
  if (p == START_BREAK)
    return loc_start_break;
  else
    return loc_end_break;
}

uint SupportBreak::getPositionStartBreak(int strand){
  if (strand == 1)
    return pos_start_break;
  return getLength() - 1 - getPositionEndBreak();
}


CandidatBreak *SupportBreak::getPotentialCandidat(uint i){
  return candidats[i];
}

pair<uint, uint> SupportBreak::getRangeEndBreak(uint strand) {
  if (hasNoEndBreak())
    throw ILLEGAL_STATE_EXCEPTION;

  if (strand == FORWARD_STRAND)
    return range_occs_fwd[pos_end_break+1];
  else
    return range_occs_rev[pos_end_break+1];

}

pair<uint, uint> SupportBreak::getRangeStartBreak(uint strand) {
  if (hasNoStartBreak())
    throw ILLEGAL_STATE_EXCEPTION;

  if (strand == FORWARD_STRAND)
    return range_occs_fwd[pos_start_break-1];
  else
    return range_occs_rev[pos_start_break-1];  
}

float SupportBreak::getScoreInsideAverages() {
  return score_inside_averages;
}

float SupportBreak::getScoreComputedIntraExon() {
  // compute the distance between the curve and the point for the same outside
  return (getThresholdForScoreIntraExon(score_outside_break)-score_during_break);
}

float SupportBreak::getScoreComputedInterExon(){
  // compute the distance between the curve and the point for the same outside
  return (getThresholdForScoreInterExon(score_outside_break)-score_during_break);
}

float SupportBreak::getScoreInsideBreak(){
  return score_during_break;
}

float SupportBreak::getScoreOutsideBreak(){
  return score_outside_break;
}

int SupportBreak::getStrandLocation(pos_location_break p){
  if (((p == START_BREAK && hasNoStartBreak()) 
       || (p == END_BREAK && hasNoEndBreak()))
      // || isRepeated()
      )
    return 0;
  if (p == START_BREAK)
    return strand_start_break;
  else
    return strand_end_break;

}

uint SupportBreak::getTagBreakLength(){
  return getPositionEndBreak()-getPositionStartBreak() + 1;
}

uint SupportBreak::getOriginalTagBreakLength(){
  return original_break_length;
}

uint SupportBreak::getTagLength(){
  return support->getRead()->getLength();
}

uint SupportBreak::getThreshold(){
  return k;
}


float SupportBreak::getThresholdForScoreInterExon(float score){
//   return (1-(log(score_during_break + M_E - (1+parameters->min_support_no_cover))
// 	     /log(score_outside_break + M_E - (1+parameters->min_support_no_cover))));
  // float thresh = 981617 - 981610 * exp(-1.0*score_outside_break/60411843) 
//     - 5.5 * exp(-1.0*score_outside_break/22.9);
  //return  0.53793 + sqrt(0.24875 * score);
  // return  0.22453 + sqrt(0.20354 * score);
  return  0.51081 + sqrt(0.16758 * score); 
}

float SupportBreak::getThresholdForScoreIntraExon(float score){
//   float thresh = 0.35+0.5 * score_outside_break;
//   float thresh = -4.063*pow(M_E,-score_outside_break/6.22) 
//     - 3004 * pow(M_E,-score_outside_break/70240) + 3009;
  //return 0.8 * pow(score_outside_break,5./6) + 0.13387;
  //return  0.812 + 0.299 * pow(score,4./5);
  //return  0.73566 + 0.16216 * pow(score,4./5);
  //return  0.83566 + 0.16216 * pow(score,4./5);
  // return  0.93330 + 0.18609 * pow(score,4./5);  
  //return  -0.42376 + 1.17887 * pow(score,1./3) + 0.06945 * pow(score,5./6);  
  // return -0.33462 + 1.16893 * pow(score,1./3) + 0.04730 * pow(score,5./6);  
  // return -0.38608 + 0.85347 * pow(score,1./3) + 0.27864 * pow(score,5./6);   
  if (score <= 5)
    return 0.19227 - 0.82749 * score + 1.51619 * pow(score,5./6); 
  else
    return -2.40850 + 2.15859 * pow(score,1./3) + 0.15670 * pow(score,5./6);
}

// bool SupportBreak::hasLongEnoughTagBreak(){
//   return (getTagBreakLength() >= parameters->min_break_length);
// }

bool SupportBreak::hasLongEnoughTagBreak() {
  return  hasNoStartBreak() || hasNoEndBreak() 
    || getThreshold() < parameters->max_bases_randomly_matched
    || getTagBreakLength() > getThreshold() - parameters->max_bases_randomly_matched;
}

bool SupportBreak::hasNoCover() {
  return score_outside_break <=  parameters->min_support_no_cover
    && score_during_break <=  parameters->min_support_no_cover;
}

bool SupportBreak::hasNoEndBreak(int strand) {
  return getNbCandidats() == 0 || isAtEnd(strand);
}

bool SupportBreak::hasNoStartBreak(int strand) {
  return getNbCandidats() == 0 || isAtStart(strand);
}

bool SupportBreak::hasPositionEndBreakOnTheRead(int strand) {
  return override_pos_end_break != 0;
}

bool SupportBreak::isAtEnd(int strand) {
  return (strand == 1 ? getPositionEndBreakOnTheRead() + 1 == getTagLength() : isAtStart(1));
}

bool SupportBreak::isAtStart(int strand) {
  return (strand == 1 ? getPositionStartBreak() == 0 : isAtEnd(1));
}

bool SupportBreak::isBiologicalEvent(){
  return ! hasNoCover() 
    && (score_during_break >= score_outside_break
        || isBiologicalIntraEvent()
	|| isBiologicalInterEvent()
	);
}

bool SupportBreak::isBiologicalIntraEvent(){
  return getScoreComputedIntraExon() <= parameters->p_value_variation_biological
    && getChrId(START_BREAK) == getChrId(END_BREAK)
    && getStrandLocation(START_BREAK) == getStrandLocation(END_BREAK) 
    && getGenomeGapLength() > 0
    && getGenomeGapLength() <= (getTagBreakLength() + parameters->max_bio_ins_del); 
}

bool SupportBreak::isBiologicalInterEvent(){
  return getScoreComputedInterExon() <= parameters->p_value_variation_biological    
    && ! isBiologicalIntraEvent();
}
    

bool SupportBreak::isChimera(){
  return (getChimeraClass() != 0);
}

bool SupportBreak::isChimeric() {
  return ! hasNoStartBreak() && ! hasNoEndBreak() && ! isNiceBreak();
}

bool SupportBreak::isDeviated(){
  return
    getScoreInsideAverages() >= 0
    || (getAverageLowInside() < 1+EPSILON
        && getAverageHighInside() < parameters->max_ambiguous_average_high
        && (isBiologicalEvent() && isSupportFalling()));
}

bool SupportBreak::isDuplicated(){
  if (candidat_chosen == NULL)
    return false;
  else
    return candidat_chosen->isDuplicated();
}

bool SupportBreak::isRepeated() {
  return repeated;
}

bool SupportBreak::isGenomeDeletion(){
  if (getTagBreakLength() <= parameters->max_bio_ins_del)
    return false;
  return (
	  getChrId(START_BREAK) == getChrId(END_BREAK)
	  && getStrandLocation(START_BREAK) == getStrandLocation(END_BREAK) 
	  && getGenomeGapLength() > 0
	  && getGenomeGapLength() < getTagBreakLength() - parameters->max_bio_ins_del
	  );
}

bool SupportBreak::isGenomeInsertion(){
  if (getTagBreakLength() > getThreshold() - 1)
    return false;
  return (
	  getChrId(START_BREAK) == getChrId(END_BREAK)
	  && getStrandLocation(START_BREAK) == getStrandLocation(END_BREAK) 
	  && getGenomeGapLength() > getTagBreakLength() + parameters->max_bio_ins_del
	  && labs(getGenomeGapLength()) <= parameters->max_splice_length
	  );
}

bool SupportBreak::isNiceBreak() {
  if (!hasLongEnoughTagBreak())
    return false;
  return (isVeryNiceBreak() || isGenomeInsertion()
          || isGenomeDeletion())
    && ! isChimera();
  // return check_nice_merge(this, this, parameters->max_splice_length);
}

bool SupportBreak::isNiceMerge(SupportBreak *sb) {
  return check_nice_merge(this, sb, false);
}

bool SupportBreak::isShortChimeric() {
  return isChimeric() && ! hasLongEnoughTagBreak();
}

bool SupportBreak::isSplice(){
  return isGenomeInsertion() 
    && (labs(getGenomeGapLength()) > parameters->max_bio_ins_del);
}

bool SupportBreak::isSupportFalling() {
  return isSupportFallingRight() || isSupportFallingLeft();
}

bool SupportBreak::isSupportFallingLeft(int strand) {
  return (strand == 1
          ? (pc_ones_left_inside >= parameters->min_perc_ones_inside
             && pc_ones_left_inside 
             >= parameters->min_ratio_support_fall * pc_ones_left_outside)
          : isSupportFallingRight());
}

bool SupportBreak::isSupportFallingRight(int strand) {
  return (strand == 1
          ? (pc_ones_right_inside >= parameters->min_perc_ones_inside
             && pc_ones_right_inside 
             >= parameters->min_ratio_support_fall * pc_ones_right_outside)
          : isSupportFallingLeft());
}

bool SupportBreak::isTagSubstitution(){
  return (getChrId(START_BREAK) == getChrId(END_BREAK)
	  && getStrandLocation(START_BREAK) == getStrandLocation(END_BREAK)
	  && getGenomeGapLength() == getTagBreakLength()); 
}

bool SupportBreak::isTagIndel(){
  bool isAGoodTagBreakLength = true;
  if (getNbTagIndels() < 0){
    isAGoodTagBreakLength = (getTagBreakLength() < getThreshold());
  }else{
    isAGoodTagBreakLength = (getThreshold() + getNbTagIndels() > getGenomeGapLength());
  }
  return (isAGoodTagBreakLength
          && getChrId(START_BREAK) == getChrId(END_BREAK)
          && getStrandLocation(START_BREAK) == getStrandLocation(END_BREAK) 
          && ! isTagSubstitution() 
          && labs(getNbTagIndels()) <= parameters->max_bio_ins_del
          && getGenomeGapLength() > 0);
}

bool SupportBreak::isVeryNiceBreak() {
  return isTagIndel() || isTagSubstitution() 
    || (isGenomeInsertion() 
        && labs(getGenomeGapLength() <= parameters->max_verynice_merge));
}

bool SupportBreak::isVeryNiceMerge(SupportBreak *sb) {
  return check_nice_merge(this, sb, true);
}

void SupportBreak::setInsideScore(){
  setInsideScore(pos_start_break, pos_end_break);
}

void SupportBreak::setInsideScore(uint pos_start, uint pos_end){
  uint nb_inside = 0;
  score_during_break = 0;
  for (uint i = pos_start ; i<=pos_end ; i++){
    if (support->getNbLocs()[i] == 0) {
      score_during_break += support->getSupport(i);
      nb_inside++;
    }
  }
  score_during_break /= nb_inside; 
}

void SupportBreak::setLocationEndBreak(uint loc_end){
  loc_end_break = loc_end;
}

void SupportBreak::setLocationStartBreak(uint loc_start){
  loc_start_break = loc_start;
}

void SupportBreak::setOutsideScore() {
  setOutsideScore(pos_start_break, pos_end_break);
}

void SupportBreak::setOutsideScore(uint pos_start, uint pos_end) {
  uint pos_s, pos_e;
  if (pos_start < parameters->support_score_window_length/2)
    pos_s = 0;
  else
    pos_s = pos_start - parameters->support_score_window_length/2;
  if (support->getLength() 
      < 1 + pos_end + parameters->support_score_window_length/2)
    pos_e = support->getLength()-1;
  else
    pos_e = pos_end + parameters->support_score_window_length/2;

  float newScore=0;
  uint nb=0;
  for (uint i = pos_s; i < pos_start; i++) {
    if (support->getNbLocs()[i] == 0) {
      nb = 0;
      newScore = 0;
    } else {
      newScore += support->getSupport(i);
      nb++;
    }
  }
  for (uint i=pos_end+1; i <= pos_e 
	 && support->getNbLocs()[i] != 0; i++) {
    newScore += support->getSupport(i);
    nb++;
  }

  if (nb > 0)
    score_outside_break = newScore / nb;
  else 
    score_outside_break = (support->getNbLocs()[pos_start]
                           + support->getNbLocs()[pos_end]) / 2.;

  //compute some metrics on support variation inside break
  computeHighAndLowAverages(pos_start, pos_end);
}

void SupportBreak::setPositionEndBreak(uint pos_end){
  pos_end_break = pos_end;
}

void SupportBreak::setPositionEndBreakOnTheRead(uint pos_end){
  override_pos_end_break = pos_end;
}

void SupportBreak::setPositionStartBreak(uint pos_start){
  pos_start_break = pos_start;
}

void SupportBreak::setRangesForward(pair<uint, uint> *p) {
  range_occs_fwd = p;
}
  
void SupportBreak::setRangesReverse(pair<uint, uint> *p) {
  range_occs_rev = p;
}

void SupportBreak::setStrandEndBreak(int strand_end){
  strand_end_break = strand_end;
}
  
void SupportBreak::setStrandStartBreak(int strand_start){
  strand_start_break = strand_start;
}

void SupportBreak::setSupport(Support *s) {
  support = s;
}

void SupportBreak::setThreshold(uint k) {
  this->k = k;
}

// PRIVATE

void SupportBreak::init() {
  strand_start_break = 0;
  loc_start_break = ~0;
  strand_end_break = 0;
  loc_end_break = ~0;
  override_pos_end_break = 0;
  score_during_break = 0;
  candidat_chosen = NULL;
  repeated = false;
  
  original_break_length = getPositionEndBreak()-getPositionStartBreak() + 1;
  // maximal allocation to save candidats during extension process
  int max_nb_candidats = parameters->max_extension_length + 1;
  // we init nb_candidats == -1 in order to distinguish init case and no-candidat case
  nb_candidats = -1;
  candidats = (CandidatBreak**)malloc(sizeof(CandidatBreak*)*max_nb_candidats);

  if (pos_start_break > 0 || pos_end_break < getLength() - 1) 
    fillWithClosestMatch();
  
  if (candidat_chosen != NULL){ 
    pos_start_break = candidat_chosen->getPosStartBreak();
    loc_start_break = candidat_chosen->getLocStartBreak();
    strand_start_break = candidat_chosen->getStrandStartBreak();
    pos_end_break = candidat_chosen->getPosEndBreak();
    loc_end_break = candidat_chosen->getLocEndBreak();
    strand_end_break = candidat_chosen->getStrandEndBreak();
  }

  // compute the score inside and outside the break for the candidat_chosen
  computeScore();
}


void SupportBreak::adjustBreakWithSplicingSites(){
  if (hasNoStartBreak() || hasNoEndBreak())
    // We can't adjust if a part of the break is mising.
    return;
  // TODO: we only consider chimeric junction of class 2
  // to avoid a merge of overlap breaks...
  // But this behaviour could be improved
  if (getChimeraClass()==2 || isSplice()){
    int overlap = getThreshold() - getTagBreakLength() - 1;
    bool found_donorsite = false;
    char nuc1=0;
    char nuc2=0;
    bool out_of_scope = false;
    uint i=0;
    if (overlap > 0){
      // first step: we remove overlap by switching left
      // (only if we can...)
      if (strand_start_break == 1) {
        if(pos_start_break > (uint)overlap){
          pos_start_break -= overlap;
          loc_start_break -= overlap;
        } else {
          pos_start_break = 0;
          loc_start_break = ~0;
          strand_start_break = 0;
          out_of_scope = true;
        }
      }else {
        if(pos_end_break + overlap < support->getLength() - 1) {
          pos_end_break += overlap;
          loc_end_break -= overlap;
        } else {
          pos_end_break = support->getLength()-1;
          loc_end_break = ~0;
          strand_end_break = 0;
          out_of_scope = true;
        }
      }
      if(!out_of_scope) {
        // second step: adjust the windows according to the donor site
        // We loop until the donor site (if there is)
        // donor site: GT in fwd strand 
        //             AC in rev strand
        while (i < (uint)overlap && !found_donorsite && !out_of_scope){
          // Reset nucleotides value to avoid conflict with previous one
          nuc1 = 0;
          nuc2 = 0;
          if (strand_start_break == 1){
            if (pos_end_break+1+i < getLength()) {
              nuc1 = support->getRead()->seq[pos_end_break+1+i];
	      if (pos_end_break+2+i < getLength())
                nuc2 = support->getRead()->seq[pos_end_break+2+i];
	      if (nuc1 == 'G'){
                if (overlap == 1) 
                  found_donorsite = true;
                else if (nuc2 == 'T')
                  found_donorsite = true;
              }
            } else {
              out_of_scope = true;
            }
          }else{
            if (pos_start_break > i) {
              nuc1 = support->getRead()->seq[pos_end_break-i];
              if (pos_start_break > i + 1)
                nuc2 = support->getRead()->seq[pos_end_break-1-i];
              if (nuc1 == 'C'){
                if (overlap == 1) 
                  found_donorsite = true;
                else if (nuc2 == 'A')
                  found_donorsite = true;
              }
            } else {
              out_of_scope = true;
            }
          }
	  if (!found_donorsite)
            i++;
        }
        // we may have gone too far in the last i++ iteration
        // and we could have gone out of scope on the other side of the splice
        // (or the chimera);
        if(strand_start_break == 1 && pos_end_break + i > getLength() - 1) {
          out_of_scope = true;
        } else if(strand_start_break == -1 && pos_start_break <= i) {
          out_of_scope = true;
        }
        // If we are out of scope, we set a no start break or no_end_break according to the strand
        if(out_of_scope) {
          if(strand_start_break == 1) {
            pos_end_break = support->getLength()-1;
            loc_end_break = ~0;
            strand_end_break = 0;
          } else {
            pos_start_break = 0;
            loc_start_break = ~0;
            strand_start_break = 0;
          }
        }
        // we adjust start_break and end_break
        else if (i>0){
          pos_start_break += i*strand_start_break;
          pos_end_break += i*strand_start_break;
          loc_start_break += i;
          loc_end_break += i;
        }
      }
    }
  }
}

void SupportBreak::adjustBoundaries(){
  //  ulong nb_disp = getGenome()->getNbLocations();
  ulong *locs_start_after, *locs_end_before;
  ulong nb_elements_rev, nb_elements_fwd;
  uint pos_start_extension = candidat_chosen->getPosStartBreak();
  uint loc_start_extension =  candidat_chosen->getLocStartBreak();
  int strand_start_extension = candidat_chosen->getStrandStartBreak();
  uint pos_end_extension = candidat_chosen->getPosEndBreak();
  uint loc_end_extension =  candidat_chosen->getLocEndBreak();
  int strand_end_extension = candidat_chosen->getStrandEndBreak();
  
  bool hasNoStartExtension = false;
  bool hasNoEndExtension = false;
  bool findStart = false;
  bool findEnd = false;
  bool findElement;    

  // a bioUndetermined candidat is notify if the extension process has not to 
  // be finished and the candidat looks like a chimera (it could be a false positive)
  if (!candidat_chosen->isNiceCandidat() && (pos_start_extension == 0 || pos_end_extension == (getLength() - 1))){
    if (pos_end_extension == (getLength() - 1)){ 
      hasNoEndExtension = true;
    }
    if (pos_start_extension == 0){
      hasNoStartExtension = true;
    }
  }
    
  // Adjust the start break
  if (!hasNoStartExtension){
    if (strand_start_extension == -1){
      while (range_occs_rev[pos_start_extension+1].first <=
	     range_occs_rev[pos_start_extension+1].second && !findStart){
	getGenome()->
	  getOccurrencesFMIndexRange(range_occs_rev[pos_start_extension+1].first,
				     range_occs_rev[pos_start_extension+1].second,
				     &locs_start_after);
	
	nb_elements_rev = range_occs_rev[pos_start_extension+1].second - range_occs_rev[pos_start_extension+1].first + 1;
	// nb_elements_rev = min(range_occs_rev[pos_start_extension+1].second - range_occs_rev[pos_start_extension+1].first + 1,
	// 				  nb_disp);
	
	findElement = false;    
	for (ulong j=0 ; j < nb_elements_rev && !findElement ; j++){
	  if (locs_start_after[j] + 1 == loc_start_extension){
	    findElement = true;
	    pos_start_extension++;
	    loc_start_extension--;
	  }
	}
	
	if (locs_start_after)
	  free(locs_start_after);
	
	if (!findElement)
	  findStart = true;
      }
    }else{
      while (range_occs_fwd[pos_start_extension+1].first <=
	     range_occs_fwd[pos_start_extension+1].second && !findStart){
	getGenome()->
	  getOccurrencesFMIndexRange(range_occs_fwd[pos_start_extension+1].first,
				     range_occs_fwd[pos_start_extension+1].second,
				     &locs_start_after);
	
	nb_elements_fwd = range_occs_fwd[pos_start_extension+1].second - range_occs_fwd[pos_start_extension+1].first + 1;
	// nb_elements_fwd = min(range_occs_fwd[pos_start_extension+1].second - range_occs_fwd[pos_start_extension+1].first + 1,
	// 			    nb_disp);
	
	findElement = false;    
	for (ulong j=0 ; j < nb_elements_fwd && !findElement ; j++){
	  if (locs_start_after[j] - 1 == loc_start_extension){
	    findElement = true;
	    pos_start_extension++;
	    loc_start_extension++;
	  }
	}
	
	if (locs_start_after)
	  free(locs_start_after);
	
	if (!findElement)
	  findStart = true;
      }
    }
  }

  // Adjust the end break
  if (!hasNoEndExtension){
    if (strand_end_extension == -1){
      while (range_occs_rev[pos_end_extension-1].first <=
	     range_occs_rev[pos_end_extension-1].second && !findEnd){
	getGenome()->
	  getOccurrencesFMIndexRange(range_occs_rev[pos_end_extension-1].first,
				     range_occs_rev[pos_end_extension-1].second,
				     &locs_end_before);
	
	nb_elements_rev = range_occs_rev[pos_end_extension-1].second - range_occs_rev[pos_end_extension-1].first + 1;
	// nb_elements_rev = min(range_occs_rev[pos_end_extension-1].second - range_occs_rev[pos_end_extension-1].first + 1,
	// 				  nb_disp);
	
	findElement = false;    
	for (ulong j=0 ; j < nb_elements_rev && !findElement ; j++){
	  if (locs_end_before[j] - 1 == loc_end_extension){
	    findElement = true;
	    pos_end_extension--;
	    loc_end_extension++;
	  }
	}
	
	if (locs_end_before)
	  free(locs_end_before);
	
	if (!findElement)
	  findEnd = true;
      }
    }else if (strand_end_extension == 1){
      while (range_occs_fwd[pos_end_extension-1].first <=
	     range_occs_fwd[pos_end_extension-1].second && !findEnd){
	getGenome()->
	  getOccurrencesFMIndexRange(range_occs_fwd[pos_end_extension-1].first,
				     range_occs_fwd[pos_end_extension-1].second,
				     &locs_end_before);
	
	nb_elements_fwd = range_occs_fwd[pos_end_extension-1].second - range_occs_fwd[pos_end_extension-1].first + 1;
	// nb_elements_fwd = min(range_occs_fwd[pos_end_extension-1].second - range_occs_fwd[pos_end_extension-1].first + 1,
	// 			    nb_disp);
	
	findElement = false;    
	for (ulong j=0 ; j < nb_elements_fwd && !findElement ; j++){
	  if (locs_end_before[j] + 1 == loc_end_extension){
	    findElement = true;
	    pos_end_extension--;
	    loc_end_extension--;
	  }
	}
	
	if (locs_end_before)
	  free(locs_end_before);
	
	if (!findElement)
	  findEnd = true;
      }
    }
  }
    
  //////////////////// UPDATE NEW BONDARIES ////////////////////
  if (!hasNoStartExtension){
    candidat_chosen->setPosStartBreak(pos_start_extension+1);
    candidat_chosen->setLocStartBreak(loc_start_extension);
  }else{
    candidat_chosen->setPosStartBreak(0);
    candidat_chosen->setLocStartBreak(~0);
    candidat_chosen->setStrandStartBreak(0);  
  }

  if (!hasNoEndExtension){
    candidat_chosen->setPosEndBreak(pos_end_extension-1);
    candidat_chosen->setLocEndBreak(loc_end_extension);
  }else{
    candidat_chosen->setPosEndBreak(getLength() - 1);
    candidat_chosen->setLocEndBreak(~0);
    candidat_chosen->setStrandEndBreak(0);
  }
  

}


bool SupportBreak::canCheckPositions(ulong start, ulong end) {
  // We do not check the positions if we have a repetition just before
  // or just after the break.
  // A repetition artificially rises the score up which does not
  // allow us to make a difference between a mutation or an error.
  ulong nb_disp = getGenome()->getNbLocations();
  return ((range_occs_fwd[start].second - range_occs_fwd[start].first+1
  	   + range_occs_rev[start].second - range_occs_rev[start].first+1) <  min(nb_disp+1,parameters->min_occ_repetition)
  	  && (range_occs_fwd[end].second - range_occs_fwd[end].first +1 
  	      + range_occs_rev[end].second - range_occs_rev[end].first +1) <  min(nb_disp+1,parameters->min_occ_repetition)
  	  );
}


bool SupportBreak::check_nice_merge(SupportBreak *sb1, SupportBreak *sb2, 
                                    bool very_nice) {
  if (sb1 == sb2)
    return (very_nice) ? sb1->isVeryNiceBreak() : sb1->isNiceBreak();
  
  SupportBreak merged(*sb1, *sb2);
  
  bool result = (very_nice) ? merged.isVeryNiceBreak() : merged.isNiceBreak();

  return result;
}

void SupportBreak::computeHighAndLowAverages(uint pos_start, uint pos_end){
  // uint nb_outside = 0;
  uint pos_start_out, pos_end_out;
  if (pos_start < parameters->support_score_window_length/2)
    pos_start_out = 0;
  else
    pos_start_out = pos_start - parameters->support_score_window_length/2;
  if (support->getLength() 
      < 1 + pos_end + parameters->support_score_window_length/2)
    pos_end_out = support->getLength()-1;
  else
    pos_end_out = pos_end + parameters->support_score_window_length/2;
  
  //computing standard deviations    
  for (uint i = pos_start_out; i < pos_start ; i++) {
    if (support->getNbLocs()[i] == 0){
      //nb_outside=0;
      pos_start_out = i+1;
    }// else{
    //   nb_outside++;
    // }
  }

  // uint support_without_loc = 0;

  // // We only consider the support which are not located
  // // because it may bias the computation.
  // for (uint i = pos_start; i<=pos_end ; i++){
  //   if (support->getNbLocs()[i] == 0
  //       || (support_without_loc > 0 
  //           && support_without_loc >= support->getSupport(i))){
  //     nb_inside++;
  //   } else
  //     support_without_loc = support->getSupport(i);
  // } 
  // for (uint i=pos_end+1; i <= pos_end_out 
  // 	 && support->getNbLocs()[i] != 0; i++) {
  //   nb_outside++;
  // }
  
  //  uint *inside_elements = new uint[nb_inside];
  uint *inside_elements = new uint[pos_end-pos_start+1];

  // nb_outside = 0;
  // nb_inside = 0;

  // support_without_loc = 0;

  // We only consider the support which are not located
  // because it may bias the computation.
  uint nb_inside = 0;
  for (uint i = pos_start; i<=pos_end ; i++){
    if (support->getNbLocs()[i] == 0
        // || (support_without_loc > 0 
        //     && support_without_loc >= support->getSupport(i))
	){
      inside_elements[nb_inside] = support->getSupport(i);
      nb_inside++;
    } // else
      // support_without_loc = support->getSupport(i);
  } 
  // Sorting the positions so that we can take
  // the lowest ones and the highest ones.
  qsort(inside_elements, nb_inside, sizeof(uint), comparUint);
  
  // Get the first and last quartiles
  median_inside = inside_elements[nb_inside/2];
  inside_quartile1 = inside_elements[nb_inside/4];
  inside_quartile4 = inside_elements[3*nb_inside/4];

  average_low_inside = 0;
  // An error affects getThreshold() factors, therefore
  // if an error occurs inside a larger break, at least
  // getThreshold() factors should be impacted.
  // On the other hand, because of random locations, 
  // some factors containing the error may be located
  // and thus may not appear in the inside support.
  uint nb_low = min(getThreshold() - parameters->nb_bases_substracted_low_average,
                    nb_inside/2);
  for (uint i = 0; i < nb_low; i++) {
    average_low_inside += inside_elements[i];
  }
  average_high_inside = 0;
  for (uint i = nb_inside/2; i < nb_inside; i++) {
    average_high_inside += inside_elements[i];
  }

  average_high_inside /= nb_inside/2 + (nb_inside%2);
  if (nb_low == 0)
    average_low_inside = average_high_inside;
  else
    average_low_inside /= nb_low;

  score_inside_averages = (getThresholdForScoreIntraExon(average_high_inside) - average_low_inside);


  pc_ones_left_inside = pc_ones_left_outside = pc_ones_right_inside = pc_ones_right_outside = 0;

  // If the low average is 1 (or almost), we have to check if these low
  // values are due to a poorly supported read or to an error.
  // If it is due to an error, it is likely that we have 1s in the support
  // at a border inside the break, but higher values outside the break.
  if (average_low_inside < 1+EPSILON) {
    if (pos_start - pos_start_out > parameters->nb_positions_check_ones) {
      pos_start_out = pos_start - parameters->nb_positions_check_ones;
    }
    if (pos_end_out - pos_end > parameters->nb_positions_check_ones) {
      pos_end_out = pos_end + parameters->nb_positions_check_ones;
    }

    for (uint i = pos_start_out; i < pos_start; i++) {
        if (support->getSupport(i) == 1)
          pc_ones_left_outside++;
    }
    pc_ones_left_outside /= (pos_start - pos_start_out);

    uint count = 0;
    for (uint i = pos_start; count < parameters->nb_positions_check_ones 
           && i <= pos_end; i++) {
      if (support->getNbLocs()[i] == 0){
        if (support->getSupport(i) == 1)
          pc_ones_left_inside++;
        count++;
      }
    }
    pc_ones_left_inside /= count;
    
    count = 0;
    for (uint i = pos_end; i > pos_start
           && count < parameters->nb_positions_check_ones; i--) {
      if (support->getNbLocs()[i] == 0){
        if (support->getSupport(i) == 1)
          pc_ones_right_inside++;      
        count++;
      }
    }
    pc_ones_right_inside /= count;

    for (uint i = pos_end+1; i <= pos_end_out; i++) {
      if (support->getSupport(i) == 1)
        pc_ones_right_outside++;
    }
    pc_ones_right_outside /= (pos_end_out - pos_end);
  }

  delete [] inside_elements;

}

void SupportBreak::fillWithClosestMatch() {
  ulong *pos_occ;
  
  if (hasNoStartBreak()) {
    // We have no localisation before the break (beginning of the tag)
    uint pos_after_break = pos_end_break+1;
    loc_start_break = ~0;
    strand_start_break = 0;
    // getGenome()->setNbLocations(1);
  
    // There are localisations on the forward strand
    if (range_occs_fwd[pos_after_break].first <=
   	range_occs_fwd[pos_after_break].second) {
      getGenome()->
   	getOccurrencesFMIndexRange(range_occs_fwd[pos_after_break].first,
   				   range_occs_fwd[pos_after_break].second,
   				   &pos_occ,0,1);
      strand_end_break = 1;
    } else {
      // Localisations are on the reverse strand
      getGenome()->
   	getOccurrencesFMIndexRange(range_occs_rev[pos_after_break].first,
   				   range_occs_rev[pos_after_break].second,
   				   &pos_occ,0,1);
      strand_end_break = -1;
    }
    loc_end_break = pos_occ[0];
    if (pos_occ)
      free(pos_occ);
  
    // getGenome()->setNbLocations(nb_disp);
  } else if (hasNoEndBreak()) {
    // We have no localisation before the break (beginning of the tag)
    loc_end_break = ~0;
    strand_end_break = 0;
    uint pos_before_break = pos_start_break-1;
  
    // getGenome()->setNbLocations(1);
  
    // There are localisations on the forward strand
    if (range_occs_fwd[pos_before_break].first <= 
   	range_occs_fwd[pos_before_break].second) {
      getGenome()->
   	getOccurrencesFMIndexRange(range_occs_fwd[pos_before_break].first,
   				   range_occs_fwd[pos_before_break].second,
   				   &pos_occ,0,1);
      strand_start_break = 1;
    } else {
      // Localisations are on the reverse strand
      getGenome()->
   	getOccurrencesFMIndexRange(range_occs_rev[pos_before_break].first,
   				   range_occs_rev[pos_before_break].second,
   				   &pos_occ,0,1);
      strand_start_break = -1;
    }
    loc_start_break = pos_occ[0];
    if (pos_occ)
      free(pos_occ);
  
    // getGenome()->setNbLocations(nb_disp);

  }
  // Both localisations before and after break are defined
  else if (!hasNoEndBreak() && !hasNoStartBreak()) {

    uint pos_end_extension = getPositionEndBreak();
    uint pos_start_extension = getPositionStartBreak();

    bool no_loc_before=false;
    bool no_loc_after=false;

    // from the break, we try to extend until max_extension_length
    // while it is not on border case or it is not the best_candidate 
    // or it is not around a repetition
    for (int i=0 ; 
	 (i <= parameters->max_extension_length)
	   && (pos_end_extension < getLength()-1
		 || pos_start_extension > 0)
	   && (!no_loc_after || !no_loc_before)
	   && ! repeated ; 
	 i++) {
           

      ///////////////////////////////////////////////////////////////////////
      /////////// CHECKING THE POSSIBILITY TO EXTEND FOR EACH SIDE //////////
      ///////////////////////////////////////////////////////////////////////
      

      // On the right, the extension is not possible if:
      // i/ it is the read's end
      if (pos_end_extension < (getLength() - 1) && ! no_loc_after)
	pos_end_extension++;      
      else
	no_loc_after = true;
      // ii/ the extension is not coming inside an other break
      if ((range_occs_fwd[pos_end_extension].first > range_occs_fwd[pos_end_extension].second)
	  && (range_occs_rev[pos_end_extension].first > range_occs_rev[pos_end_extension].second)) {
	no_loc_after = true;
	pos_end_extension--;
      }  

      // On the right, the extension is not possible if:
      // i/ it is the read's start
      if (pos_start_extension > 0 && ! no_loc_before)
	pos_start_extension--;
      else
	no_loc_before = true;
      // ii/ the extension is not coming inside an other break
      if ((range_occs_fwd[pos_start_extension].first > range_occs_fwd[pos_start_extension].second)
	  && (range_occs_rev[pos_start_extension].first > range_occs_rev[pos_start_extension].second)){
	no_loc_before = true;
	pos_start_extension++;
      }

      ///////////////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////////////
      //////////// RUNNING UPDATELOCATION AND CHECKING CANDIDATS ////////////
      ///////////////////////////////////////////////////////////////////////

      if (!no_loc_before || !no_loc_after) {
        if (! canCheckPositions(pos_start_extension, pos_end_extension)) {
	  repeated = true;
        }
        // only it is not around a repetition
        else{
	  updateLocation(pos_start_extension, pos_end_extension);
        }
      }
           
      ///////////////////////////////////////////////////////////////////////
      
    }
    

    ///////////////////////////////////////////////////////////////////////
    /////////////////// CH0OSE THE BEST CANDIDAT///////////////////////////
    ///////////////////////////////////////////////////////////////////////

    if (nb_candidats > 0){
      candidat_chosen = new CandidatBreak(*candidats[getBestCandidat()]);

      // a candidat is necessarily chosen
      if (candidat_chosen == NULL)
  throw ILLEGAL_STATE_EXCEPTION;
      

      // a candidat is chosen and (max_extension_length == 0), so we just must adjust PosStartBreak and posEndBreak
      if (parameters->max_extension_length == 0){
	candidat_chosen->setPosStartBreak(pos_start_extension+1);
	candidat_chosen->setPosEndBreak(pos_end_extension-1);   
      }
      // a candidat is chosen, so its boundaries are updated 
      else
	adjustBoundaries();
    }
    // nb_candidats == 0 if there are directly a repetion at the first step of pos_start and pos_end breaks both.    
    else{
      nb_candidats = 0;
    }
  }
}


 /**
  * Return the position of the best_candidat :
  * - nice_candidat (at least a splice)
  * - single consistent
  * - the minimal break on the genome
  * And notify a duplicated candidat if the table contains several good candidats or (several not godd candidats)
  */
 int SupportBreak::getBestCandidat(CandidatBreak **potential_candidats, int nb_potential_candidats, bool after_extension){
   int i_chosen = 0;
   for (int i = 1 ; i < nb_potential_candidats ; i++){
    // current candidat is a chimera
    if (!potential_candidats[i_chosen]->isNiceCandidat()){
      // new candidat is at least a splice, so it is better
      if (potential_candidats[i]->isNiceCandidat()){
	i_chosen = i;
      }
      // at least two different chimeras
      else{
	// in case of after extension, a same chimera is automatically extended  
	// but only if the new candidat is not duplicated or the two candidats are duplicated both
	if (after_extension && (!potential_candidats[i]->isDuplicated() 
				//|| potential_candidats[i_chosen]->isDuplicated()
				)
	    ){
	  i_chosen = i;
	}
	// normal process in case of not after extension
	else{
	  if (!isSameCandidat(potential_candidats,i,i_chosen)){
	    potential_candidats[i_chosen]->setDuplicated(true);
	    potential_candidats[i]->setDuplicated(true);
	
	    // new candidat is also a chimera but it is consistent with the 
	    // single map, so it is better 
	    if (!potential_candidats[i_chosen]->isSingleConsistent() && potential_candidats[i]->isSingleConsistent()){
	      i_chosen = i;
	    }
	    // two different chimeras (both single or both not), so: 
	    // In case of after the extension process, the good one depends of getBreakLength
	    else if (!potential_candidats[i_chosen]->isSingleConsistent() || potential_candidats[i]->isSingleConsistent()) {
	      if (after_extension){
		if (potential_candidats[i]->getReadBreakLength() < getThreshold() &&  hasLongEnoughTagBreak()){
		  i_chosen = i;
		}
	      }
	      // else we try to do a choice
	      else{
		// we keep the candidat with the smallest single_distance
		if (potential_candidats[i]->isSingleConsistent() 
		    && potential_candidats[i]->getSingleDistance() < potential_candidats[i_chosen]->getSingleDistance()) {
		  i_chosen = i;
		}
		// or we keep the candidat with the smallest gap on the genome
		else if ((
			  !potential_candidats[i]->isSingleConsistent() 
			  || (potential_candidats[i]->getSingleDistance() == potential_candidats[i_chosen]->getSingleDistance())
			  )
			 && 
			 (
			  labs((long long int) potential_candidats[i]->getGenomeGapLength()-potential_candidats[i]->getReadBreakLength()) 
			  < labs((long long int)(potential_candidats[i_chosen]->getGenomeGapLength()-potential_candidats[i_chosen]->getReadBreakLength()))
			  )
			 ){
		  i_chosen = i;
		}
	      }
	    }
	  }
	}
      } // end at least two different chimeras
    }
    // current candidat is at least a splice
    else{
      // new candidat is also at least a splice
      if (potential_candidats[i]->isNiceCandidat()){
	// new candidat is also a "NiceCandidat" but it is consistent with the 
	// single map, so it is better 
	if (!potential_candidats[i_chosen]->isSingleConsistent() && potential_candidats[i]->isSingleConsistent()){
	  i_chosen = i;
 	}
	// two different "NiceCandidat" (both single or both not), so it is a duplicated cause
	else if (!potential_candidats[i_chosen]->isSingleConsistent() || potential_candidats[i]->isSingleConsistent()) {
	  if (!isSameCandidat(potential_candidats,i,i_chosen)){
	     potential_candidats[i_chosen]->setDuplicated(true);
	     potential_candidats[i]->setDuplicated(true);
	     
	     // if both single consistency, we keep the candidat with the smallest single_distance
	     if (potential_candidats[i]->isSingleConsistent() 
		 && potential_candidats[i]->getSingleDistance() < potential_candidats[i_chosen]->getSingleDistance()) {
	       i_chosen = i;
	     }
	     // we keep the candidat with the smallest gap on the genome
	     else if ((
		       !potential_candidats[i]->isSingleConsistent() 
		       || (potential_candidats[i]->getSingleDistance() == potential_candidats[i_chosen]->getSingleDistance())
		       )
		      && 
		      (
		       (potential_candidats[i]->getGenomeGapLength()-potential_candidats[i]->getReadBreakLength()) 
		       < (potential_candidats[i_chosen]->getGenomeGapLength()-potential_candidats[i_chosen]->getReadBreakLength())
		       )
		      ){
	       i_chosen = i;
	     }
	  }
	}
      }
    }
  }
  // in case of after extension, we have not a good candidat and we have extend to the bitter end or to the bitter start
  // so we keep the last candidat to easily check no_start and no_end break in adjust process 
  if (after_extension && !potential_candidats[i_chosen]->isNiceCandidat() 
      && (potential_candidats[nb_potential_candidats-1]->getPosEndBreak() == (getLength() - 1) || potential_candidats[nb_potential_candidats-1]->getPosStartBreak() == 0)
      ) {
    i_chosen = nb_potential_candidats-1;
  }
  
  return i_chosen;
}

int SupportBreak::getBestCandidat(){
  return getBestCandidat(this->candidats, this->nb_candidats, true);
}

bool SupportBreak::isSameCandidat(CandidatBreak **potential_candidats,ulong pos_i, ulong pos_j){
  return isSameCandidatEnd(potential_candidats, pos_i, pos_j) && isSameCandidatStart(potential_candidats, pos_i, pos_j);
}

bool SupportBreak::isSameCandidatEnd(CandidatBreak **potential_candidats,ulong pos_i, ulong pos_j){
  long long int diff_end = labs((long long int) potential_candidats[pos_j]->getPosEndBreak() - potential_candidats[pos_i]->getPosEndBreak());
  // Two candidats are not the same if there are differents strand
  if (potential_candidats[pos_i]->getStrandEndBreak() != potential_candidats[pos_j]->getStrandEndBreak())
    return false;
  // or different chr  
  if (potential_candidats[pos_i]->getChrIdEndBreak() != potential_candidats[pos_j]->getChrIdEndBreak())
    return false;
  // or the same with an extension step (or not if the extension is not longer possible before or after the break)
  if (
      (
       (potential_candidats[pos_i]->getPosEndBreak() != potential_candidats[pos_j]->getPosEndBreak() 
	|| potential_candidats[pos_i]->getLocEndBreak() != potential_candidats[pos_j]->getLocEndBreak()
	)
       && labs((long long int) potential_candidats[pos_i]->getLocEndBreak() - potential_candidats[pos_j]->getLocEndBreak()) != diff_end
       )
      )
    return false;
  // else two candidats are the same
  return true;
}

bool SupportBreak::isSameCandidatStart(CandidatBreak **potential_candidats,ulong pos_i, ulong pos_j){
  long long int diff_start = labs((long long int) potential_candidats[pos_j]->getPosStartBreak() - potential_candidats[pos_i]->getPosStartBreak());
  // Two candidats are not the same if there are differents strand
  if (potential_candidats[pos_i]->getStrandStartBreak() != potential_candidats[pos_j]->getStrandStartBreak())
    return false;
  // or different chr  
  if (potential_candidats[pos_i]->getChrIdStartBreak() != potential_candidats[pos_j]->getChrIdStartBreak())
    return false;
  // or the same with an extension step (or not if the extension is not longer possible before or after the break)
  if (
      (
       (potential_candidats[pos_i]->getPosStartBreak() != potential_candidats[pos_j]->getPosStartBreak() 
	|| potential_candidats[pos_i]->getLocStartBreak() != potential_candidats[pos_j]->getLocStartBreak()
	) 
       && labs((long long int) potential_candidats[pos_i]->getLocStartBreak() - potential_candidats[pos_j]->getLocStartBreak()) != diff_start
       ) 
      )
    return false;
  // else two candidats are the same
  return true;
}

void SupportBreak::saveCandidates(uint pos_before_break, uint pos_after_break,
				  ulong *pos_occ_start, uint nb_loc_start, 
				  int strand_start, ulong *pos_occ_end, 
				  uint nb_loc_end, int strand_end, CandidatBreak ***potential_candidats, ulong &nb_potential, ulong &nb_potential_max) {
    
  uint i, j;
  
  if (nb_potential == 0)
    *potential_candidats = (CandidatBreak**)malloc(sizeof(CandidatBreak*)*nb_potential_max); 
  
  uint loc_occ_single = ~0;
  int strand_occ_single = 0;
  uint pos_occ_single = ~0;
  
  // get single_loc features
  if (support->getBestLocation() != NULL){
    loc_occ_single = getGenome()->getAbsolutePosition(support->getBestLocation());
    strand_occ_single = support->getBestLocation()->getStrand();
    pos_occ_single = support->getPositionOfLocation();
  }
  
  ///////// TEMPORARY SAVE ALL CANDIDATS FOR AN EXTENSION STEP //////////

  i = 0; j = 0;
  while (i < nb_loc_start && j < nb_loc_end) {
    if (nb_potential >= nb_potential_max) {
      nb_potential_max *= 2;
      *potential_candidats = (CandidatBreak**)realloc(*potential_candidats, sizeof(CandidatBreak*)*nb_potential_max);
      
    }
    CandidatBreak *temp_candidat = new CandidatBreak(pos_before_break, pos_occ_start[i], strand_start, pos_after_break, pos_occ_end[j], strand_end, this, parameters);
    temp_candidat->checkSingleCorrespondence(pos_occ_single, loc_occ_single, strand_occ_single);
    (*potential_candidats)[nb_potential] = temp_candidat;
    nb_potential++;
    if (i < (nb_loc_start - 1))
      i++;
    else{ 
      j++;
      i=0;
    }
  }
}



void SupportBreak::updateLocation(uint pos_before_break,
				  uint pos_after_break){
  
  // maximum number of locations defined in parameter by user
  //ulong nb_disp = getGenome()->getNbLocations();
    
  ulong *pos_occ_start_fwd, *pos_occ_end_fwd;
  ulong *pos_occ_start_rev, *pos_occ_end_rev;
  uint nb_loc_start_fwd=0, nb_loc_end_fwd=0;
  uint nb_loc_start_rev=0, nb_loc_end_rev=0;
 
  CandidatBreak **potential_candidats;
  ulong nb_potential_candidats=0;
  ulong max_potential_candidats = parameters->min_occ_repetition;
  
  bool has_start_forward_strand = false; 
  bool has_start_reverse_strand = false;
  bool has_end_forward_strand = false; 
  bool has_end_reverse_strand = false;

  // Locating at start forward
  if (range_occs_fwd[pos_before_break].first <=
      range_occs_fwd[pos_before_break].second){
    has_start_forward_strand = true;
    getGenome()->
      getOccurrencesFMIndexRange(range_occs_fwd[pos_before_break].first,
				 range_occs_fwd[pos_before_break].second,
				 &pos_occ_start_fwd);
    nb_loc_start_fwd = range_occs_fwd[pos_before_break].second - range_occs_fwd[pos_before_break].first+1;
  }

  // Locating at end forward
  if (range_occs_fwd[pos_after_break].first <= 
      range_occs_fwd[pos_after_break].second) {
    has_end_forward_strand = true;
    getGenome()->
      getOccurrencesFMIndexRange(range_occs_fwd[pos_after_break].first,
				 range_occs_fwd[pos_after_break].second,
				 &pos_occ_end_fwd);
    nb_loc_end_fwd = range_occs_fwd[pos_after_break].second - range_occs_fwd[pos_after_break].first+1;
  }

  // Locating at start reverse
  if (range_occs_rev[pos_before_break].first <= 
      range_occs_rev[pos_before_break].second){
    has_start_reverse_strand = true;
    getGenome()->
      getOccurrencesFMIndexRange(range_occs_rev[pos_before_break].first,
				 range_occs_rev[pos_before_break].second,
				 &pos_occ_start_rev);
    nb_loc_start_rev = range_occs_rev[pos_before_break].second - range_occs_rev[pos_before_break].first+1;
  }

  // Locating at end reverse
  if (range_occs_rev[pos_after_break].first <= 
      range_occs_rev[pos_after_break].second) {
    has_end_reverse_strand = true;
    getGenome()->
      getOccurrencesFMIndexRange(range_occs_rev[pos_after_break].first,
				 range_occs_rev[pos_after_break].second,
				 &pos_occ_end_rev);
    nb_loc_end_rev = range_occs_rev[pos_after_break].second - range_occs_rev[pos_after_break].first+1;
  }

  ///// A GOOD CANDIDAT IS PROBABLE IF START/END ARE ON THE SAME STRAND ////
  if ((has_start_forward_strand && has_end_forward_strand)
      || (has_start_reverse_strand && has_end_reverse_strand)){
    //save candidats on forward strand
    if (has_start_forward_strand && has_end_forward_strand){
      saveCandidates(pos_before_break, pos_after_break,
		     pos_occ_start_fwd, nb_loc_start_fwd,1, 
		     pos_occ_end_fwd, nb_loc_end_fwd,1,
		     &potential_candidats, nb_potential_candidats,max_potential_candidats); 
    }
    //save candidats on reverse strand
    if (has_start_reverse_strand && has_end_reverse_strand){
      saveCandidates(pos_before_break, pos_after_break,
		     pos_occ_start_rev, nb_loc_start_rev,-1, 
		     pos_occ_end_rev, nb_loc_end_rev,-1,
		     &potential_candidats, nb_potential_candidats, max_potential_candidats);  
    }
  }
  
  // test if there is a good candidat (at least a splice)
  bool isIntermediarySplice = false;
  if (nb_potential_candidats > 0){
    int i_chosen = getBestCandidat(potential_candidats, nb_potential_candidats);
    if (i_chosen != 0){
      delete potential_candidats[0];
      potential_candidats[0] = new CandidatBreak(*potential_candidats[i_chosen]);
    }
    for (uint i=1; i < nb_potential_candidats; i++) {
      delete potential_candidats[i];
    }
    nb_potential_candidats = 1;
    isIntermediarySplice = potential_candidats[0]->isNiceCandidat();
  }
  
  /////////////////// WE SAVE SOME CHIMERAS //////////////////
  // get other chimeras if there is not at least an intermediary splice 
  if (nb_potential_candidats == 0 || !isIntermediarySplice){ 
    if (has_start_forward_strand && has_end_reverse_strand){
      saveCandidates(pos_before_break, pos_after_break,
		     pos_occ_start_fwd, nb_loc_start_fwd,1, 
		     pos_occ_end_rev, nb_loc_end_rev,-1,
		     &potential_candidats, nb_potential_candidats, max_potential_candidats);
    }  
    if (has_start_reverse_strand && has_end_forward_strand){
      saveCandidates(pos_before_break, pos_after_break,
		     pos_occ_start_rev, nb_loc_start_rev,-1, 
		     pos_occ_end_fwd, nb_loc_end_fwd,1,
		     &potential_candidats, nb_potential_candidats, max_potential_candidats);  
    }
  }

  ////////////////////////////////////////////////////////////////////////
  /////////////////// SAVE THE BEST CANDIDAT AT THIS STEP /////////////////
  ////////////////////////////////////////////////////////////////////////

  if (nb_potential_candidats > 0){
    int i_chosen = 0;
    // a new candidat is chosen only if it is not an at least intermediary splice
    if (!isIntermediarySplice && nb_potential_candidats > 1){
      i_chosen = getBestCandidat(potential_candidats, nb_potential_candidats);
    }
    // we have init nb_candidats == -1 in order to distinguish init case and no-candidat case
    if (nb_candidats == -1)
      nb_candidats=0;

    candidats[nb_candidats] = new CandidatBreak(*potential_candidats[i_chosen]);
    nb_candidats++;
  }
  
  /////////////////// DELETE THE TEMPORARY TABLE  /////////////////////////  
 
  if (has_start_forward_strand) 
    free(pos_occ_start_fwd);
  if (has_end_forward_strand) 
    free(pos_occ_end_fwd);
  if (has_start_reverse_strand)
    free(pos_occ_start_rev);
  if (has_end_reverse_strand)
    free(pos_occ_end_rev);
  
  for (uint i=0; i < nb_potential_candidats; i++) {
    delete potential_candidats[i];
  }
  if (nb_potential_candidats>0)
    free(potential_candidats);
}


ostream &operator<<(ostream &os, SupportBreak &sb) {
  os << "[break={" << sb.getPositionStartBreak() << "," 
     << sb.getPositionEndBreak() << "};gap={" 
     << sb.getPositionStartBreak() << "," 
     << sb.getPositionEndBreak() << "};loc={";
  if (sb.getLocationStartBreak())
    os << *sb.getLocationStartBreak();
  else
    os << "undef";
  os << ",";
  if (sb.getLocationEndBreak())
    os << *sb.getLocationEndBreak();
  else
    os << "undef";
  os << "}]" << endl;
  return os;
}
