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
#include "tagInfo.h"
#include "../const.h"
#include "../libSSA/utils.h"
#include "../utils.h"
#include "Cigar.h"
#include <limits.h>
#include <string>
#include <sstream>
#include <assert.h>     /* assert */

TagInfo::TagInfo(LocateOnGenome *genome, Support *s, Read *r)
  :code(0),genome(genome), nb_each_type(NULL),info_each_type(NULL), nb_causes(0),
   support(s), current_break(0),
  read(r) {

  if (support->isNone())
    addNone();
  if (support->isDuplicate())
    addDuplication();
  if (support->isMultiple())
    addMultiple();
  if (support->isSingle()){
    addSingle();
    if (support->hasRepeat()) 
      addRepetition(support->getStartPosRepeat(),
                    support->getEndPosRepeat());
    else if (support->isContinuous()) {
      if (! support->isAlmostNormal()){
  	addNormal();
      }else{
  	addAlmostNormal();
      }
    }
  }

  chrpos[START_BREAK] = (ChrPosition **)calloc(support->getNbBreaks(), 
                                               sizeof(ChrPosition*));
  chrpos[END_BREAK] = (ChrPosition **)calloc(support->getNbBreaks(), 
                                             sizeof(ChrPosition*));
  is_good_location = true;
  }


TagInfo::~TagInfo() {
  if (nb_each_type != NULL) {
    for (uint i=0; i < INFO_NB_TYPES; i++) {
      for (uint j=0; j < nb_each_type[i]; j++) {
        delete info_each_type[i][j];
      }
      if (nb_each_type[i])
        free(info_each_type[i]);
    }
    free(info_each_type);
    free(types_in_order);
    free(nb_causes_per_break);
  }
  for (uint j = 0; j < NB_POS_BREAK; j++) {
    for (uint i=0; i < support->getNbBreaks(); i++) 
      if (chrpos[j][i]) 
        delete (chrpos[j][i]);
    free(chrpos[j]);
  }

  if (nb_each_type != NULL)
    free(nb_each_type);
}

void TagInfo::addGenericElement(byte type, GenericInfo *gi) {
  if (nb_each_type == NULL) {
    nb_each_type = (byte *)calloc(INFO_NB_TYPES,sizeof(byte));
    info_each_type = (GenericInfo ***)malloc(sizeof(GenericInfo**)*INFO_NB_TYPES);
    types_in_order = (byte *)malloc(sizeof(byte) * NB_INFO_STORED);
    nb_causes_per_break = (byte *)calloc(NB_INFO_BREAKS+1, sizeof(byte));
    nb_causes = 0;
  }

  if (nb_each_type[type] == 0) 
    info_each_type[type] = (GenericInfo **)malloc(NB_ELEMENTS_PER_TYPE * sizeof(gi));
  else if (nb_each_type[type] % NB_ELEMENTS_PER_TYPE == 0)
    info_each_type[type] = (GenericInfo **)realloc(info_each_type[type]
                                                   , (nb_each_type[type]
                                                      +NB_ELEMENTS_PER_TYPE)*sizeof(gi));
    
  if (type != INFO_REPETITION && type != INFO_DUPLICATION) { 
    if (current_break > 0 && current_break % NB_INFO_BREAKS == 0
	&& nb_causes_per_break[current_break] == 0) {
      nb_causes_per_break = (byte *)realloc(nb_causes_per_break,
                                            (current_break + NB_INFO_BREAKS+1) * sizeof(byte));
      for (uint i=current_break ; i <= current_break + NB_INFO_BREAKS; i++)
        nb_causes_per_break[i] = 0;
    }
    nb_causes_per_break[current_break]++;

    if (nb_causes > 0 && nb_causes % NB_INFO_STORED == 0) {
      types_in_order = (byte *) realloc(types_in_order, 
                                        (nb_causes + NB_INFO_STORED) * sizeof(byte));
    }

    gi->setBreakId(current_break);
    gi->setCauseId(nb_causes);

    types_in_order[nb_causes++] = type;
  }

#ifdef DBG
  cerr << "Adding element " <<(int)type<<" for tag " << num << endl;
#endif

  gi->setReadLength(read->getLength());
  info_each_type[type][nb_each_type[type]] = gi;
  nb_each_type[type]++;
}

byte TagInfo::getNbGenericElement(uint type) {
  if (nb_each_type == NULL)
    return 0;
  return nb_each_type[type];
}

GenericInfo **TagInfo::getGenericInfos(uint type) {
  if (nb_each_type == NULL || nb_each_type[type] == -2)
    return NULL;
  return info_each_type[type];
}

void TagInfo::changeGenericElement(byte type, byte newType, int oldMask, 
                                   int newMask,uint i, GenericInfo *newGi) {
  GenericInfo *oldGi = info_each_type[type][i];
  uint causeId = oldGi->getCauseId();

  // We change the type of the element in types in order
  types_in_order[causeId] = newType;

  // Set read length breakId and causeId from the old element
  newGi->setReadLength(oldGi->getReadLength());
  newGi->setBreakId(oldGi->getBreakId());
  newGi->setCauseId(oldGi->getCauseId());

  // Delete the old Generic Element
  delete oldGi;
  memmove(&(info_each_type[type][i]), &(info_each_type[type][i+1]), (nb_each_type[type] - i - 1)*sizeof(oldGi));
  if (nb_each_type[type] == 1)
    free(info_each_type[type]);
  nb_each_type[type]--;

  // update code: we have a new type
  code |= newMask;
  // And we may remove the previous one (if there is no more)
  if (nb_each_type[type] == 0)
    code ^= oldMask;

  // Insert the new Generic Element in the right place of InfoEachType
  // Update info_each_type size if needed
  if (nb_each_type[newType] == 0) {
    info_each_type[newType] = (GenericInfo **)malloc(NB_ELEMENTS_PER_TYPE * sizeof(newGi));
    info_each_type[newType][nb_each_type[newType]] = newGi;
  } else { 
    if (nb_each_type[newType] % NB_ELEMENTS_PER_TYPE == 0)
    info_each_type[newType] = (GenericInfo **)realloc(info_each_type[newType]
                                                   , (nb_each_type[newType]
                                                      +NB_ELEMENTS_PER_TYPE)*sizeof(newGi));

    // Count elements of newType in array types in order to get the place to insert the new element 
    uint nbNewType = 0;
    for(uint j=0; j<causeId; j++) {
      if(types_in_order[j] == newType) {
        nbNewType++;
      }
    }
    // Move the element of new type in order to insert the new generic info in the right place
    memmove(&(info_each_type[newType][nbNewType+1]), &(info_each_type[newType][nbNewType]), (nb_each_type[newType] - nbNewType)*sizeof(newGi));
    // Delete the element duplicated by the memmove
    //if (nb_each_type[newType]  > 1)
    //  delete info_each_type[newType][nbNewType];
    info_each_type[newType][nbNewType] = newGi;
  }
  nb_each_type[newType]++;
}

void TagInfo::changeGenericElement(byte type, byte newType, int oldMask, 
                                   int newMask, uint break_id, uint cause_num, 
                                   GenericInfo *newGi) {
  for (uint i = 0; i < nb_each_type[type]; i++) {
    if (info_each_type[type][i]->getBreakId() == break_id) {
      // We found a cause in the same break
      if (cause_num == 0) {
        // Additionnally it is the same cause in this break.
        // We found what we were looking for!
        changeGenericElement(type, newType, oldMask, newMask, i, newGi);
        return;
      } else {
        cause_num--;
      }
    }
  }
}

int TagInfo::getCode() {
  return code;
}

float TagInfo::getAverageLowInside(uint i) {
  return support->getBreak(i)->getAverageLowInside();
}

float TagInfo::getAverageHighInside(uint i) {
  return support->getBreak(i)->getAverageHighInside();
}

SupportBreak *TagInfo::getBreak(uint i, int strand) {
  if (strand == 1)
    return support->getBreak(i);
  return support->getBreak(getNbBreaks() - i - 1);
}

byte TagInfo::getCauseType(uint i, int strand) {
  if (strand == 1)
    return types_in_order[i];
  else
    return types_in_order[getNbCauses() - i - 1];
}

ChrPosition *TagInfo::getChrPosEndBreak(uint i) {
  bool new_value = ! chrpos[END_BREAK][i] 
    || getBreak(i)->getPositionLocation(END_BREAK) != chrpos[END_BREAK][i]->getRelativePosition();
  if (new_value) {
    if (chrpos[END_BREAK][i]) 
      delete chrpos[END_BREAK][i];
    chrpos[END_BREAK][i] = support->getBreak(i)->getChrPosition(END_BREAK);
  }
  return chrpos[END_BREAK][i];
}

ChrPosition *TagInfo::getChrPosStartBreak(uint i) {
  bool new_value = ! chrpos[START_BREAK][i] 
    || getBreak(i)->getPositionLocation(START_BREAK) != chrpos[START_BREAK][i]->getRelativePosition();
  if (new_value) {
    if (chrpos[START_BREAK][i]) 
      delete chrpos[START_BREAK][i];
    chrpos[START_BREAK][i] = support->getBreak(i)->getChrPosition(START_BREAK);
  }
  return chrpos[START_BREAK][i];
}

uint TagInfo::getCurrentBreak() {
  return current_break;
}

int TagInfo::getEndPosRepeat() {
  return support->getEndPosRepeat();
}

gap_size_t TagInfo::getGenomeGapLength(uint i) {
  return support->getBreak(i)->getGenomeGapLength();
}

// bool TagInfo::getIsExtended(uint i) {
//   return is_extended[i];
// }

uint *TagInfo::getLocalisations() {
  return support->getNbLocs();
}

ChrPosition *TagInfo::getLocation() {
  return support->getBestLocation();
}

uint TagInfo::getLocationEndBreak(uint i) {
  return support->getBreak(i)->getPositionLocation(END_BREAK);
}

uint TagInfo::getLocationStartBreak(uint i) {
  return support->getBreak(i)->getPositionLocation(START_BREAK);
}

uint TagInfo::getNbBreaks() {
  return support->getNbBreaks();;
}

int TagInfo::getNbCandidats(uint i) {
  return support->getBreak(i)->getNbCandidats();
}

uint TagInfo::getNbDuplicate() {
  return support->getNbDuplicate();
}

uint TagInfo::getNbMultiple() {
  return support->getNbMultiple();
}

uint TagInfo::getNbSingle() {
  return support->getNbSingle();
}

uint TagInfo::getPositionEndBreak(uint i, int strand) {
  if (strand == 1)
    return support->getBreak(i)->getPositionEndBreak();
  return support->getBreak(getNbBreaks() - i - 1)->getPositionEndBreak(strand);
}

uint TagInfo::getPositionOfLocation(){
  return support->getPositionOfLocation();
}

uint TagInfo::getPositionStartBreak(uint i, int strand) {
  if (strand == 1)
    return support->getBreak(i)->getPositionStartBreak();
  return support->getBreak(getNbBreaks() - i - 1)->getPositionStartBreak(strand);
}

Read *TagInfo::getRead() {
  return read;
}

float TagInfo::getScoreComputedInterExon(uint i) {
  return support->getBreak(i)->getScoreComputedInterExon();
}

float TagInfo::getScoreComputedIntraExon(uint i) {
  return support->getBreak(i)->getScoreComputedIntraExon();
}

float TagInfo::getScoreInsideAverages(uint i) {
  return support->getBreak(i)->getScoreInsideAverages();
}

float TagInfo::getScoreInsideBreak(uint i) {
  return support->getBreak(i)->getScoreInsideBreak();
}

float TagInfo::getScoreOutsideBreak(uint i) {
  return support->getBreak(i)->getScoreOutsideBreak();
}

int TagInfo::getStartPosRepeat() {
  return support->getStartPosRepeat();;
}

uint *TagInfo::getSupport() {
  return support->getSupport();
}

Support *TagInfo::getSupportObject() {
  return support;
}

uint TagInfo::getSupportLength() {
  return support->getLength();
}

ulong TagInfo::getTagPosition() {
  return getLocation()->getRelativePosition();
}

char *TagInfo::getTagChromosome() {
  return getLocation()->getChrPosition();
}

uint TagInfo::getTagBreakLength(uint i) {
  return support->getBreak(i)->getTagBreakLength();
}

bool TagInfo::getTagStrand() {
  return getLocation()->getStrand();
}

uint TagInfo::getThreshold() {
  return support->getThreshold();
}

/* bool functions */
bool TagInfo::isAlmostNormal() {
  return (code & MASK_ALMOST_NORMAL);
}

bool TagInfo::isDuplicated(uint i) {
  return support->getBreak(i)->isDuplicated();
}

bool TagInfo::isExplainable() {
  return (code & (MASK_ALMOST_NORMAL | MASK_NORMAL | MASK_DUPLICATION | MASK_REPETITION 
		  | MASK_SNP | MASK_SEQ_ERR | MASK_SPLICE
		  | MASK_INTRA_TRANSPLICING | MASK_INTER_TRANSPLICING
		  | MASK_BIOLOGICAL_TAG_INDEL | MASK_BIOLOGICAL_UNDETERMINED));
}

bool TagInfo::isNormal() {
  return (code & MASK_NORMAL);
}

bool TagInfo::isMultiple() {
  return (code & MASK_MULTIPLE);
}

bool TagInfo::isNone() {
  return (code & MASK_NONE);
}

bool TagInfo::isRepeated(uint i) {
  return support->getBreak(i)->isRepeated();
}

bool TagInfo::isSingle() {
  return (code & MASK_SINGLE);
}

bool TagInfo::isSupportFallingLeft(uint i, int strand) {
  return support->getBreak(i)->isSupportFallingLeft(strand);
}

bool TagInfo::isSupportFallingRight(uint i, int strand) {
  return support->getBreak(i)->isSupportFallingRight(strand);
}

bool TagInfo::hasRepetition() {
  return (code & MASK_REPETITION);
}

bool TagInfo::isDuplication() {
  return (code & MASK_DUPLICATION);
}

bool TagInfo::hasSNP() {
  return (code & MASK_SNP);
}

bool TagInfo::hasBioTagIndel() {
  return (code & MASK_BIOLOGICAL_TAG_INDEL);
}

bool TagInfo::hasBioUndetermined() {
  return (code & MASK_BIOLOGICAL_UNDETERMINED);
}


bool TagInfo::hasSeqErr() {
  return (code & MASK_SEQ_ERR);
}


bool TagInfo::hasSplice() {
  return (code & MASK_SPLICE);
}

bool TagInfo::hasSpliceIntraChr() {
  return (code & MASK_INTRA_TRANSPLICING);
}

bool TagInfo::hasSpliceInterChr() {
  return (code & MASK_INTER_TRANSPLICING);
}

bool TagInfo::hasSpliceNoCover() {
  return (code & MASK_SPLICE_NO_COVER);
}

bool TagInfo::hasNothing() {
  return code == 0;
}

bool TagInfo::hasUndeterminedError() {
  return (code & MASK_UNDETERMINED_ERROR);
}

/* getNb* commands */

byte TagInfo::getNbBioUndetermined() {
  return getNbGenericElement(INFO_BIOLOGICAL_UNDETERMINED);
}

byte TagInfo::getNbRepetition() {
  return getNbGenericElement(INFO_REPETITION);
}

byte TagInfo::getNbSNP() {
  return getNbGenericElement(INFO_SNP);
}

byte TagInfo::getNbBioTagIndel() {
  return getNbGenericElement(INFO_BIOLOGICAL_TAG_INDEL);
}

uint TagInfo::getNbCauses() {
  return nb_causes;
}

byte TagInfo::getNbCausesInBreak(uint i, int strand) {
  if (strand == 1)
    return nb_causes_per_break[i];
  else
    return nb_causes_per_break[getNbBreaks() - i -1];
}

byte TagInfo::getNbSeqErr() {
  return getNbGenericElement(INFO_SEQ_ERR);
}

byte TagInfo::getNbSplice() {
  return getNbGenericElement(INFO_SPLICE);
}

byte TagInfo::getNbSpliceIntra() {
  return getNbGenericElement(INFO_INTRASPLICE);
}

byte TagInfo::getNbSpliceInter() {
  return getNbGenericElement(INFO_INTERSPLICE);
}

byte TagInfo::getNbSpliceNoCover() {
  return getNbGenericElement(INFO_SPLICE_NO_COVER);
}

byte TagInfo::getNbUndeterminedError() {
  return getNbGenericElement(INFO_UNDETERMINED_ERROR);
}

/* getInfos* commands */

DuplicationInfo **TagInfo::getInfosDuplication() {
  return (DuplicationInfo **) getGenericInfos(INFO_DUPLICATION);
}

RepetitionInfo **TagInfo::getInfosRepetition() {
  return (RepetitionInfo **)getGenericInfos(INFO_REPETITION);
}

SNPInfo **TagInfo::getInfosSNP() {
  return (SNPInfo **)getGenericInfos(INFO_SNP);
}

BioIndelInfo **TagInfo::getInfosBioTagIndel() {
  return (BioIndelInfo **)getGenericInfos(INFO_BIOLOGICAL_TAG_INDEL);
}

BioUndeterminedInfo **TagInfo::getInfosBioUndetermined() {
  return (BioUndeterminedInfo **)getGenericInfos(INFO_BIOLOGICAL_UNDETERMINED);
}

SeqErrInfo **TagInfo::getInfosSeqErr() {
  return (SeqErrInfo **)getGenericInfos(INFO_SEQ_ERR);
}

SpliceInfo **TagInfo::getInfosSplice() {
  return (SpliceInfo **)getGenericInfos(INFO_SPLICE);
}

SpliceIntraInfo **TagInfo::getInfosSpliceIntra() {
  return (SpliceIntraInfo **)getGenericInfos(INFO_INTRASPLICE);
}

SpliceInterInfo **TagInfo::getInfosSpliceInter() {
  return (SpliceInterInfo **)getGenericInfos(INFO_INTERSPLICE);
}

SpliceNoCoverInfo **TagInfo::getInfosSpliceNoCover() {
  return (SpliceNoCoverInfo **)getGenericInfos(INFO_SPLICE_NO_COVER);
}

UndeterminedErrorInfo **TagInfo::getInfosUndeterminedError() {
  return (UndeterminedErrorInfo **)getGenericInfos(INFO_UNDETERMINED_ERROR);
}

/* add commands */

void TagInfo::addDuplication() {
  code |= MASK_DUPLICATION;
  //  addGenericElement(INFO_DUPLICATION, new DuplicationInfo());
#ifdef DBG
  cerr << "Adding Duplication for tag " << num << endl;
#endif
}

void TagInfo::addAlmostNormal() {
  code |= MASK_ALMOST_NORMAL;
#ifdef DBG
  cerr << "Adding almostNormal for tag " << num << endl;
#endif
}

void TagInfo::addMultiple() {
  code |= MASK_MULTIPLE;
#ifdef DBG
  cerr << "Adding Multiple for tag " << num << endl;
#endif
}


void TagInfo::addNormal() {
  code |= MASK_NORMAL;
#ifdef DBG
  cerr << "Adding Normal for tag " << num << endl;
#endif
}

void TagInfo::addNone() {
  code |= MASK_NONE;
#ifdef DBG
  cerr << "Adding None for tag " << num << endl;
#endif
}

void TagInfo::addRepetition(int position_in, int position_out) {
  code |= MASK_REPETITION;
  addGenericElement(INFO_REPETITION, new RepetitionInfo( position_in, position_out));
}

void TagInfo::addSingle() {
  code |= MASK_SINGLE;
}

void TagInfo::addSNP(char tag_nuc, error_context snpID ) {
  ChrPosition chrPos;
  int shift = 0;
  char genome_nuc = 0;

  pair<ChrPosition, uchar *> result = retrieveChrPosAndNuc(snpID, 1, shift);

  if (result.second) {
    genome_nuc = result.second[0];
    free(result.second);
  }
  chrPos = result.first;

  // if (genome_nuc == tag_nuc) {
  //   cerr << "Warning: Discarding SNP at position " 
  //        << getPositionEndBreak(current_break) - shift << " in read "
  //        << support->getTagNum() << endl;
  //   const char mymsg[] = "Should be a SNP. But doesn't make sense!";
  //   uint len = strlen(mymsg);
  //   char *msg = new char [len+1];
  //   strcpy(msg, mymsg);
  //   msg[len] = 0;
  //   addGenericElement(INFO_UNDETERMINED_ERROR, new UndeterminedErrorInfo(msg));
  //   return ;
  // }

  code |= MASK_SNP;
  addGenericElement(INFO_SNP, new SNPInfo(chrPos, 
					  getPositionEndBreak(current_break) - shift,
					  getScoreComputedIntraExon(current_break),
					  isDuplicated(current_break),
					  tag_nuc, genome_nuc));
}

void TagInfo::addSNP(uint pos, char base_read, char base_genome, 
                     ChrPosition &chrPos, float score) {
  code |= MASK_SNP;
  addGenericElement(INFO_SNP, 
                    new SNPInfo(chrPos, pos, score, isDuplicated(current_break),
				base_read, base_genome));
}

void TagInfo::addBioTagIndel(uint nb_ins, uint nb_del) {
  code |= MASK_BIOLOGICAL_TAG_INDEL;

  ChrPosition chrPos = *getChrPosStartBreak(current_break);

  if (chrPos.getStrand() == -1) {
    if (nb_ins > 0)
      chrPos -= 1;
    else 
      chrPos -= nb_del + 1;
  } else {
    if (nb_ins > 0)
      chrPos = *getChrPosEndBreak(current_break) - 1;
    else 
      chrPos = *getChrPosEndBreak(current_break) - nb_del + 1;
  }
  if (nb_ins > 0)
    addGenericElement(INFO_BIOLOGICAL_TAG_INDEL, 
		      new BioIndelInfo(chrPos,
				       getPositionEndBreak(current_break) -nb_ins + 1,
				       nb_ins, nb_del, getScoreComputedIntraExon(current_break),isDuplicated(current_break)));
  else
    addGenericElement(INFO_BIOLOGICAL_TAG_INDEL, 
		      new BioIndelInfo(chrPos,
				       getPositionEndBreak(current_break) + nb_del,
				       nb_ins, nb_del, getScoreComputedIntraExon(current_break),isDuplicated(current_break)));
  
}


void TagInfo::addBioUndetermined(ChrPosition *chrPos, uint position,
    BioUndeterminedInfo::Type type,
    const char *message, ...) {
  code |= MASK_BIOLOGICAL_UNDETERMINED;
  char *error = new char[MAX_SIZE_MESSAGE_UNDETERMINED_ERROR];
  va_list args;
  va_start(args, message);
  vsprintf(error, message, args);
  addGenericElement(INFO_BIOLOGICAL_UNDETERMINED, new BioUndeterminedInfo(chrPos, position, type, error,support->getParameters()->detailed_sam));
  va_end(args);
}

void TagInfo::addSeqErr(error_context errorID, uint position, uint nb_ins, 
			uint nb_del, float score, 
			uint genome_seq_length,
			char *tag_seq, uint tag_seq_length) {
  if (tag_seq == NULL) 
    tag_seq_length =~0;
  ChrPosition chrPos;
  uchar *genome_seq = NULL;
  int shift=0;

  if (errorID == UNKNOWN_START_MISSING)
    chrPos = *getChrPosEndBreak(current_break);
  else if (errorID == UNKNOWN_END_MISSING)
    chrPos = *getChrPosStartBreak(current_break);
  else {
    pair<ChrPosition, uchar *> result = retrieveChrPosAndNuc(errorID, 
                                                             genome_seq_length,
                                                             shift);
    genome_seq = result.second;
    
    // if (genome_seq) {
    //   if (tag_seq && ! strcmp((char *)result.second, tag_seq)) {
    //     cerr << "Warning: Discarding error at position " << position << " in read "
    //          << support->getTagNum() << endl;
    //     delete [] tag_seq;
    //     free(result.second);
    //     const char mymsg[] = "Should be an error. But doesn't make sense!";
    //     uint len = strlen(mymsg);
    //     char *msg = new char [len+1];
    //     strcpy(msg, mymsg);
    //     msg[len] = 0;
    //     addGenericElement(INFO_UNDETERMINED_ERROR, new UndeterminedErrorInfo(msg));

    //     return ;
    //   }
    // }
    chrPos = ChrPosition(result.first);
  }
  if (genome_seq_length == 0)
    genome_seq_length = ~0;

  code |= MASK_SEQ_ERR;
  addGenericElement(INFO_SEQ_ERR, new SeqErrInfo(chrPos,
                                                 position, nb_ins, nb_del,
                                                 score, isDuplicated(current_break),
                                                 (char *)genome_seq, tag_seq,
                                                 genome_seq_length, tag_seq_length));
}


void TagInfo::addSplice(uint gap_length) {
  code |= MASK_SPLICE;

  ChrPosition pos = *getChrPosStartBreak(current_break);
  if (pos.getStrand() == 1) {
    pos += getThreshold() - 1;
  } else {
    pos = *getChrPosEndBreak(current_break) + getThreshold() - 1;
  }

  addGenericElement(INFO_SPLICE,
		    new SpliceInfo(pos, getPositionEndBreak(current_break), gap_length,isDuplicated(current_break)));
}

void TagInfo::addSpliceNoCover(uint gap_length) {
  code |= MASK_SPLICE_NO_COVER;

  ChrPosition pos = *getChrPosStartBreak(current_break);
  if (pos.getStrand() == 1) {
    pos += getThreshold() - 1;
  } else {
    pos = *getChrPosEndBreak(current_break) + getThreshold() - 1;
  }

  addGenericElement(INFO_SPLICE_NO_COVER,
		    new SpliceNoCoverInfo(pos, getPositionEndBreak(current_break), 
                                          gap_length,isDuplicated(current_break)));
}

void TagInfo::addSpliceIntra(uint gap_length) {
  code |= MASK_INTRA_TRANSPLICING;

  ChrPosition pos = *getChrPosStartBreak(current_break);
  if (pos.getStrand() == 1) {
    pos += getThreshold() - 1;
  } else {
    pos = *getChrPosEndBreak(current_break) + getThreshold() - 1;
  }

  addGenericElement(INFO_INTRASPLICE,
		    new SpliceIntraInfo(pos, getPositionEndBreak(current_break), 
                                        gap_length,isDuplicated(current_break)));
}

void TagInfo::addSpliceInter(float chimera_score,const string chim_score_info) {
  code |= MASK_INTER_TRANSPLICING;

  ChrPosition chrPos1 = *getChrPosStartBreak(current_break), chrPos2;
  if (chrPos1.getStrand() == 1) {
    chrPos1 += getThreshold() - 1;
  }
  chrPos2 = *getChrPosEndBreak(current_break);
  
  if (chrPos2.getStrand() == -1){
    chrPos2 += getThreshold() + 1;
  }

  // only for the data simulated otherwise there are many FP
  // the reason is because all chimera are built on the strand +
  // so the reading is inversed if strand(chrPos1) == -1
  //   if (chrPos1.getStrand() == -1) {
  //     ChrPosition tmp = chrPos1;
  //     chrPos1 = chrPos2;
  //     chrPos2 = tmp;
  //   }
  
  addGenericElement(INFO_INTERSPLICE, 
		    new SpliceInterInfo(chrPos1, getPositionEndBreak(current_break), 
					chrPos2, isDuplicated(current_break), 
					support->getBreak(current_break)->getChimeraClass(),
					chimera_score, chim_score_info));
}


void TagInfo::addUndeterminedError(UndeterminedErrorInfo::Type type,const char *format, ...) {
  code |= MASK_UNDETERMINED_ERROR;
  char *error = new char[MAX_SIZE_MESSAGE_UNDETERMINED_ERROR];
  va_list args;
  va_start(args, format);
  vsprintf(error, format, args);
  addGenericElement(INFO_UNDETERMINED_ERROR, new UndeterminedErrorInfo(type,error,support->getParameters()->detailed_sam));
  va_end(args);
}

void TagInfo::removeSpliceInter(uint i,BioUndeterminedInfo::Type type, const char* info) {
  if (i >= 0 && i < getNbSpliceInter()) {
    SpliceInterInfo *spliceInter = getInfosSpliceInter()[i];
    ChrPosition *chrPos = new ChrPosition(spliceInter->getChrPosition());
    uint position = spliceInter->getPosition();

    char *error = new char[MAX_SIZE_MESSAGE_UNDETERMINED_ERROR];
    strcpy(error,info);
    
    changeGenericElement(INFO_INTERSPLICE,INFO_BIOLOGICAL_UNDETERMINED,
                         MASK_INTER_TRANSPLICING, MASK_BIOLOGICAL_UNDETERMINED,
                         i,new BioUndeterminedInfo(chrPos, position, type, error,
                         support->getParameters()->detailed_sam));
  }
}

SamLine *TagInfo::getSamLine(ChrPosition *loc, ulong pos_of_original_loc_in_read) {
  SamLine *sam = new SamLine();

  // Set Read Group ID
  sam->addOptionalField("RG","1");

  if(!loc) {
    loc = getLocation();
    pos_of_original_loc_in_read = getSupportObject()->getPositionOfLocation();
  }

  //uint flag = 0;
  ChrPosition *loc_to_display = loc ; // Location to be displayed

  bool unmapped = ! loc || isNone(); // || (isMultiple() && !support->getParameters()->treat_multiple); 

  int strand = 1;
  int bases5_not_in_chromosome = 0; // Bases that are mapped 5' but not
                                    // actually in the chromsome (must ignore
                                    // the 5' bases)
  ostringstream string_stream;
  string qname;
  string rname;
  string sequence = read->seq;

  if (!read->name) {
    qname = intToString(read->id);
  } else {
    qname = read->name;
  }
  sam->setQname(qname);
  sam->setSeq(sequence);

  if(read->qual) {
    string quality = read->qual;
    sam->setQual(quality);
  }

  if (unmapped) {
    sam->setSegmentUnmapped();
  } else {
    int first_nice_break = -1;  // First nice break on strand +1 (5' end)
    int last_nice_break = -1;   // Last nice break on strand -1 (5' end)
    int cause_type;
    int num_cause = 0;

    // Look for the first nice break 
    for (int i = 0; i < (int)getNbBreaks(); i++) {
      cause_type = getCauseType(num_cause);
      num_cause += getNbCausesInBreak(i);

      // Is it a nice break? In that case, update the variables
      if ((getBreak(i)->isNiceBreak() || cause_type == INFO_INTERSPLICE)
          && cause_type != INFO_UNDETERMINED_ERROR
          && cause_type != INFO_BIOLOGICAL_UNDETERMINED) {
        if (first_nice_break == -1 && getBreak(i)->getStrandLocation(START_BREAK) == 1)
          first_nice_break = i;
        if (getBreak(i)->getStrandLocation(END_BREAK) == -1)
          last_nice_break = i;
      }
    }

    // Define the location of the first nice break
    if (first_nice_break > -1 || last_nice_break > -1) {
      uint pos5_nice;           // 5' position of the 5'-most nice break
      uint end5;               // last valid position starting from pos5_nice
                                // and going to the 5' end.
      if (first_nice_break > -1)
        strand = 1;
      else if (last_nice_break > -1)
        strand = -1;
      
      // Define some variables so that we don't need to take care of strand
      // afterwards
      pos_location_break type_loc = (strand == 1) ? START_BREAK : END_BREAK;
      int id_nice_break = (strand == 1) ? first_nice_break : last_nice_break;

      // We take the location of the first nice
      // break and we shift it to the left (5') as much as we can
      loc_to_display = getBreak(id_nice_break)->getChrPosition(type_loc);
      pos5_nice = getBreak(id_nice_break)->getPositionStartBreak(strand) - 1;
      end5 = (((strand == 1 && id_nice_break > 0) 
               || (strand == -1 && id_nice_break < (int)getNbBreaks() - 1))
              ? getBreak(id_nice_break - strand)->getPositionEndBreak(strand) + 1
              : 0);
      if (end5 > pos5_nice)
        end5 = pos5_nice;

      if (loc_to_display->getRelativePosition() < (pos5_nice - end5)) {
        // Oops we are before the chromosome!
        bases5_not_in_chromosome = (pos5_nice - end5) - loc_to_display->getRelativePosition();
        // Go back a little!
        end5 += bases5_not_in_chromosome;
      }
      *loc_to_display -= (pos5_nice - end5);
    }

    if (loc_to_display->getStrand() == -1) {
      sam->setSeqReverseComplemented();
      //flag = 16;
      strand = -1;
    }
  }

  if (!unmapped) {
    // String for what is matched, what is inserted, what is deleted
    // M : alignment match (sequence match or mismatch
    // I : insertion
    // D : deletion
    // N : intron -> 0N special case for telling there is a chimera
    // = : seq. match
    // X : seq. mismatch
    // S : soft clipping (unaligned)
    byte current_cause_type[INFO_NB_TYPES];
    uint current_cause = 0;
    SeqErrInfo **infos = getInfosSeqErr();
    int nb_suspicious_breaks = 0;    /* number of breaks where the mapping
                                        may not be very accurate
                                     */
    for (uint i=0; i < INFO_NB_TYPES; i++)
      current_cause_type[i]= (strand == 1) ? 0 : getNbGenericElement(i)-1;

    bool ended = false;
    Cigar cigar(support->getParameters()->use_x_in_cigar);
    //Cigar cigar;

    bool started_for_real = false; /* put at true when we have something else
                                      than no start break or undetermined
                                   */
    uint last_match_pos = 0;      // Last time a match was seen.

    // Do we have at least one nice break?
    if (loc != loc_to_display) {

      // TODO: improve the test in order to update NH field
      if (strcmp(loc->getChrPosition(),loc_to_display->getChrPosition()) != 0 || loc->getStrand() != loc_to_display->getStrand())
	is_good_location = false;
      
      
      for (uint i=0; i < getNbBreaks(); i++) {
        SupportBreak *sb = getBreak(i, strand);

        for (uint j = 0; j < getNbCausesInBreak(i, strand); j++) {
          byte cause = getCauseType(current_cause, strand);
          bool is_undetermined = 
            cause == INFO_BIOLOGICAL_UNDETERMINED
            || cause == INFO_UNDETERMINED_ERROR
            || (cause == INFO_SEQ_ERR 
                && infos[current_cause_type[cause]]->getNbIns() == (uint) ~0
                && infos[current_cause_type[cause]]->getNbDel() == (uint) ~0);

          if (getPositionStartBreak(i, strand) == 0 
              || (is_undetermined && !started_for_real)) {
            last_match_pos = getPositionEndBreak(i, strand)+1;
          } else {
            // First time we arrive here.
            // Write down the S we may have
            if (! started_for_real && last_match_pos > 0) {
              cigar.append(last_match_pos, 'S');
            }

            started_for_real = true;

            if (j > 0 
                && ((cause == INFO_SNP 
                     && getCauseType(current_cause-1, strand) == INFO_SNP)
                    || (cause == INFO_SEQ_ERR 
                        && getCauseType(current_cause-1, strand) == INFO_SEQ_ERR))) {
              // Double SNP or double error
              // We need to take into account the bases that are matched
              // between the two mutations
              cigar.append(getPositionEndBreak(i, strand)
                        - (getPositionStartBreak(i, strand) 
                           + sb->getThreshold()),
                        'M');
            } else {
              // Write down the numbers of match we have.

              if (is_undetermined 
                  // either a too large insertion in the read
                  // that is not compatible with the gap length
                  // -> doesn't make sense.
                  && ((sb->getNbGenomeIndels() < 0
                       && - sb->getNbGenomeIndels() 
                       > (int) (getPositionEndBreak(i, strand)
                                - getPositionStartBreak(i, strand)))
                      // or too large genome gap that cannot be considered
                      // as a splice
                      || (sb->getNbGenomeIndels() > 0
                          && sb->getNbGenomeIndels() > (int)getSupportObject()->getParameters()->max_splice_length)
                      // or we are dealing with a kind of chimera
                      || sb->isChimera())){
                // Since these cases are not coherent, 
                // we tell to ignore the read until the end.
                cigar.append(getPositionStartBreak(i, strand) 
                          + sb->getThreshold() -2
                          - last_match_pos + 1, 'M');
                cigar.append(read->getLength()
                          - (getPositionStartBreak(i, strand) 
                             + sb->getThreshold() -1), 
                          'S');
                ended = true;
                // Leave the two outer loops
                i = getNbBreaks();
                break;
              }
              // Large insertion in read  are not like usual insertions
              // they should be considered like chimera and therefore
              // like large... deletions!
              // same for insertions that lead to a negative gap length on the
              // genome
              bool large_insertion = sb->getNbGenomeIndels() <= 0 
                && (((- sb->getNbGenomeIndels() 
                     >= (int)(getPositionEndBreak(i, strand) 
                         - getPositionStartBreak(i, strand)+1))
                    || (sb->getGenomeGapLength() <= 0))
                    || (sb->isChimera() && sb->isBiologicalInterEvent()));

              uint min_last_pos = 
                min(getPositionStartBreak(i, strand) + sb->getThreshold() -2,
                    sb->hasNoEndBreak(strand)
                    ? (uint)read->getLength() - 1 
                    : getPositionEndBreak(i, strand) 
                    - (sb->getNbGenomeIndels() < 0 && ! large_insertion
                       ? - sb->getNbGenomeIndels() : 1));

              uint nb_match = min_last_pos      
                // Beware short gaps!
                - last_match_pos + 1;
              // In case of deletion, we may have one more match, in fact
              if ((sb->getNbGenomeIndels() > 0 || large_insertion)
                  && min_last_pos < (getPositionStartBreak(i, strand) 
                                     + sb->getThreshold()
                                     -2))
                nb_match++;
              cigar.append(nb_match, 'M');
            }


            // Take care of gaps that can be larger than expected
            // Of course we don't consider large gaps that contain several causes
            // nor end breaks
            if (getNbCausesInBreak(i, strand) == 1 
                && ! sb->hasNoEndBreak(strand)) {
              uint expected_max_gap_length = sb->getThreshold();
              if (sb->getNbGenomeIndels() > 0)
                expected_max_gap_length--;
              else if (sb->getNbGenomeIndels() < 0) 
                expected_max_gap_length += -sb->getNbGenomeIndels() - 1;
              if (expected_max_gap_length < sb->getTagBreakLength()) {
                if (is_undetermined)
                  nb_suspicious_breaks++;
                cigar.append(sb->getTagBreakLength() - expected_max_gap_length, 
                          'X');
              }
            }

            if (! is_undetermined) {
              // SNP, splice, indel, seq err (sub + indel), chimera
              cigar.append(getGenericInfos(cause)[current_cause_type[cause]]->getCigarInfo());
            } else {
              // Bio undetermined, undetermined, seq err

              if (sb->hasNoEndBreak(strand)) {
                // end of read
                ended = true;
                cigar.append(read->getLength() - 1
                          - (getPositionStartBreak(i, strand) 
                             + sb->getThreshold()-1) + 1, 'S');
              } else {
                if (sb->getNbGenomeIndels() > 0) {
                  // We have a kind of splice
                  cigar.append(sb->getNbGenomeIndels(), 'N');
                } else if (sb->getNbGenomeIndels() < 0 
                         && sb->getGenomeGapLength() > 0)
                  cigar.append(-sb->getNbGenomeIndels(), 'I');
                else {
                  // Well, well, well, ... don't know what to say about
                  // that gap. Unfortunately there is no "unknown" operator
                  // in CIGAR. We use X instead...
                  cigar.append(1, 'X');
                }
              }
            }
            if (strand == -1)
              getGenericInfos(cause)[current_cause_type[cause]]
                ->setStrandedPosition(sb->getPositionEndBreak(strand));
          }

          current_cause_type[cause] += strand;
          current_cause++;
        }
        if (! ended)
          last_match_pos = getPositionEndBreak(i, strand) + 1;
      }

      // And treat all the remaining bases that may have mapped
      // at the end of the read after the breaks (but be careful to ignored
      // breaks)
      if (! ended) {
        cigar.append(read->getLength() - last_match_pos,
                  'M');
      }
    } else {
      // If we only had undetermined breaks, we have no reliable way
      // of identifying the read's location.
      // Therefore, we are falling back to the location found 
      // in Support (chosen in the longest run of single locations).
      // We will set the cigar so that it corresponds to this fact.

      // In that case loc == loc_to_display. Separate them, so that by
      // modifying one, we don't modify the other
      loc_to_display = new ChrPosition(*loc);

      // Check how far we can extend that to the left and to the right
      // (BTW choose the leftmost position)
      uint pos_shift[2];
      // Index of the left (5') end and the right (3') end
      int left_idx = (strand == 1) ? 0 : 1;
      int right_idx = (strand == 1) ? 1: 0;

      pos_shift[0]= pos_of_original_loc_in_read;
      pos_shift[1] = pos_of_original_loc_in_read;

      // If we have no break the leftmost and rightmost positions
      // are quite easy to find.
      if (getNbBreaks() == 0) {
        pos_shift[0] = 0;
        pos_shift[1] = getSupportLength() - 1;
      } else {
        while (pos_shift[left_idx] > 0 
               && pos_shift[left_idx] < getSupportLength() - 1
               && getLocalisations()[pos_shift[left_idx]] > 0) 
          pos_shift[left_idx] -= strand;
        if (getLocalisations()[pos_shift[left_idx]] == 0)
          pos_shift[left_idx] += strand;

        while (pos_shift[right_idx] > 0 
               && pos_shift[right_idx] < getSupportLength() - 1
               && getLocalisations()[pos_shift[right_idx]] > 0) 
          pos_shift[right_idx] += strand;
        if (getLocalisations()[pos_shift[right_idx]] == 0)
          pos_shift[right_idx] -= strand;
      }
    
      // Take the leftmost (5') position
      if (loc_to_display->getRelativePosition() < strand * (pos_of_original_loc_in_read 
                                                            - pos_shift[left_idx])) {
        // Oops we are before the start of the chromosome !
        pos_shift[left_idx] = pos_of_original_loc_in_read - loc_to_display->getRelativePosition() * strand;
      }
      *loc_to_display -= strand * (pos_of_original_loc_in_read 
                                   - pos_shift[left_idx]);
      pos_of_original_loc_in_read = pos_shift[left_idx];

      int bases[2] = {0,0};
      bases[0] = pos_shift[0];
      if (read->getLength() -1
          > (pos_shift[1] + getThreshold() - 1))
        bases[1] = read->getLength() -1
          - (pos_shift[1] + getThreshold() - 1);
      
      // set CIGAR
      if (bases[left_idx])
        cigar.append(bases[left_idx], 'S');
      // The match is only considered between the leftmost and rightmost
      // positions
      cigar.append((read->getLength() 
                 - bases[left_idx]
                 - bases[right_idx]), 'M');
      if (bases[right_idx]) 
        cigar.append(bases[right_idx], 'S');

    }

    // If we had bases to ignore in 5'
    if (bases5_not_in_chromosome) {
      if (cigar[0].type == 'S') {
        assert(cigar.count() >= 2);
        cigar[0].nb += bases5_not_in_chromosome;
        cigar[1].nb -= bases5_not_in_chromosome;
      } else {
        // We have a match from the beginning. We need to add a S at the start!
        cigar[0].nb -= bases5_not_in_chromosome;
        // Add the S
        cigar.prepend(bases5_not_in_chromosome, 'S');
      }
    }

    // Join duplicates, remove empty operators, ...
    cigar.filter();

    //sam->setFlag(flag);
    rname = loc_to_display->getChrPosition();
    sam->setRname(rname,loc_to_display->getChrId());
    sam->setPos(loc_to_display->getRelativePosition() + 1);
    sam->setMapQ(unmapped? 0 : 254 - (getNbBreaks() > 0 
                 ? 204.*nb_suspicious_breaks/getNbBreaks()
                 : 0));
    sam->setCigar(cigar);

    // XS field for splicing event (cufflinks or stringTie)
    if (hasSplice() || hasSpliceNoCover() || hasSpliceIntraChr()){
      if (strand == 1)
	sam->addOptionalField("XS",'+');
      else
	sam->addOptionalField("XS",'-');
    }
    
    // Number of causes
    sam->addOptionalField("XC",getNbCauses());
    
    // Repetition, XX field
    if(hasRepetition()) {
      string_stream.str("");
      for (uint i = 0; i < nb_each_type[INFO_REPETITION]; i++) {
	if(string_stream.str() != "")
	  string_stream << ";";
	info_each_type[INFO_REPETITION][i]->samOutput(string_stream, strand);
      }
      sam->addOptionalField("XX",string_stream.str());
    }
  }

  // k-mer chosen
  // if (pos_of_loc_cause_in_read != (uint)~0) {
  //   os << "XL:Z:" << loc_cause << "\t";
  //   os << "XP:i:" << pos_of_loc_cause_in_read << "\t";
  // }
  if (loc) {
    string_stream.str("");
    string_stream << *loc;
    sam->addOptionalField("XO",string_stream.str());
    // XQ was reversed when strand==-1
    if (sam->isSeqReverseComplemented()) {
      sam->addOptionalField("XQ",getSupportLength() - 1 - getSupportObject()->getPositionOfLocation());
    }else{
      sam->addOptionalField("XQ",getSupportObject()->getPositionOfLocation());
    }
  }

  
  // reverse complemente sequence and quality
  if (sam->isSeqReverseComplemented()) {
    sam->reverseComplementeSeq();
    sam->reverseQual();
  }
  
  // Unique, duplicate, multiple
  sam->addOptionalField("XU",isSingle());
  sam->addOptionalField("XD",isDuplication());
  sam->addOptionalField("XM",isMultiple());
  uint xn_value = 0;
  if (isNormal()) {
    xn_value = 1;
  } else if (isAlmostNormal()) {
    xn_value = 2;
  }
  sam->addOptionalField("XN",xn_value);

  
  // Display causes in optional field XE if there is any
  // Events are reversed when strand == -1
  if(getNbCauses() > 0 && !unmapped) {
    byte *current_each_type = (byte *)calloc(INFO_NB_TYPES, sizeof(byte));
    uint this_break_id = 0;
    uint current_nb_causes_this_break = 0;
    int strand_info_type = strand;
    string_stream.str("");
    if (sam->isSeqReverseComplemented()) {
      for (int i = getNbCauses() -1 ; i >= 0; i--) {
	GenericInfo *info_type = info_each_type[types_in_order[i]][current_each_type[types_in_order[i]]];
	if (types_in_order[i] == INFO_SNP){
	  strand_info_type = ((SNPInfo *) info_type)->getChrPosition().getStrand();
	}else if (types_in_order[i] == INFO_BIOLOGICAL_TAG_INDEL){
	  strand_info_type = ((BioIndelInfo *) info_type)->getChrPosition().getStrand();
	}
	current_nb_causes_this_break++;
    	// add separator in case of multiple events
    	if(string_stream.str() != "")
    	  string_stream << ";";
    	string_stream << i << ":" << this_break_id << ":";
    	info_type->samOutput(string_stream, strand_info_type);
    	current_each_type[types_in_order[i]]++;
    	if (nb_causes_per_break[this_break_id] == current_nb_causes_this_break) {
    	  this_break_id++;
    	  current_nb_causes_this_break = 0;
    	}
      }
    }else{
      for (uint i = 0; i < getNbCauses(); i++) {
	GenericInfo *info_type = info_each_type[types_in_order[i]][current_each_type[types_in_order[i]]];
	if (types_in_order[i] == INFO_SNP){
	  strand_info_type = ((SNPInfo *) info_type)->getChrPosition().getStrand();
	}else if (types_in_order[i] == INFO_BIOLOGICAL_TAG_INDEL){
	  strand_info_type = ((BioIndelInfo *) info_type)->getChrPosition().getStrand();
	}
	current_nb_causes_this_break++;
	// add separator in case of multiple events
	if(string_stream.str() != "")
	  string_stream << ";";
	string_stream << i << ":" << this_break_id << ":";
	info_type->samOutput(string_stream, strand_info_type);
	current_each_type[types_in_order[i]]++;
	if (nb_causes_per_break[this_break_id] == current_nb_causes_this_break) {
	  this_break_id++;
	  current_nb_causes_this_break = 0;
	}
      }
    }
    sam->addOptionalField("XE",string_stream.str());
    free(current_each_type);
  }

  // Breaks are reversed when strand == -1
  if (getSupportObject()->getParameters()->detailed_sam) {
    if(getNbBreaks() > 0) {
      string_stream.str("");
      if (sam->isSeqReverseComplemented()) {
      	for (int i = getNbBreaks() -1 ; i >= 0 ; i--) {
      	  if(string_stream.rdbuf()->in_avail() > 0)
      	    string_stream << ";";
      	  string_stream << i 
      			<< ":is_duplicated=" << getBreak(i, strand)->isDuplicated() 
      			<< ";genome_indels=" << getBreak(i, strand)->getNbGenomeIndels()
      			<< ";score_intra=" << getBreak(i, strand)->getScoreComputedIntraExon()
      			<< ";score_inter=" << getBreak(i, strand)->getScoreComputedInterExon()
      			<< ";deviation=" << getBreak(i, strand)->getScoreInsideAverages()
      			<< ";falling_left=" << getBreak(i, strand)->isSupportFallingLeft(strand)
      			<< ";falling_right=" << getBreak(i, strand)->isSupportFallingRight(strand)
      			<< ";inside_first_quartile=" 
      			<< getBreak(i, strand)->getInsideQuartile1()
      			<< ";inside_last_quartile="
      			<< getBreak(i, strand)->getInsideQuartile4()
      			<< ";inside_score=" << getBreak(i, strand)->getInsideScore()
      			<< ";outside_score=" << getBreak(i, strand)->getOutsideScore()
      			<< ";average_low_inside=" << getBreak(i, strand)->getAverageLowInside()
      			<< ";average_high_inside=" << getBreak(i, strand)->getAverageHighInside()
      			<< ";has_no_start_break=" << getBreak(i, strand)->hasNoStartBreak(strand)
      			<< ";has_no_end_break=" << getBreak(i, strand)->hasNoEndBreak(strand)
      			<< ";is_deviated=" << getBreak(i, strand)->isDeviated()
      			<< ";is_nice_break=" << getBreak(i, strand)->isNiceBreak()
      			<< ";is_very_nice_break=" << getBreak(i, strand)->isVeryNiceBreak()
      			<< ";pos_start_break=" << getPositionStartBreak(i, strand)
      			<< ";pos_end_break=" << getPositionEndBreak(i, strand);
      	}
      }else{
	for (uint i = 0; i < getNbBreaks(); i++) {
	 if(string_stream.rdbuf()->in_avail() > 0)
	   string_stream << ";";
	 string_stream << i 
		       << ":is_duplicated=" << getBreak(i, strand)->isDuplicated() 
		       << ";genome_indels=" << getBreak(i, strand)->getNbGenomeIndels()
		       << ";score_intra=" << getBreak(i, strand)->getScoreComputedIntraExon()
		       << ";score_inter=" << getBreak(i, strand)->getScoreComputedInterExon()
		       << ";deviation=" << getBreak(i, strand)->getScoreInsideAverages()
		       << ";falling_left=" << getBreak(i, strand)->isSupportFallingLeft(strand)
		       << ";falling_right=" << getBreak(i, strand)->isSupportFallingRight(strand)
		       << ";inside_first_quartile=" 
		       << getBreak(i, strand)->getInsideQuartile1()
		       << ";inside_last_quartile="
		       << getBreak(i, strand)->getInsideQuartile4()
		       << ";inside_score=" << getBreak(i, strand)->getInsideScore()
		       << ";outside_score=" << getBreak(i, strand)->getOutsideScore()
		       << ";average_low_inside=" << getBreak(i, strand)->getAverageLowInside()
		       << ";average_high_inside=" << getBreak(i, strand)->getAverageHighInside()
		       << ";has_no_start_break=" << getBreak(i, strand)->hasNoStartBreak(strand)
		       << ";has_no_end_break=" << getBreak(i, strand)->hasNoEndBreak(strand)
		       << ";is_deviated=" << getBreak(i, strand)->isDeviated()
		       << ";is_nice_break=" << getBreak(i, strand)->isNiceBreak()
		       << ";is_very_nice_break=" << getBreak(i, strand)->isVeryNiceBreak()
		       << ";pos_start_break=" << getPositionStartBreak(i, strand)
		       << ";pos_end_break=" << getPositionEndBreak(i, strand); 
	}
      }
      sam->addOptionalField("XB",string_stream.str());
    }


    // We used to reverse p_support and p_loc, but there was a time where
    // it was a problem (bug #16502, see J. Audoux the asshole), so we drew back and do not reverse them.
    // Later, it was an other problem (bug #18524) because pos_junction was always reversed but not the p_loc anymore...
    // In fact, we have to make a choice: either we must reverse all fieds in order to be "sam consistent"
    // or none. We choose the first one (ie reverse all) and we update unit tests.
    // p_loc and p_support are reversed when strand == -1

    string_stream.str("");
    string_stream << "p_support=";
    if (sam->isSeqReverseComplemented()) {
      for (uint j=getSupportLength() -1; j > 0; j--) {
        string_stream << getSupport()[j] << "," ;
      }
      string_stream << getSupport()[0] << ";p_loc=";
      for (uint j=getSupportLength() -1; j > 0; j--) {
        string_stream <<  getLocalisations()[j] << "," ;
      }
      string_stream << getLocalisations()[0];
    } else {
      for (uint j=0; j < getSupportLength()-1; j++) {
        string_stream << getSupport()[j] << "," ;
      }
      string_stream << getSupport()[getSupportLength()-1] << ";p_loc=";
      for (uint j=0; j < getSupportLength()-1; j++) {
        string_stream <<  getLocalisations()[j] << "," ;
      }
      string_stream << getLocalisations()[getSupportLength()-1];
    }
    sam->addOptionalField("XR",string_stream.str());
  }
  
  if (loc_to_display != loc)
    delete loc_to_display;

  return sam;
}

vector<SamLine*> *TagInfo::getChimericAlignements(SamLine *base_line, uint cigar_index, uint splice_id, uint left_softclip, uint *right_softclip) {
  if(cigar_index >= base_line->getCigar().count()) {
    // If there is no more 0N (chimeras)
    return new vector<SamLine*>();
  } else {
    // Clone base line
    SamLine *chimeric_alignement = new SamLine(*base_line);
    // Create a new cigar for this slice
    Cigar new_cigar(support->getParameters()->use_x_in_cigar);
    //Cigar new_cigar;
    
    // Create the cigar for this slice util we found an other chimera
    while(cigar_index < base_line->getCigar().count() && getInfosSpliceInter()[0]->getCigarInfo() != base_line->getCigar().get(cigar_index)) {
      new_cigar.append(base_line->getCigar().get(cigar_index));
      cigar_index++;
    }
    
    uint nb_positions_read_before_softclipping = new_cigar.getNbPositionsRead();

    // RECURSION
    vector<SamLine*> *sam_lines = getChimericAlignements(base_line,cigar_index+1,splice_id+1,left_softclip+nb_positions_read_before_softclipping,right_softclip);

    // Add right softclip calculated from previous child
    if(*right_softclip > 0) {
      new_cigar.append(*right_softclip,'S');
    }
    // Update right softclip with length of this slice
    *right_softclip+=nb_positions_read_before_softclipping;


    // If the strand has change we need to reverse seq and qual to always
    // be on strand forward only when this is not the primary alignement
    if(splice_id > 0) {
      chimeric_alignement->setChimericAlignement();
      // If this is not the primary line
      // Update chromosome and start position
      // with chimera info
      string new_rname = getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getChrPosition();
      int new_rid = getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getChrId();
      chimeric_alignement->setRname(new_rname,new_rid);
      if(getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getStrand() == -1) {
        uint pos = getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getRelativePosition();
        // We start from this end pos and we go back until we reach the last Cigar element
        for(uint i=0; i < new_cigar.count(); i++) {
          if(new_cigar.isReferenceBasedElement(i)) {
            pos -= new_cigar[i].nb;
          }
        }
        chimeric_alignement->setPos(pos);
      } else {
        chimeric_alignement->setPos(getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getRelativePosition());
      }
      new_cigar.prepend(left_softclip,'S');
      if(getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getStrand() == -1 && !chimeric_alignement->isSeqReverseComplemented()) {
        chimeric_alignement->setSeqReverseComplemented();
        chimeric_alignement->reverseComplementeSeq();
        chimeric_alignement->reverseQual();
        new_cigar.reverse();
      } else if(getInfosSpliceInter()[splice_id-1]->getChromosomeDest().getStrand() == 1 && chimeric_alignement->isSeqReverseComplemented()) {
        chimeric_alignement->unsetSeqReverseComplemented();
        chimeric_alignement->reverseComplementeSeq();
        chimeric_alignement->reverseQual();
        new_cigar.reverse();
      }
      // We remove CRAC extend fields
      chimeric_alignement->removeAllUserOptionalFields();
    }
    //
    // set cigar
    chimeric_alignement->setCigar(new_cigar);

    // Add the line to the vector
    sam_lines->insert(sam_lines->begin(),chimeric_alignement);
    return sam_lines;
  }
}


vector<SamLine*> *TagInfo::getSamLines() {
  //sam_alignements sam_alignements;
  vector<SamLine*> *sam_lines;
  // Multiple alignements
  // we consider several alignments only if there is no break (otherwise it is a little bit tricky)
  if ((
	    (support->getParameters()->treat_multiple > 0 && (support->isMultiple() || support->isDuplicate())) 
	    //|| (support->isDuplicate() && !support->isMultiple())
	    )
	   && support->isContinuous()
     && !support->isSingle()) {
    sam_lines = new vector<SamLine*>;
    uint pos_loc = support->getPositionOfLocation();
    ChrPosition **locs = support->getLocations();
    uint nb_locs = min(support->getNbLocsMax(), support->getParameters()->treat_multiple);
    for(ulong i = 0; i < nb_locs ; i++) {
      SamLine *multiple_alignement = getSamLine(locs[i], pos_loc);
      // If this is not the primary line
      // We set the proper flag and remove the sequence and qual from the
      // sam_line to reduce file size (as it recommended in SAM specifications)
      if(i>0) {
        multiple_alignement->setSecondaryAlignement();
        multiple_alignement->setSeq("*");
        multiple_alignement->setQual("*");
      }
      sam_lines->push_back(multiple_alignement);
      //      delete locs[i];
    }        
    // if (locs!=NULL)
    //   free(locs);

    // Set flags for "Next hit" : CC, CP
    vector<SamLine*>::iterator it = sam_lines->begin() ;
    vector<SamLine*>::iterator next_it = sam_lines->begin() ;
    next_it++;
    for (; next_it != sam_lines->end(); ++next_it) {
      if((*next_it)->getRname() == (*it)->getRname()) {
        (*it)->addOptionalField("CC","=");
      } else {
        (*it)->addOptionalField("CC",(*next_it)->getRname());
      }
      (*it)->addOptionalField("CP",(*next_it)->getPos());
      it++;
    }
  }
  // Hack for chimeric alignement
  // Since SAM format now supports chimeric alignement
  // we need to provide one SAM line for each alignement
  else {
    SamLine *base_line = getSamLine();
    if(base_line->getCigar().getNbChimeras() > 0) {
      uint right_softclip = 0;
      sam_lines = getChimericAlignements(base_line, 0, 0, 0, &right_softclip);
      delete base_line;

      // Generating SA optional field
      ostringstream string_stream;
      for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
        string_stream.str(""); // reset buffer
        for (vector<SamLine*>::iterator it2 = sam_lines->begin() ; it2 != sam_lines->end(); ++it2) {
          // We do not print the current line info to SA field
          // that would be redoundant
          if(it != it2) {
            // Add a semicolon separator if this is not the first alignement
            if (string_stream.str().size() > 0) {
              string_stream << ";";
            }
            // print canonical alignemnt with coma separated fields
            string_stream << (*it2)->getRname() << ',' << (*it2)->getPos() << ',';
            if((*it2)->isSeqReverseComplemented()) {
              string_stream << "-";
            } else {
              string_stream << "+";
            }
            string_stream << ',' << (*it2)->getCigar() << ',' <<(*it2)->getMapQ() << ',' << (*it2)->getCigar().getEditDistance();
          }
        }
        (*it)->addOptionalField("SA",string_stream.str());
      }
    } else {
      sam_lines = new vector<SamLine*>;
      sam_lines->push_back(base_line);
    }
  }


  // HACK: we check that the given alignements are not going out of scope (ie.
  // longer than the chromosome length).
  // If it does we softclip the end of the alignement
  for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
    SamLine *sam_line = *it;
    // We do nothing if the read is not mapped
    if(!sam_line->isSegmentUnmapped()) {
      uint chr_length = genome->getGenomeInfo()->getChrLength(sam_line->getRid());
      uint start = sam_line->getPos();
      Cigar cigar = sam_line->getCigar();
      // We first check if the read is going out of scope before we
      // do some time processing operations
      if(start + cigar.getReferenceAlignementLength() > chr_length) {
        Cigar new_cigar(support->getParameters()->use_x_in_cigar);
        bool find_end = false;
        for(uint i=0; i < cigar.count(); i++) {
          // If we have reach the end of the chromose and the cigar operator correspond
          // to the read sequence we place a softclip or hardclip if it was the case
          if(find_end) {
            if(cigar.isQueryBasedElement(i)) {
              if(cigar[i].type == 'H') {
                new_cigar.append(cigar[i].nb, 'H');
              } else {
                new_cigar.append(cigar[i].nb, 'S');
              }
            }
          } else {
            if(cigar.isReferenceBasedElement(i) && start + cigar[i].nb > chr_length) {
              if(cigar.isQueryBasedElement(i)) {
                // Only append if the value is positive
                if(chr_length > start)
                  new_cigar.append(chr_length - start + 1,cigar[i].type);
                // Append a soft-clip for the remaining length
                new_cigar.append(start + cigar[i].nb - chr_length - 1, 'S');
              }
              find_end = true;
            } else {
              new_cigar.append(cigar[i]);
            }
          }
        }
        // Join duplicates, remove empty operators, ...
        new_cigar.filter();
        // We set the newly created cigar
        sam_line->setCigar(new_cigar);
      }
      // Then we check that the alignment is not longer than the read itself
      if(sam_line->getSeq().compare("*") != 0 && sam_line->getSeq().length() != sam_line->getCigar().getQueryAlignementLength()) {
        // Display an error message to report a new bug
        cerr << "TagInfo error: CIGAR and sequence length are inconsistent, read is considered unmapped." << endl;
        cerr << "Read " << sam_line->getSeq() << 
          ", sequence length " << sam_line->getSeq().length() << 
          " vs " << sam_line->getCigar().getQueryAlignementLength() << " from CIGAR" << endl;
        cerr << "Please submit a new bug report to crac-bugs@lists.gforge.inria.fr" << endl;
        sam_line->unsetEachSegmentsMapped();
        sam_line->setSegmentUnmapped();
        sam_line->setUnknownPos();
        sam_line->setUnknownRname();
        sam_line->removeAllUserOptionalFields();
      }
    }
  }
    

  // Set optional fields
  uint i = 1;
  for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
    (*it)->addOptionalField("NM",(*it)->getCigar().getEditDistance());
    //(*it)->addOptionalField("IH",sam_lines->size());
    (*it)->addOptionalField("HI",i);
    if (isNone()) {
      (*it)->addOptionalField("NH",0);
    } else if (is_good_location) {
       (*it)->addOptionalField("NH",support->getNbLocsMax());
    } else {
      // TODO
      // Since we have change the location k-mer from a chromosome (or strand)
      // it means that we have at least 2 different hits. But there is maybe more...
      // How do we choose a nice value for NH field?
      (*it)->addOptionalField("NH",2);
      (*it)->setMapQ(0);
    }
    i++;
  }

  return sam_lines;
}

ostream &TagInfo::samOutput(ostream &os) {
  vector<SamLine*> *sam_lines = getSamLines();
  for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
    (*it)->writeLine(os);
    delete (*it);
  }
  delete sam_lines;
  return os;
}

void TagInfo::setCurrentBreak(uint i) {
  current_break = i;
}

ostream &operator<<(ostream &os, TagInfo &info) {
  //   uint current_break = info.getCurrentBreak();
  //   if (info.getNbBreaks() > 0 && current_break < info.getNbBreaks()) {
  //     os << "gap_indel=" << (int) (info.getTagBreakLength(current_break) 
  //                                  - info.getThreshold())
  //        << " score_intra=" << info.getScoreComputedIntraExon(current_break)
  //        << " score_inter=" << info.getScoreComputedInterExon(current_break)
  //        << " deviation=" << info.getScoreInsideAverages(current_break)
  //        << " falling_left=" << info.isSupportFallingLeft(current_break)
  //        << " falling_right=" << info.isSupportFallingRight(current_break)
  //        << " inside_first_quartile=" 
  //        << info.getBreak(current_break)->getInsideQuartile1()
  //        << " inside_last_quartile="
  //        << info.getBreak(current_break)->getInsideQuartile4() 
  //        << " " << endl;
  //   }

  if (info.getLocation() != NULL){
    os << *info.getLocation() << " ";
    os << "pos_location=" << info.getPositionOfLocation() << " ";
  }

  os << info.getRead()->seq << " " ;
  
  for (uint j=0; j < info.getSupportLength()-1; j++) {
    os << info.getSupport()[j] << "," ;
  }
  os << info.getSupport()[info.getSupportLength()-1] << " ";
  for (uint j=0; j < info.getSupportLength()-1; j++) {
    os <<  info.getLocalisations()[j] << "," ;
  }
  os << info.getLocalisations()[info.getSupportLength()-1];
  return os;
}


// PRIVATE
pair<ChrPosition, uchar *> TagInfo::retrieveChrPosAndNuc(error_context type,
                                                         uint nb_nuc,
                                                         int &shift) {
  uint absolute_genome_loc;
  ChrPosition chrPos = *getChrPosStartBreak(current_break);
  uchar *dna = NULL;
  if (type == SECOND_SUBSTITUTION) {
    shift = (int)  getTagBreakLength(current_break) - getThreshold();
  }
  if (chrPos.getStrand() == -1) {
    if (type != FIRST_SUBSTITUTION || type == INSERTION
        || type == DELETION) {
      chrPos += -1;
      absolute_genome_loc = getLocationStartBreak(current_break) - 1;
    } else {
      chrPos = *getChrPosEndBreak(current_break) + getThreshold();
      absolute_genome_loc = getLocationEndBreak(current_break) + getThreshold();
    }
  } else {
    chrPos = *getChrPosEndBreak(current_break) - shift - 1;
    absolute_genome_loc = getLocationEndBreak(current_break) - shift - 1;
  }
  if (type == DELETION)
    shift = -1;
  if (type != INSERTION && nb_nuc > 0) {
    if (nb_nuc <= support->getParameters()->max_bases_retrieved)
      dna = genome->getGenomeSubstring(absolute_genome_loc,
                                       nb_nuc, chrPos.getStrand());
    else
      cerr << "Warning substring of length " 
           << nb_nuc << " not retrieved from genome (read " 
           << read->id << ")" << endl;
  }
  return pair<ChrPosition, uchar *>(chrPos, dna);
}

