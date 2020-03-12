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

#include <Classifier.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include "utils.h"

using namespace std;

/****************************************
*         SingleReadClassifier          *
****************************************/

SingleReadClassifier::SingleReadClassifier(Read *r, 
                      LocateOnGenome *genome, 
                      ReadIndex *tags, 
                      Parameters *parameters):
  read(r),genome(genome),tags(tags),parameters(parameters),pair_support(NULL){
}

void SingleReadClassifier::init() {
  uint *S = tags->getSupport(read);
  this->suppo = new Support(parameters, S, read, genome, tags, pair_support);
  this->taginfo = new TagInfo(genome, suppo, read);
}

SingleReadClassifier::~SingleReadClassifier() {
  if (taginfo != NULL)
    delete taginfo;
  if (suppo != NULL)
    delete suppo;
  if (read != NULL)
    delete read;
  if (sam_lines != NULL) {
    for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
      delete (*it);
    }
    delete sam_lines;
  }
}

void SingleReadClassifier::classify() {
  char *tag = read->seq;
  int nb_tag_indel;
  int nb_genome_indel;

  if (! suppo->isContinuous() 
      && ( suppo->isSingle() || suppo->isDuplicate() || suppo->isMultiple() )
      ) {
    // all cases with at least an interruption and where the read is not classify "single" or "duplicate"
    // we have to do an other pass to find the interruption
    for (uint i= 0 ; i < suppo->getNbBreaks() ; i++){
      taginfo->setCurrentBreak(i);
      SupportBreak *suppoBreak = suppo->getBreak(i);
      nb_tag_indel = suppoBreak->getNbTagIndels();
      nb_genome_indel = suppoBreak->getNbGenomeIndels();

      // We check if there is a N character causing the break
      bool found_N = false;
      for (uint j = suppoBreak->getPositionStartBreak() + suppoBreak->getThreshold() -1; j <= suppoBreak->getPositionEndBreak(); j++) {
        if(tag[j] == 'N' || tag[j] == 'n') {
          found_N = true;
          break;
        }
      }
      if(found_N) {
        taginfo->addUndeterminedError(UndeterminedErrorInfo::N_LETTERS,"Break involves a N letter(s)");
        //taginfo->addSeqErr(FIRST_SUBSTITUTION,
        //       suppoBreak->getPositionEndBreak(),0,0,
        //       suppoBreak->getScoreComputedIntraExon(),
        //       1, 
        //       getTagAtEndBreak(suppoBreak, tag, 1), 1);
        continue;
      }

      // if (suppoBreak->isDuplicated()){
      // 	cout << "DUPLICATION" << endl;
      // }

      if (! suppoBreak->hasLongEnoughTagBreak()) {
	  taginfo->addUndeterminedError(UndeterminedErrorInfo::BREAK_TOO_SMALL,"The break in the tag is too small (break length: %d)",
					suppoBreak->getTagBreakLength());
      } else if (suppoBreak->hasNoCover()) {
	if (suppoBreak->isGenomeInsertion()
	    && parameters->max_bio_ins_del <= labs(nb_genome_indel)
	    && parameters->max_splice_length >= labs(nb_genome_indel)) {
	  // Splice with little support	
	  // Todo: change the category
	  taginfo->addSpliceNoCover(nb_genome_indel);
	} else {
	  taginfo->addUndeterminedError(UndeterminedErrorInfo::NO_COVER,"The read (or the part of read) has no cover, we cannot deduce anything! (score in break:%.2f, score outside break:%.2f)",
					suppoBreak->getScoreInsideBreak(), suppoBreak->getScoreOutsideBreak());
	}
      }
      // Biological event
      else if (suppoBreak->isBiologicalEvent()) {
        // In this case we have a combination of biological reason and error
        // But we can't really say anything on that.
        if (suppoBreak->isDeviated()
	    // a splice deviated is considered as a good candidat
	    && !suppoBreak->isGenomeInsertion()) {
          taginfo->addBioUndetermined(NULL,
                                      suppoBreak->getPositionEndBreak(),
                                      BioUndeterminedInfo::COMBINATION_ERROR_BIOLOGICAL,
                                      "Combination of error and biological reason (deviation score: %.2f)",
                                      suppoBreak->getScoreInsideAverages());
        }  else if (suppoBreak->isRepeated()
	    // a splice deviated is considered as a good candidat
		    && !suppoBreak->isGenomeInsertion()) {
	// ambiguous cases (no limit start or end)
	  taginfo->addBioUndetermined(NULL, suppoBreak->getPositionStartBreak(),
        BioUndeterminedInfo::REPETITION_AROUND_THE_BREAK,
        "Repetition around the break");
	} 
	else if (suppoBreak->hasNoStartBreak()){ 
	  taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
				      suppoBreak->getPositionEndBreak(),
              BioUndeterminedInfo::NO_START_BREAK,
              "ambiguous case : no start break");
          if (parameters->deep_snp_search && !suppoBreak->hasNoEndBreak()) {
            perform_deep_snp_search(suppoBreak, i);
          }
	}	
	else if (suppoBreak->hasNoEndBreak()){
	  taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
				      suppoBreak->getPositionStartBreak(),
              BioUndeterminedInfo::NO_END_BREAK,
              "ambiguous case : no end break");
          if (parameters->deep_snp_search && !suppoBreak->hasNoStartBreak()) {
            perform_deep_snp_search(suppoBreak, i);
          }
	}      
	// unambiguous cases
	else {
	  // Check if the candidat is ambiguous or not
	  if (parameters->no_ambiguity && suppoBreak->isDuplicated()){
      BioUndeterminedInfo::Type bioundetermined_type;
      string message = "Biological event reclassication because of #no-ambiguity option. Probably several matches of the candidat on the genome.";
      if(suppoBreak->isTagSubstitution()) {
        message += " #snp";
        bioundetermined_type = BioUndeterminedInfo::AMBIGUOUS_SNP;
      } else if(suppoBreak->isTagIndel() || suppoBreak->isGenomeDeletion() ||
                (suppoBreak->isGenomeInsertion() && parameters->max_bio_ins_del > labs(nb_genome_indel))) {
        message += " #indel";
        bioundetermined_type = BioUndeterminedInfo::AMBIGUOUS_INDEL;
      } else if(suppoBreak->isGenomeInsertion() && parameters->max_splice_length > labs(nb_genome_indel)) {
        message += " #splice";
        bioundetermined_type = BioUndeterminedInfo::AMBIGUOUS_SPLICE;
      } else if(suppoBreak->isChimera()) {
        message += " #chimera";
        bioundetermined_type = BioUndeterminedInfo::AMBIGUOUS_CHIMERA;
      }
	    taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					suppoBreak->getPositionStartBreak(),
					bioundetermined_type,
          message.c_str());
	  }else{
	    // SNP 
	    if (suppoBreak->isTagSubstitution() 
		&& suppoBreak->isBiologicalIntraEvent()){
	      //char nuc;
	      // two substituions
	      if (suppoBreak->getTagBreakLength() > suppo->getThreshold()) {
		// nuc = tag[suppoBreak->getPositionEndBreak()
		// 	  -(suppoBreak->getTagBreakLength() - suppo->getThreshold())];
		//if (nuc != 'N'){
		  taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()
				      -(suppoBreak->getTagBreakLength() - suppo->getThreshold())],
				  SECOND_SUBSTITUTION);
		//}else{
		//  uint shift = suppoBreak->getTagBreakLength() - suppo->getThreshold();
		//  taginfo->addSeqErr(SECOND_SUBSTITUTION
		//		     ,suppoBreak->getPositionEndBreak()-shift,0,0,
		//		     suppoBreak->getScoreComputedIntraExon(),
		//		     1,
		//		     getTagAtEndBreak(suppoBreak, tag, 1, -shift), 1);
		//}
	      }
	      //nuc = tag[suppoBreak->getPositionEndBreak()];
	      // one substitution
	      //if (nuc != 'N'){
		taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()],
				FIRST_SUBSTITUTION);
	      //}else{
		//taginfo->addSeqErr(FIRST_SUBSTITUTION,
		//		   suppoBreak->getPositionEndBreak(),0,0,
		//		   suppoBreak->getScoreComputedIntraExon(),
		//		   1, 
		//		   getTagAtEndBreak(suppoBreak, tag, 1), 1);
	  //    }
	    } // end if snp
	    // case biological tag ins/del
	    else if (suppoBreak->isTagIndel()
		     && suppoBreak->isBiologicalIntraEvent()){
	      // snp (insertion) at potsition j-1 into the tag (a snp is one sub, one ins or one del) 
	      if (nb_tag_indel == 1){ 
		taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()], 
				INSERTION);
	      }
	      // insertion at position j-1 into the tag 
	      else if (nb_tag_indel > 0){
		taginfo->addBioTagIndel(nb_tag_indel,0);
	      }
	      // snp (deletion) at position j-1 into the tag 
	      else if (nb_tag_indel == -1){ 
		taginfo->addSNP(0, DELETION);
	      }
	      //deletion of lg ins_del_tag at position j-1 into the tag
	      else {
		taginfo->addBioTagIndel(0,labs(nb_tag_indel));
	      } 
	    } // end if tag indel
	    // case biological deletion in the genome
	    else if (suppoBreak->isGenomeDeletion()
		     && suppoBreak->isBiologicalIntraEvent()){
	      taginfo->addBioTagIndel(labs(nb_genome_indel),0);
	    } // end if genomeDeletion
	    // case biological insertion in the genome
	    else if (suppoBreak->isGenomeInsertion()){
	      // case insertion genome < max_bio_ins_del
	      if (parameters->max_bio_ins_del > labs(nb_genome_indel)){
		taginfo->addBioTagIndel(0,labs(nb_genome_indel));
	      } 
	      // case splice
	      else if (parameters->max_splice_length >= labs(nb_genome_indel)){
		taginfo->addSplice(nb_genome_indel);
	      } else {
		taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					    suppoBreak->getPositionStartBreak(),
              BioUndeterminedInfo::UNEXPECTED_CASE,
              "This case was not expected.");
	      }
	    } // end if genomeInsertion
	    // case inter/intra transplicing (chimera)
	    else if (suppoBreak->isChimera()){
        string message = "";
	      pair<float,string> stringent_info = getChimeraScore(suppoBreak);
	      // We re-classify chimera when chim_score < DEFAULT_CHIMERA_SCORE
	      // If we are using stringent chimera option, the chim_score must be greater than parameters->min_stringent_chimera_score 
	      if(stringent_info.first < parameters->min_chimera_score
		 || (parameters->stringent_chimera && stringent_info.first < parameters->min_stringent_chimera_score)){
          message = concatStringAndFloat("Chimera reclassication because of stringent-chimera option. Score is: ",stringent_info.first);
          message += " and features are: ";
          message += stringent_info.second;
		  taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					      suppoBreak->getPositionStartBreak(),
					      BioUndeterminedInfo::WEAK_CHIMERA,
                message.c_str());
	      }else{
		taginfo->addSpliceInter(stringent_info.first,stringent_info.second);
	      }
	    } // end if chimera
	    // is continuous but no known biological reasons
	    // different cases of biological undetermined reasons
	    else{
	      if (suppoBreak->getTagBreakLength() >= suppo->getThreshold())
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
              BioUndeterminedInfo::BREAK_TOO_LARGE,
              "The break in the tag is too large (%d). There may have several biological causes",
              suppoBreak->getTagBreakLength());
	      // else if (suppoBreak->getBreakLength() > suppo->getThreshold())
	      //   taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
	      // 				  suppoBreak->getPositionEndBreak(),
	      // 				  "The break is too large (%d). There may have several biological causes",
	      // 				  suppoBreak->getBreakLength());
	      else if (suppoBreak->getTagBreakLength() < parameters->min_break_length)
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
              BioUndeterminedInfo::BREAK_TOO_SMALL,
              "The break is too short (%d).",
              suppoBreak->getTagBreakLength());
	      else if (! suppoBreak->isBiologicalIntraEvent() 
		       && suppoBreak->isBiologicalInterEvent()) {
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
              BioUndeterminedInfo::LOW_SUPPORT_INTRA_EVENT,
              "Support is ok for an inter event but not for an intra.");
	      } else {
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
              BioUndeterminedInfo::UNEXPECTED_CASE,
              "This case was not expected.");
	      }
	    }
	  } 
	}
      } // end if (isBiologicalEvent)
      // --------> start of erronneous reason
      else{
	if (suppoBreak->isRepeated()) {
	  // In this case we have generally a few positions which are not located
	  // But we can't really say anything on that.
	  taginfo->addUndeterminedError(UndeterminedErrorInfo::REPETITION_AROUND_THE_BREAK,
          "pos=%d-%d Repetitions around the break",
					suppoBreak->getPositionStartBreak(),
					suppoBreak->getPositionEndBreak());
	}
	// ambiguous cases (no limit start or end)
	else if (suppoBreak->hasNoStartBreak()){
	  taginfo->addSeqErr(UNKNOWN_START_MISSING,
	  		     suppoBreak->getPositionEndBreak(),~0,~0,
	  		     suppoBreak->getScoreComputedIntraExon());	    
	}else if(suppoBreak->hasNoEndBreak()){
	  uint pos = suppoBreak->getPositionStartBreak() + suppo->getThreshold() - 1;
	  taginfo->addSeqErr(UNKNOWN_END_MISSING,
			     pos >= read->getLength() ? read->getLength()-1 : pos,~0,~0,
			     suppoBreak->getScoreComputedIntraExon());	    
	}
	// seq error ---------------------------> (unambiguous choice)
	else if (suppoBreak->isTagSubstitution()){
	  //two errors
	  if (suppoBreak->getTagBreakLength() > suppo->getThreshold()) {
	    uint shift = suppoBreak->getTagBreakLength() - suppo->getThreshold();
	    taginfo->addSeqErr(SECOND_SUBSTITUTION
			       ,suppoBreak->getPositionEndBreak()-shift,0,0,
			       suppoBreak->getScoreComputedIntraExon(),
			       1,
			       getTagAtEndBreak(suppoBreak, tag, 1, -shift), 1);
	  }
	  // one error
	  taginfo->addSeqErr(FIRST_SUBSTITUTION,
			     suppoBreak->getPositionEndBreak(),0,0,
			     suppoBreak->getScoreComputedIntraExon(),
			     1, 
			     getTagAtEndBreak(suppoBreak, tag, 1), 1);
	  
	}
	// case ins/del error into tag
	else if (suppoBreak->isTagIndel()) {
	  // insertion of lg ins_del_tag 
	  if (nb_tag_indel >= 0){
	    taginfo->addSeqErr(INSERTION,
			       suppoBreak->getPositionEndBreak(),
			       nb_tag_indel,0,
			       suppoBreak->getScoreComputedIntraExon(),
                               0,
			       getTagAtEndBreak(suppoBreak, tag, nb_tag_indel),
			       nb_tag_indel);
	  }
	  //suppression of lg ins_del_tag 
	  else{
	    taginfo->addSeqErr(DELETION
			       ,suppoBreak->getPositionEndBreak(),
			       0,labs(nb_tag_indel),
			       suppoBreak->getScoreComputedIntraExon(),
			       labs(nb_tag_indel), NULL, ~0);
	  } 
	}	
	//insertion of lg ins_del_genome 
	else if (suppoBreak->isGenomeInsertion()) {
	  if (parameters->max_bio_ins_del < labs(nb_genome_indel)
	      && parameters->max_splice_length >= labs(nb_genome_indel)
              // The following condition is there just to select
              // splices with no cover.
              // Otherwise it may be a false positive due to duplications
              // and an error
	      && suppoBreak->getScoreOutsideBreak() <= parameters->max_support_out_no_cover) {
	    // Splice with little support	
	    taginfo->addSpliceNoCover(nb_genome_indel);
	  }else {
            ChrPosition *chrPos = suppoBreak->getLocationEndBreak();
            taginfo->addUndeterminedError(UndeterminedErrorInfo::LARGE_INDEL,
                "position %d. Too many indels to be an error (probably a splice gap=%d, location=%s|%d,%u).", 
                suppoBreak->getPositionEndBreak(), 
                nb_genome_indel,suppoBreak->getChr(END_BREAK),
                suppoBreak->getStrandLocation(END_BREAK),
                chrPos->getRelativePosition());
            //delete chrPos;
          }
	}
	//suppresion of lg labs(ins_del_genome) 
	else if (suppoBreak->isGenomeDeletion()) { 
	  taginfo->addSeqErr(INSERTION,
                             suppoBreak->getPositionEndBreak()
			     ,labs(nb_genome_indel),0,
			     suppoBreak->getScoreComputedInterExon(),
			     0,
			     getTagAtEndBreak(suppoBreak, tag, labs(nb_genome_indel)),
			     labs(nb_genome_indel));
	}
	// is not continuous but no erroneous reason
	else{
	  taginfo->addUndeterminedError(UndeterminedErrorInfo::UNEXPECTED_CASE,
        "position = %d -> No location neither biological nor erroneous cause found (ins_del_tag: %d and ins_del_genome: %d, chr_before_break: %d and  chr_after_break: %d, strand_before_break: %d and strand_after_break: %d)"
					,suppoBreak->getPositionEndBreak()
					,nb_tag_indel
					,nb_genome_indel
					,suppoBreak->getChrId(START_BREAK)
					,suppoBreak->getChrId(END_BREAK)
					,suppoBreak->getStrandLocation(START_BREAK)
					,suppoBreak->getStrandLocation(END_BREAK));
	}
      }
    }
    taginfo->setCurrentBreak(suppo->getNbBreaks());
  }

  sam_lines = taginfo->getSamLines();
}

Read *SingleReadClassifier::getRead() {
  return read;
}

TagInfo *SingleReadClassifier::getTagInfo() {
  return taginfo;
}

Support *SingleReadClassifier::getSupport() {
  return suppo;
}

void SingleReadClassifier::setPairSupport(Support *pair_support){
  this->pair_support = pair_support;
}

vector<SamLine*> *SingleReadClassifier::getSamLines() {
  return sam_lines;
}

ostream &SingleReadClassifier::samOutput(ostream &os) {
  for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
    (*it)->writeLine(os);
  }
  return os;
}

int SingleReadClassifier::samOutput(samFile *out, const bam_hdr_t *h) {
  for (vector<SamLine*>::iterator it = sam_lines->begin() ; it != sam_lines->end(); ++it) {
    (*it)->writeBamRecord(out,h);
  }
  return 0;
}

void SingleReadClassifier::writeOutputs(ostream *snp
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
      , ostream *single) {

  uint i = read->id;

  // First step : the mapping output
  if (taginfo->hasNothing()){
    if (nothing != NULL)
      *nothing << i << " " << *taginfo <<  endl;
  }else{
    if (single != NULL)
      if (taginfo->isSingle())
        *single << i << " " << *taginfo << endl;
    if (duplication != NULL)
      if (taginfo->isDuplication())
        *duplication << i << " " << *taginfo <<  endl ;
    if (multiple != NULL)
      if (taginfo->isMultiple())
        *multiple << i << " " << *taginfo << endl;
    if (none != NULL)
      if (taginfo->isNone())
        *none << i << " " << *taginfo << endl;
  }
    
  //Second step : the classification output
  // if it is an erroneous tag then we uniquely write in the output errors file
  if (taginfo->getNbSeqErr() > 0) {
    if (seqErr != NULL)
      for (uint j=0; j < taginfo->getNbSeqErr(); j++) 
        *seqErr << i << " " << *taginfo->getInfosSeqErr()[j] << " " << *taginfo <<  endl ;
  }
  //     else{
  if (normal != NULL)
    if (taginfo->isNormal())
      *normal << i << " " << *taginfo <<  endl;
  if (almostNormal != NULL)
    if (taginfo->isAlmostNormal())
      *almostNormal << i << " " << *taginfo <<  endl;
  if (snp  != NULL)
    for (uint j=0; j < taginfo->getNbSNP(); j++) 
      *snp << i << " " << *taginfo->getInfosSNP()[j] << " " << *taginfo <<  endl;
  if (bioTagIndel  != NULL)
    for (uint j=0; j < taginfo->getNbBioTagIndel(); j++) 
      *bioTagIndel << i << " " << *taginfo->getInfosBioTagIndel()[j] << " " << *taginfo <<  endl;
  if (repetition != NULL)
    for (uint j=0; j < taginfo->getNbRepetition() ; j++) 
      *repetition << i << " " << *taginfo->getInfosRepetition()[j] << " " << *taginfo <<  endl;
  if (undetermined != NULL)
    for (uint j=0; j < taginfo->getNbUndeterminedError(); j++) 
      *undetermined << i << " " << *taginfo->getInfosUndeterminedError()[j] << " " << *taginfo <<  endl;
  if (splice != NULL)
    for (uint j=0; j < taginfo->getNbSplice(); j++)
      *splice << i << " " << *taginfo->getInfosSplice()[j] << " " << *taginfo <<  endl ;
  if (spliceNoCover != NULL)
    for (uint j=0; j < taginfo->getNbSpliceNoCover(); j++)
      *spliceNoCover << i << " " << *taginfo->getInfosSpliceNoCover()[j] << " " << *taginfo <<  endl;
  if (spliceInter != NULL)
    for (uint j=0; j < taginfo->getNbSpliceInter(); j++)
      *spliceInter << i << " " << *taginfo->getInfosSpliceInter()[j] << " " << *taginfo <<  endl;
  if (bioUndermined != NULL)
    for (uint j=0; j < taginfo->getNbBioUndetermined(); j++)
      *bioUndermined << i << " " << *taginfo->getInfosBioUndetermined()[j] << " " << *taginfo << endl ;
}

void SingleReadClassifier::updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable, uint *nb_single, uint *nb_duplication, uint *nb_multiple, uint *nb_none) {
  int code = taginfo->getCode();
  for (uint j=0; j < NB_MASKS; j++) {
    if (code & (1<<j))
      (*nb_classes)[j]++;
  }
  if (taginfo->isExplainable())
    (*nb_explainable)++;
  if (code == 0)
    (*nb_classes)[NB_MASKS]++;

  // We save NH information for statistics
  uint nb_locs = getSupport()->getNbLocsMax();
    if (nb_locs == 0)
      (*nb_none)++;
    else if (nb_locs == 1)
      (*nb_single)++;
    else if (nb_locs <= parameters->max_localisation_duplication)
      (*nb_duplication)++;
    else
      (*nb_multiple)++;
}

/// PRIVATE ///
char *SingleReadClassifier::getTagAtEndBreak(SupportBreak *suppoBreak, 
                            char *tag,
				                    uint length, 
                            uint shift) {
  if (length > parameters->max_bases_retrieved) 
    return NULL;

  uint pos = suppoBreak->getPositionEndBreak();
  char *dna = new char[length+1];
  dna[length] = 0;
  strncpy(dna, &tag[pos+shift-length+1], length);
  return dna;
}

void SingleReadClassifier::perform_deep_snp_search(SupportBreak *sbreak, uint break_num) {
  
  if (! sbreak->isBiologicalEvent() 
      || sbreak->getScoreComputedIntraExon() > parameters->p_value_variation_biological)
    return;
  // For now, we currently search substitutions.

  uint loc_break;               // location of the k-mer before or after the break
  uint loc_snp;                 // location of the putative SNP
  // number of nucleotides for the comparison (the minimum is specified here, 
  // but that can be larger)
  uint nb_nuc = parameters->number_nucleotides_snp_comparison; 
  
  uchar *test_qmer;      // test q-mer from the read that (may) contain
  // the SNP (the SNP is at the end of the q-mer in the case of a start break
  // and at the start of the q-mer in the case of an end break)
  uint loc_test_qmer;           /* putative location of the test q-mer on the 
                                 * genome */
  uchar *real_content;      // real DNA from the genome, including the SNP that should correspond to test_qmer
  uint pos_snp_in_qmer;         // position of the putative SNP in the test q-mer
  uint pos_of_test_qmer;        // position of the modified q-mer in the read
  int strand = 0;

  // Variables used to update the break's informations
  uint pos_start_break;
  uint pos_end_break;
  uint loc_start_break;
  uint loc_end_break;

  //const char *read = taginfo->getTag();
  uint read_length = read->getLength();
  if (sbreak->hasNoStartBreak() 
      && sbreak->getPositionEndBreak() >= nb_nuc ) {
    
    nb_nuc = max(nb_nuc, sbreak->getPositionEndBreak());
    test_qmer = new uchar[nb_nuc+2];

    // End break's most likely location
    // SNP should be one position before (on fwd strand)
    pos_end_break = sbreak->getPositionEndBreak();
    loc_end_break = loc_break = taginfo->getLocationEndBreak(break_num);
    strand = taginfo->getChrPosEndBreak(break_num)->getStrand();
    if (strand == -1) {
      loc_end_break += taginfo->getThreshold() - nb_nuc;
    }
    pos_of_test_qmer = sbreak->getPositionEndBreak() - nb_nuc;
    pos_start_break = pos_of_test_qmer + 1;
    strncpy((char *)test_qmer,
            &read->seq[pos_of_test_qmer],
            nb_nuc+1);
    pos_snp_in_qmer = nb_nuc;
    if (strand == 1) {
      loc_test_qmer = loc_break - (1 + nb_nuc);
      loc_snp = loc_break - 1;
      loc_start_break = loc_test_qmer; // starts at the same position as the test-qmer
    } else {
      loc_test_qmer = loc_break + taginfo->getThreshold();
      loc_snp = loc_test_qmer;
      loc_start_break = loc_test_qmer + 1; // ends one position before the test q-mer

    }
  } else if (sbreak->hasNoEndBreak()
             && sbreak->getPositionStartBreak() + taginfo->getThreshold() + nb_nuc <= read_length) {
    nb_nuc = read_length - (sbreak->getPositionStartBreak() + taginfo->getThreshold());
    test_qmer = new uchar[nb_nuc+2];

    loc_break = taginfo->getLocationStartBreak(break_num);
    pos_of_test_qmer = sbreak->getPositionStartBreak() + taginfo->getThreshold() - 1;
    pos_end_break = pos_of_test_qmer;
    pos_start_break = pos_end_break - nb_nuc + 1;
    strand = taginfo->getChrPosStartBreak(break_num)->getStrand();
    loc_start_break = loc_break + strand*(pos_start_break - sbreak->getPositionStartBreak());
    if (strand == -1)
      // If we are on the reverse strand, we must correct the position since the threshold changed
      // (on the reverse strand, the position depends on the threshold)
      loc_start_break += taginfo->getThreshold() - nb_nuc; 

    strncpy((char *)test_qmer,
            &read->seq[pos_of_test_qmer],
            nb_nuc+1);
    pos_snp_in_qmer = 0;

    if (strand == 1) {
      loc_test_qmer = loc_break + taginfo->getThreshold();
      loc_snp = loc_test_qmer;
      loc_end_break = loc_test_qmer + 1; // starts one position after the snp
    } else {
      loc_test_qmer = loc_break - 1 - nb_nuc;
      loc_snp = loc_break - 1;
      loc_end_break = loc_test_qmer; // ends at the same pos as the q-mer with the SNP
    }
  } else {
    return;
  }

  test_qmer[nb_nuc + 1]=0;
  real_content = genome->getGenomeSubstring(loc_test_qmer, nb_nuc+1, strand);
  
  snp_type type = SNP_UNKNOWN;
  bool match = true;
  for (uint i=0, j=0; match && i <= nb_nuc && j <= nb_nuc; i++, j++) {
    if (real_content[i] != test_qmer[j]) {
      if (i == j && i == pos_snp_in_qmer)
        type = SNP_SUBSTITUTION;
      else
        match = false;
    }
  }

  ChrPosition *pchrPos = (genome->getChrPos(loc_snp, strand));
  ChrPosition chrPos(*pchrPos);
  delete pchrPos;
  if (type == SNP_SUBSTITUTION && match) {
    taginfo->changeGenericElement(INFO_BIOLOGICAL_UNDETERMINED, INFO_SNP,
                                  MASK_BIOLOGICAL_UNDETERMINED, MASK_SNP,
                                  break_num, 0,
                                  new SNPInfo(chrPos, pos_of_test_qmer + pos_snp_in_qmer, 
                                              sbreak->getScoreComputedIntraExon(),
                                              taginfo->isDuplicated(break_num),
                                              test_qmer[pos_snp_in_qmer], 
                                              real_content[pos_snp_in_qmer]));
    // Update break's informations
    sbreak->setPositionStartBreak(pos_start_break);
    sbreak->setPositionEndBreak(pos_end_break);
    sbreak->setThreshold(nb_nuc);
    sbreak->setLocationStartBreak(loc_start_break);
    sbreak->setLocationEndBreak(loc_end_break);
    sbreak->setStrandStartBreak(strand);
    sbreak->setStrandEndBreak(strand);
    sbreak->computeScore(pos_start_break, pos_end_break);
  }
  delete [] test_qmer;
  free(real_content);
}

pair<float,string> SingleReadClassifier::getChimeraScore(SupportBreak *sbreak) {

  float score_chimera = 1;
  std::ostringstream score_detailed;
  
  // Check the size of break with higher expectations than usual
  // Theoric size : k - 1
  float score_break_length = (float)sbreak->getOriginalTagBreakLength()/((sbreak->getThreshold()-1)*1.0);
  score_detailed << "bl=" << score_break_length;

  // Check if the chimera is ambiguous
  float score_is_duplicated = 0;
  if (!sbreak->isDuplicated()){
    score_is_duplicated = 1;
  }
  score_detailed << ",id=" << score_is_duplicated;

  //Compute the coefficient of variation of the support in the break
  float score_coefficient_variation_inside = sbreak->getCoefficientVariationInside(); 
  score_detailed << ",cb=" << score_coefficient_variation_inside;

  //Compute the coefficient of variation of the support ion the read
  float score_coefficient_variation = getSupport()->getCoefficientVariation(); 
  score_detailed << ",cv=" << score_coefficient_variation;
  
  
  //Score chimera
  // Extraction of rules from random forest model (inTRees R package)
  // (version 2.0)
  // len freq    err    
  // [1,] "3" "0.333" "0.045"
  // [2,] "4" "0.114" "0.036"
  // [3,] "2" "0.09"  "0.185"
  // [4,] "4" "0.014" "0.153"
  // [5,] "1" "0.577" "0.081"
  // condition                                                      pred   
  // [1,] "bl>0.9761905 & cb<=0.6316135 & id>0.5"                        "TRUE" 
  // [2,] "bl>0.9285715 & cb<=0.234149 & cb>0.04478385 & id>0.5"         "TRUE" 
  // [3,] "bl>0.8333335 & cb<=0.01267865"                                "TRUE" 
  // [4,] "bl>0.8333335 & bl<=0.9285715 & cb<=0.01625305 & cb>0.0119848" "TRUE" 
  // [5,] "bl<=0.9761905"                                                "FALSE"
  // impRRF              
  // [1,] "1"                 
  // [2,] "0.0705675317455272"
  // [3,] "0.0501957113882427"
  // [4,] "0.0189422300867082"
  // [5,] "0.0111771619880209"
  //
  // (version 2.1)
  // len freq    err
  // [1,] "3" "0.333" "0.045"
  // [2,] "4" "0.114" "0.036"
  // [3,] "2" "0.483" "0.124"
  // [4,] "5" "0.208" "0.152"
  // [5,] "2" "0.012" "0"
  // [6,] "2" "0.03"  "0"
  //       pred    impRRF
  // [1,] "TRUE"  "1"                 
  // [2,] "TRUE"  "0.0705675317455272"
  // [3,] "TRUE"  "1"
  // [4,] "TRUE"  "0.0330776033257221"
  // [5,] "FALSE" "0.0718883203646134"
  // [6,] "FALSE" "0.164602339595228"
  // condition
  // [1,] "bl>0.9761905 & cb<=0.6316135 & id>0.5"
  // [2,] "bl>0.9285715 & cb<=0.234149 & cb>0.04478385 & id>0.5"
  // [3,] "bl>0.9761905 & id>0.5"
  // [4,] "bl>0.9285715 & cb>0.04063885 & cb<=0.526224 & cv>0.2247855 & id>0.5"
  // [5,] "cb<=0.0999676 & cb>0.0984417"
  // [6,] "cb>0.0648697 & cb<=0.07375785"

  
  if (score_break_length>0.9761905 && score_coefficient_variation_inside<=0.6316135 && score_is_duplicated==1) {
    //Then TRUE
    score_chimera = ((((1.0-score_break_length)+score_coefficient_variation_inside)/2.0));
  }   
  else if (score_break_length>0.9285715 && score_coefficient_variation_inside<=0.234149 && score_coefficient_variation_inside>0.04478385 && score_is_duplicated==1) {
    //Then TRUE
    score_chimera = (((1.0-score_break_length)+score_coefficient_variation_inside)/2.0);
  }
  else if (score_break_length>0.9761905 && score_is_duplicated==1) {
    //Then TRUE
    score_chimera = ((1.0-4*score_break_length)/2.0);
  }   
  else if (score_break_length>0.9285715 && score_coefficient_variation_inside<=0.526224 && score_coefficient_variation_inside>0.04063885 && score_coefficient_variation>0.2247855 && score_is_duplicated==1) {
    //Then TRUE
    score_chimera = (((1.0-score_break_length)+score_coefficient_variation_inside)/2.0);
  }
  else {
    //Then FALSE
    if (score_is_duplicated == 1)
      score_chimera = (1.0-2*pow(score_break_length,2));    
  }
  
  //Reverse value for min stringent chimera
  score_chimera = 1-score_chimera;

  //Some extreme value are superior to 1
  if(score_chimera>1) {
    score_chimera=1;
  }
  
  return make_pair(score_chimera,score_detailed.str());
}

/****************************************
*       PairedEndReadClassifier         *
****************************************/

PairedEndReadClassifier::PairedEndReadClassifier(Read *r1, Read *r2,
                         LocateOnGenome *genome, 
                         ReadIndex *tags, 
                         Parameters *parameters):
  r1(r1),r2(r2),tags(tags),genome(genome),parameters(parameters) {
    this->paired_end_chimera = 0;
    classifier1 = new SingleReadClassifier(r1,genome,tags,parameters);
    classifier2 = new SingleReadClassifier(r2,genome,tags,parameters);
}

void PairedEndReadClassifier::init() {
  classifier1->init();
  classifier2->setPairSupport(classifier1->getSupport());
  classifier2->init();
}

PairedEndReadClassifier::~PairedEndReadClassifier() {
  delete classifier1;
  delete classifier2;
}

void PairedEndReadClassifier::classify() {
  classifier1->classify();
  classifier2->classify();

  TagInfo *taginfo1 = classifier1->getTagInfo();
  TagInfo *taginfo2 = classifier2->getTagInfo();

  // reclassement of the chimera(s) in the first read
  uint i = 0;
  while (i < taginfo1->getNbSpliceInter()) {
    if(!chimeraPairedEndCheck(taginfo1->getInfosSpliceInter()[i],classifier2)) {
      taginfo1->removeSpliceInter(i,BioUndeterminedInfo::DISCORDANT_PAIR,"Chimera reclassification because paired-end informations does not match the chimera.");
    } else {
      i++;
    }
  }

  // reclassement of the chimera(s) in the second read
  i = 0;
  while (i < taginfo2->getNbSpliceInter()) {
    if(!chimeraPairedEndCheck(taginfo2->getInfosSpliceInter()[i],classifier1)) {
      taginfo2->removeSpliceInter(i,BioUndeterminedInfo::DISCORDANT_PAIR,"Chimera reclassification because paired-end informations does not match the chimera.");
    } else {
      i++;
    }
  }
  
  // Looking for paired-end chimera
  if(taginfo1->isSingle() && taginfo2->isSingle() && !taginfo1->hasSpliceInterChr() && !taginfo2->hasSpliceInterChr()) {
    ChrPosition *paire1 = taginfo1->getLocation();
    ChrPosition *paire2 = taginfo2->getLocation();

    if(paire1 != NULL && paire2 != NULL) {

      if(strcmp(paire1->getChrPosition(), paire2->getChrPosition()) != 0) {
        paired_end_chimera = 1;
      } else {
        if(paire1->getStrand() != paire2->getStrand()) {
          long long int relativeDist = labs((long long int)paire1->getRelativePosition() - (long long int)paire2->getRelativePosition());
          if((paire1->getStrand() == 1) && (paire1->getRelativePosition() > paire2->getRelativePosition() && relativeDist > taginfo2->getSupportLength())) {
            paired_end_chimera = 3;
          } else if((paire2->getStrand() == 1) && (paire2->getRelativePosition() > paire1->getRelativePosition() && relativeDist > taginfo1->getSupportLength())) {
            paired_end_chimera = 3;
          } else if(relativeDist > parameters->max_splice_length) {
            paired_end_chimera = 2;
          }
        } else { // Paired-end read should be located on opposite strand according to the protocol
          paired_end_chimera = 4;
        }
      }
    }
  }
}

ostream &PairedEndReadClassifier::samOutput(ostream &os) {
  vector<SamLine*> sam_lines = getSamLines();

  for (vector<SamLine*>::iterator it = sam_lines.begin() ; it != sam_lines.end(); ++it) {
    (*it)->writeLine(os);
  }

  return os;
}

int PairedEndReadClassifier::samOutput(samFile *out, const bam_hdr_t *h) {
  vector<SamLine*> sam_lines = getSamLines();

  for (vector<SamLine*>::iterator it = sam_lines.begin() ; it != sam_lines.end(); ++it) {
    (*it)->writeBamRecord(out,h);
  }
  return 0;
}

vector<SamLine*> PairedEndReadClassifier::getSamLines() {
  vector<SamLine*> sam_lines;

  vector<SamLine*> *sam_lines1 = classifier1->getSamLines();
  vector<SamLine*> *sam_lines2 = classifier2->getSamLines();

  SamLine *primary_line1 = (*sam_lines1)[0];
  SamLine *primary_line2 = (*sam_lines2)[0];

  postProcessSamLines(*sam_lines1,*primary_line1,*primary_line2,true);
  postProcessSamLines(*sam_lines2,*primary_line2,*primary_line1,false);

  sam_lines.insert(sam_lines.end(), sam_lines1->begin(), sam_lines1->end());
  sam_lines.insert(sam_lines.end(), sam_lines2->begin(), sam_lines2->end());

  return sam_lines;
}

void PairedEndReadClassifier::setPairedEndOptionalFields(SamLine &line, const SamLine &paired_line) {
  ostringstream string_stream;

  //line.addOptionalField("R2",paired_line.getSeq());
  string_stream << paired_line.getCigar();
  line.addOptionalField("MC",string_stream.str());
  line.addOptionalField("MQ",paired_line.getMapQ());
}

void PairedEndReadClassifier::postProcessSamLines(vector<SamLine*> &sam_lines, SamLine &primary_line, SamLine &paired_primary_line, bool is_first_taginfo) {

  setPairedEndOptionalFields(primary_line, paired_primary_line);
  setPairedEndOptionalFields(paired_primary_line, primary_line);
  
  // If PE read is mapped we print paired-end additional fields
  if(!paired_primary_line.isSegmentUnmapped()) {
    ostringstream additionalInfos;
    additionalInfos << "loc:" << getSecondTagInfo()->isSingle() << ':' << getSecondTagInfo()->isDuplication() << ':' << getSecondTagInfo()->isMultiple();
    // If there is a chimera between both reads, we print that
    // piece of information
    if(this->hasPairedEndChimera()) {
      additionalInfos << ";chimera:" << *getFirstTagInfo()->getLocation() << ':' << *getSecondTagInfo()->getLocation();
    }
    primary_line.addOptionalField("XP",additionalInfos.str());
  }

  // We add informations about paired-end in the SAM lines
  // Setting flag bits for paired-end
  for (vector<SamLine*>::iterator it = sam_lines.begin() ; it != sam_lines.end(); ++it) {
    // Set FLAG
    if(is_first_taginfo) {
      (*it)->setFirstSegmentInTheTemplate();
    } else {
      (*it)->setLastSegmentInTheTemplate();
    }

    if ((is_first_taginfo &&
	 !getFirstTagInfo()->isNone() &&
	 (getFirstTagInfo()->hasSplice() || getFirstTagInfo()->hasSpliceNoCover() || getFirstTagInfo()->hasSpliceIntraChr()))
	||
	(!is_first_taginfo &&
	 !getSecondTagInfo()->isNone() &&
	 (getSecondTagInfo()->hasSplice() || getSecondTagInfo()->hasSpliceNoCover() || getSecondTagInfo()->hasSpliceIntraChr()))
	){
      // update XS field
      switch (parameters->paired_end_orientation) {
      case FORWARD_REVERSE:
	if ((!(*it)->isSeqReverseComplemented() && (*it)->isFirstSegmentInTheTemplate())
	    || ((*it)->isSeqReverseComplemented() && !(*it)->isFirstSegmentInTheTemplate()))
	  (*it)->addOptionalField("XS",'+');
	else
	  (*it)->addOptionalField("XS",'-');
	break;
      case REVERSE_FORWARD:
	if ((!(*it)->isSeqReverseComplemented() && !(*it)->isFirstSegmentInTheTemplate())
	    || ((*it)->isSeqReverseComplemented() && (*it)->isFirstSegmentInTheTemplate()))
	  (*it)->addOptionalField("XS",'+');
	else
	  (*it)->addOptionalField("XS",'-');
	break;
      case FORWARD_FORWARD:
	if ((*it)->isSeqReverseComplemented())
	  (*it)->addOptionalField("XS",'-');
	else
	  (*it)->addOptionalField("XS",'+');
	break;
      default:
	cerr << "Wrong paired_end_orientation value" << endl;
	exit(1);
      }
    }
    // Manage multiple alignments
    (*it)->setTemplateHavingMultipleSegments();
    if(primary_line.isSegmentUnmapped() || paired_primary_line.isSegmentUnmapped()) {
      (*it)->unsetEachSegmentsMapped();
    } else {
      (*it)->setEachSegmentsMapped();
    } 
    if(paired_primary_line.isSegmentUnmapped())
      (*it)->setNextSegmentUnmapped();
    if(paired_primary_line.isSeqReverseComplemented())
      (*it)->setNextSeqReverseComplemented();
    // Set Pnext (field 7)
    (*it)->setPnext(paired_primary_line.getPos());
    // Set Rnext (field 8)
    if((*it)->getRname() == paired_primary_line.getRname()) {
      (*it)->setRnext("=",paired_primary_line.getRid());
    } else {
      (*it)->setRnext(paired_primary_line.getRname(),paired_primary_line.getRid());
    }
    // Set TLEN (field 9)
    // if both segment are mapped to the same reference
    if(!primary_line.isSegmentUnmapped() && !paired_primary_line.isSegmentUnmapped()
        && primary_line.getRname() == paired_primary_line.getRname()) {
      int tlen = 0;
      if(primary_line.getPos() < paired_primary_line.getPos()) {
        tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getReferenceAlignementLength()) - primary_line.getPos();
        // Confused definition of tlen, should we compute a genomic distance
        // or an exonic distance?
        //tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getNbDeletions())
        //  - (primary_line.getPos() + primary_line.getCigar().getNbPositionsRead());
      } else if (primary_line.getPos() > paired_primary_line.getPos()) {
        tlen = paired_primary_line.getPos() - (primary_line.getPos() + primary_line.getCigar().getReferenceAlignementLength());
        //tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getNbPositionsRead())
        //  - (primary_line.getPos() + primary_line.getCigar().getNbDeletions());
      }
      (*it)->setTlen(tlen);
      //primary_line1->setTlen(primary_line1->getPos() - (primary_line2->getPos() + primary_line2->getCigar().getReferenceAlignementLength())
      //primary_line2->setTlen(primary_line2->getPos() - (primary_line1->getPos() + primary_line1->getCigar().getReferenceAlignementLength())
    }
    //(*it)->writeLine(os);
  }
}

void PairedEndReadClassifier::writeOutputs(ostream *snp
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
      , ostream *single) {

  classifier1->writeOutputs(snp,
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

  classifier2->writeOutputs(snp,
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
}

void PairedEndReadClassifier::updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable, uint *nb_single, uint *nb_duplication, uint *nb_multiple, uint *nb_none) {
  classifier1->updateStatistics(nb_classes, nb_explainable, nb_single, nb_duplication, nb_multiple, nb_none);
  classifier2->updateStatistics(nb_classes, nb_explainable, nb_single, nb_duplication, nb_multiple, nb_none);
  if(this->hasPairedEndChimera()) {
    (*nb_classes)[posBitInMask(MASK_PAIRED_END_CHIMERA)]++;
  }
}

TagInfo *PairedEndReadClassifier::getFirstTagInfo() {
  return classifier1->getTagInfo();
}

TagInfo *PairedEndReadClassifier::getSecondTagInfo() {
  return classifier2->getTagInfo();
}

bool PairedEndReadClassifier::hasPairedEndChimera() {
  return paired_end_chimera > 0;
}

void PairedEndReadClassifier::writePairedEndChimera(ostream *pairedEndChimera) {
  // If there is a chimera between the paired-end tags
  if(this->hasPairedEndChimera()) {
    *pairedEndChimera << paired_end_chimera << " " << classifier1->getRead()->id << " " << *getFirstTagInfo() << " " << classifier2->getRead()->id << " " << *getSecondTagInfo() << endl;
  }
}

/// PRIVATE ///
bool PairedEndReadClassifier::chimeraPairedEndCheck(SpliceInterInfo *chimera,
                              SingleReadClassifier *pairedEndClassifier) {
  TagInfo *pairedEndTag = pairedEndClassifier->getTagInfo();
  ChrPosition P1 = chimera->getChrPosition(); // Position of the chimera before the break
  ChrPosition P2 = chimera->getChromosomeDest(); // Position of the chimera before the break
  ChrPosition *pairePos = NULL;
  
  if (pairedEndTag->isSingle())
    pairePos = pairedEndTag->getLocation(); // Position of the paired-end read only if it is the once (otherwise it is crap!!!)

  bool isAValidChimera = false;

  // If the paired-end read is located
  if (pairePos != NULL) {
    bool chimera_class2 = true;
    if (pairedEndTag->hasSpliceInterChr()){
      SpliceInterInfo **chimeraInfos = pairedEndTag->getInfosSpliceInter();
      for (uint i=0 ; i < pairedEndTag->getNbSpliceInter() && chimera_class2 ; i++){
	chimera_class2 = (chimeraInfos[i]->getChimeraClass() == 2);
      }
    }
    // if it only has chimeras of class 2
    if (chimera_class2){
      // ------------------------------------
      // 1. We use P1 position of the chimera
      // ------------------------------------
      // a) P1(+) Case
      //                                        ^ (P1)
      //                                 -------|\\\\\|--------- strand : +
      // ======================================================= (ADN)
      // -------------------- (paired-end read)                  strand : -
      if (P1.getStrand() == 1) {
	if (P1.getStrand() != pairePos->getStrand() &&
	    strcmp(P1.getChrPosition(), pairePos->getChrPosition()) == 0 &&
	    P1.getRelativePosition() >= pairePos->getRelativePosition() &&
	    (P1.getRelativePosition() - pairePos->getRelativePosition()) <= parameters->max_splice_length) {
	  isAValidChimera = true;
	}
      } else if (!isAValidChimera) {
	// b) P1(-) Case
	//                  (paired-end read) -------------------- strand : +
	// ======================================================= (ADN)
	// -------|\\\\\|---------                                 strand : -
	//              ^ (P1)
	if (P1.getStrand() != pairePos->getStrand() &&
	    strcmp(P1.getChrPosition(), pairePos->getChrPosition()) == 0 &&
	    P1.getRelativePosition() <= pairePos->getRelativePosition() &&
	    (pairePos->getRelativePosition() - P1.getRelativePosition()) <= parameters->max_splice_length) {
	  isAValidChimera = true;
	}
      }
      // ------------------------------------
      // 2. We use P2 position of the chimera
      // ------------------------------------
      // a) P2(+) Case
      //              ^ (P2)
      // -------|\\\\\|---------                                 strand : +
      // ======================================================= (ADN)
      //                  (paired-end read) -------------------- strand : -
      if (P2.getStrand() == 1) {
	if (P2.getStrand() != pairePos->getStrand() &&
	    strcmp(P2.getChrPosition(), pairePos->getChrPosition()) == 0 &&
	    P2.getRelativePosition() <= pairePos->getRelativePosition() &&
	    (pairePos->getRelativePosition() - P2.getRelativePosition()) <= parameters->max_splice_length) {
	  isAValidChimera = true;
	}
      } else if (!isAValidChimera) {
	// b) P2(-) Case
	// -------------------- (paired-end read)                  strand : +
	// ======================================================= (ADN)
	//                                 -------|\\\\\|--------- strand : +
	//                                        ^ (P2)
	if (P2.getStrand() != pairePos->getStrand() &&
	    strcmp(P2.getChrPosition(), pairePos->getChrPosition()) == 0 &&
	    P2.getRelativePosition() >= pairePos->getRelativePosition() &&
	    (P2.getRelativePosition() - pairePos->getRelativePosition()) <= parameters->max_splice_length) {
	  isAValidChimera = true;
	}
      }
    } else {
      // no assumption can be made about the paired-end tag, the chimera is not discarded
      isAValidChimera = true;
    }
  }

  return isAValidChimera;
}
