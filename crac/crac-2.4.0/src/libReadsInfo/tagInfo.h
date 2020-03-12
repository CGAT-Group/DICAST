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

#ifndef TAG_INFO_H
#define TAG_INFO_H
#include <iostream>
#include <stdarg.h>
#include "SNPInfo.h"
#include "BioTagIndelInfo.h"
#include "BioUndeterminedInfo.h"
#include "DuplicationInfo.h"
#include "SpliceInfo.h"
#include "SpliceInterInfo.h"
#include "SpliceNoCoverInfo.h"
#include "SpliceIntraInfo.h"
#include "SeqErrInfo.h"
#include "repetitionInfo.h"
#include "undeterminedErrorInfo.h"
#include "genericInfo.h"
#include "../libSSA/chrPosition.h"
#include "../libSSA/locateOnGenome.h"
#include "../types.h"
#include "../Support.h"
#include "../SupportBreak.h"
#include "../Parameters.h"
#include "samLine.h"
#include "Read.h"

using namespace std;

class TagInfo {
 private:
  int code;
  LocateOnGenome *genome;
  byte *nb_each_type;
  GenericInfo ***info_each_type; /* Stores the info on breaks depending on 
                                  * the type it corresponds to */
  uint nb_causes;        /* Number of causes found in the read 
                          * INFO_DUPLICATION and INFO_REPETITION are not stored
                          * are not taken into account
                          */
  byte *types_in_order; /* Stores in which order the add* function have been 
                           called.
                           Each cell stores a type (INFO_SNP, INFO_SEQ_ERR, ...).
                           INFO_DUPLICATION and INFO_REPETITION are not stored
                           This array has nb_causes cells.
                        */
  byte *nb_causes_per_break;    /* Number of causes in each break
                                 * INFO_DUPLICATION and INFO_REPETITION are not stored
                                 * are not taken into account
                                 * This array has at least 
                                 * getNbBreaks() cells */

  Support *support;
  uint current_break;		/* The current break id processed */

  Read *read; /* The read object */

  ChrPosition **chrpos[NB_POS_BREAK];

  void addGenericElement(byte type, GenericInfo *gi);

  byte getNbGenericElement(uint type);
  GenericInfo **getGenericInfos(uint type);
  bool is_good_location;

  
 public:
  /**
   * The support is not copied but is deleted when this is deleted
   */
  TagInfo(LocateOnGenome *genome, Support *s, Read *r);

  ~TagInfo();

  float getAverageHighInside(uint);
  float getAverageLowInside(uint);
  SupportBreak *getBreak(uint, int = 1);
  /**
   * @return the type of the i-th cause 
   */
  byte getCauseType(uint i, int=1);
  ChrPosition *getChrPosEndBreak(uint);
  ChrPosition *getChrPosStartBreak(uint);
  uint getCurrentBreak();
  float getDeviationInsideBreak(uint);
  int getEndPosRepeat();
  gap_size_t getGenomeGapLength(uint);
/*   bool getIsExtended(uint); */
  uint *getLocalisations();
  ChrPosition *getLocation();
  uint getLocationEndBreak(uint);
  uint getLocationStartBreak(uint);
  uint getNbBreaks();
  int getNbCandidats(uint i);
  uint getNbCauses();
  byte getNbCausesInBreak(uint, int = 1);
  uint getNbDuplicate();
  uint getNbMultiple();
  uint getNbSingle();
  uint getPositionOfLocation();
  Read *getRead();
  int getStartPosRepeat();
  uint *getSupport();
  Support *getSupportObject();
  uint getSupportLength();
  /* uint getPositionEndBreak(uint); */
  uint getPositionEndBreak(uint, int=1);
  /* uint getPositionStartBreak(uint); */
  uint getPositionStartBreak(uint, int=1);
  float getScoreComputedIntraExon(uint);
  float getScoreComputedInterExon(uint);
  float getScoreInsideBreak(uint);
  float getScoreInsideAverages(uint);
  float getScoreOutsideBreak(uint);
  char *getTagChromosome();
  uint getTagBreakLength(uint);
  ulong getTagPosition();
  bool getTagStrand();
  uint getThreshold();

  bool isAlmostNormal();
  bool isDuplicated(uint);
  bool isExplainable();
  bool isNormal();
  bool isMultiple();
  bool isNone();
  bool isRepeated(uint);
  bool isSingle();
  bool isSupportFallingLeft(uint, int=1);
  bool isSupportFallingRight(uint, int=1);
  bool isDuplication();
  bool hasSNP();
  bool hasBioTagIndel();
  bool hasBioUndetermined();
  bool hasSeqErr();
  bool hasRepetition();
  bool hasSplice();
  bool hasSpliceIntraChr();
  bool hasSpliceInterChr();
  bool hasSpliceNoCover();
  bool hasNothing();
  bool hasUndeterminedError();
  int getCode();

  DuplicationInfo **getInfosDuplication();

  RepetitionInfo **getInfosRepetition();
  byte getNbRepetition();
  
  SNPInfo **getInfosSNP();
  byte getNbSNP();
  
  BioIndelInfo **getInfosBioTagIndel();
  byte getNbBioTagIndel();

  BioUndeterminedInfo **getInfosBioUndetermined();
  byte getNbBioUndetermined();

  SeqErrInfo **getInfosSeqErr();
  byte getNbSeqErr();

  UndeterminedErrorInfo **getInfosUndeterminedError();
  byte getNbUndeterminedError();

  SpliceInfo **getInfosSplice();
  byte getNbSplice();

  SpliceIntraInfo **getInfosSpliceIntra();
  byte getNbSpliceIntra();

  SpliceInterInfo **getInfosSpliceInter();
  byte getNbSpliceInter();

  SpliceNoCoverInfo **getInfosSpliceNoCover();
  byte getNbSpliceNoCover();

  /// COMMANDS ///

  void addDuplication();

  void addAlmostNormal();
  void addNormal();
  void addMultiple();
  void addNone();

  /**
   * position_in: position in the tag where the repetition starts (default: -1
   *              unkown)
   * position_out: position in the tag where the repetition ends.
   */
  void addRepetition(int position_in=-1, int position_out=-1);

  void addSingle();

  /**
   * tag_nuc: nucleotide found in the tag, or 0 if unknown (or deletion)
   * genome_nuc: nucleotide found in the genome, or 0 if unknown (or insertion)
   * snpID: identify which kind of SNP (depending on in which context it occurs)
   */
  void addSNP(char tag_nuc=0,
              error_context snpID=FIRST_SUBSTITUTION);

  /**
   * Another version of addSNP, less easy to use but that allows to give all
   * the necessary informations without letting the other addSNP version
   * computing them.
   * @param pos : position of the SNP in the read
   * @param base_read : nucleotide  from the read (or 0 if it doesn't exist)
   * @param base_genome: corresponding nucleotide on the genome 
   *                     (or 0 if it doesn't exist) 
   * @param chrPos : location of the SNP on the genome
   * @param score: score of the SNP.
   */
  void addSNP(uint pos, char base_read, char base_genome, ChrPosition &chrPos,
              float score);

  /**
   * nb_ins: number of insertions
   * nb_del: number of deletions
   */
  void addBioTagIndel(uint nb_ins, uint nb_del);

  /**
   * We know that there is a biological cause, but since
   * it is at the beginning or at the end of the support, it is impossible
   * to determined accurately the biological cause.
   */
  void addBioUndetermined(ChrPosition *, uint, BioUndeterminedInfo::Type, const char *, ...);
  //void addBioUndetermined(ChrPosition *, uint, );
  
  /**
   * position: position in the tag where the sequence error is
   * nb_ins: number of insertions
   * nb_del: number of deletions
   */
  void addSeqErr(error_context , uint position, uint nb_ins, uint nb_del,
		 float score,
		 uint genome_seq_length=~0,
		 char *tag_seq=NULL, uint tag_seq_length=~0);

  
  /**
   * position: position where the splice is
   * gap_length: length of the intron
   */
  void addSplice(uint gap_length);
  void addSpliceIntra(uint gap_length);
  void addSpliceNoCover(uint gap_length);

  /**
   * chimera_score: a score compute for each chimeric break 
   * according to several algorithmic filters
   */
  void addSpliceInter(float chimera_score,const string chim_score_info);

  
  void addUndeterminedError(UndeterminedErrorInfo::Type,const char *format, ...);

  /**
   * Change a GenericInfo with an other one.
   * The new GenericInfo will be inserted at the same place as the old one,
   * Indeed, it will apply to the same break.
   * This function is useful when we are reclassifying breaks afterwards.
   *
   * @param type the type of the old genericInfo to replace
   * @param newType the new type that will replace the current type
   * @param oldMask: old mask corresponding to the old cause to be replaced
   *                  (see constant whose name starts with MASK_ in const.h)
   * @param newMask: new mask corresponding to the new cause replacing the new one
   * @param i the index of the GenericInfo in info_each_type[type] array.
   * @param newGi the new GenericElement to insert in array info_each_type[newType]
   *
   * @post getNbGenericElement(type) == old getNbGenericElement(type) - 1
   *    && getNbGenericElement(newType) == old getNbGenericElement(newType) + 1
   *    && (old getNbGenericElement(type) == 1 ==> getCode() & type == 0)
   *    && getCode() & newType == newType
   */
  void changeGenericElement(byte type, byte newType, int oldMask, int newMask, 
                            uint i, GenericInfo *newGi);

  /**
   * Same function as above.
   * The difference resides in the parameter given.
   * On that version, you must provide the break ID, the cause lies in, and additionally
   * the number of a cause of that type in this break (starts at 0).
   *
   * @param type the type of the old genericInfo to replace
   * @param newType the new type that will replace the current type
   * @param oldMask: old mask corresponding to the old cause to be replaced
   *                  (see constant whose name starts with MASK_ in const.h)
   * @param newMask: new mask corresponding to the new cause replacing the new one
   * @param break_id The ID of the break where we must change something
   * @param cause_num The number of the cause (starts at 0) we want to change, having
   *                  the same type, in break break_id.
   * @param newGi the new GenericElement to insert in array info_each_type[newType]
   */
  void changeGenericElement(byte type, byte newType, int oldMask, int newMask, 
                            uint break_id, uint cause_num, GenericInfo *newGi);


  /**
   * Remove a splice Inter and add a bioUndetermined info instead.
   *
   * @param i the indice of the spliceInter if Info_each_type array
   * @param info a string containing the reason of the reclassification
   */
  void removeSpliceInter(uint i,BioUndeterminedInfo::Type, const char* info);

  /**
   * SAM output for the given read
   * XU: Is it single located?
   * XD: It it duplicated?
   * XM: Is it multiple?
   * XN: Is it normal? (1: normal, 2: almost, 0: no)
   * XL: location taken to minimise the probability of FP.
   * XP: position in the read of XL
   * XO: Location taken in Support as a single location
   * XQ: position in the read of XO
   * XC: Number of causes in the read
   * XE: Detail of a cause
   * XB: Large details for a cause
   */
  ostream &samOutput(ostream &os);

  /*
   * Return a set of sam line that corresponds to the alignement of the read.
   * Multiple line for one read can occur in two situations:
   * 1. the read has multiple locations on the genome
   * 2. the read has a chimeric alignement
   * The primary line of the read is always returned if first postion
   * of the vector. 
   */
  vector<SamLine*> *getSamLines();

  /**
   * Tells which break we are working on.
   * The given integer is the break number.
   * @post Calling one of addSNP, addGenomeIndel, ... will affect one of those causes
   * to the current break
   */
  void setCurrentBreak(uint);

  friend ostream &operator<<(ostream &os, TagInfo &);

 private:

  pair<ChrPosition, uchar *> retrieveChrPosAndNuc(error_context type,
                                                  uint nb_nuc,
                                                  int &shift);


  /**
   * This method constructs from a single SamLine object a set of SamLines
   * that corresponds to the chimericAlignement of the read if there is one.
   * This method is recursive and called by the "getSamLine()" method of this
   * very Class.
   */
  vector<SamLine*> *getChimericAlignements(SamLine *base_line, uint cigar_index, uint splice_id, uint left_softclip, uint *right_softclip);

  /**
   * Return a prototype of SAMLine that can be used to create
   * multiple alignement for one tagInfo.
   * @loc is the location of a chosen kmer in the read which represents
   * the read location
   * @pos_loc the kmer position in the read
   */
  SamLine *getSamLine(ChrPosition *loc = NULL, ulong pos_loc = 0);


  
};

ostream &operator<<(ostream &os, TagInfo &);


#endif
