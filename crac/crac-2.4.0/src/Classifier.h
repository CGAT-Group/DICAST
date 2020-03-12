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

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <iostream>
#include "ReadIndex.h"
#include "Read.h"
#include "libReadsInfo/tagInfo.h"
#include "Support.h"
#include "const.h"
#include "SupportBreak.h"
#include <config.h>

class Classifier {
protected:
  string output_sam;

public:
  virtual ~Classifier(){}

  /*
   * Classify read(s) depending of the implementation 
   * of Classifier
   */
  virtual void classify()=0;

  /*
   * Init process :
   * - creating support profil
   * - creating loc profile
   * /!\ Heavy treatment!!!
   *
   * -> Must be called before calling classify()!
   * @param pair_support the support of the mate
   */
  virtual void init() = 0;

  /**
   * Write the line(s) in SAM file
   */
  virtual ostream &samOutput(ostream &os)=0;

  /**
   * Write the line(s) in SAM file
   */
  virtual int samOutput(samFile *out, const bam_hdr_t *h)=0;


  /**
   * Print specific information about classification in
   * each file
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single)=0;

  /**
   * Update statitics about classification
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable, uint *nb_single, uint *nb_duplication, uint *nb_multiple, uint *nb_none)=0;

};

class SingleReadClassifier : public Classifier {
private:
  Read *read;
  LocateOnGenome *genome;
  ReadIndex *tags;
  Parameters *parameters;
  Support *suppo;
  Support *pair_support;
  TagInfo *taginfo;
  vector<SamLine*> *sam_lines;

public:
  /**
   * Constructor of SingleReadClassifier
   *
   * @param r the read to classify
   * @param genome the LocateOnGenome object
   * @param tags the ReadIndex object
   * @param parameters the parameters object
   */
  SingleReadClassifier(Read *r, LocateOnGenome *genome, ReadIndex *tags, Parameters *parameters);

  /*
   * Init process :
   * - creating support profil
   * - creating loc profile
   * /!\ Heavy treatment!!!
   *
   * -> Must be called before calling classify()!
   */
  virtual void init();

  /*
   * Destructor
   * - delete the read object that was passed in parameter of the constructor
   *   (/!\ no copy have been made)
   */
  ~SingleReadClassifier();

  /**
   * Classify the read contained if the classifier
   */
  virtual void classify();

  /**
   * Print the SAM line of the read in the file in arguement
   *
   * @param os the stream where the SAM line will be written
   */
  virtual ostream &samOutput(ostream &os);

  /**
   * Write the line(s) in SAM file
   */
  virtual int samOutput(samFile *out, const bam_hdr_t *h);

  /**
   * Print informations about the classification of the read
   * in the related stream
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single);

  /**
   * Update the statistics about classification with the current classified read
   *
   * @param nb_classes the array containing statistics about each type of classification
   * @param nb_explainable the number of read explainable
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable, uint *nb_single, uint *nb_duplication, uint *nb_multiple, uint *nb_none);

  /**
   * @return the read handled by the classifier
   */
  Read *getRead();

  /**
   * @return the tagInfo that has been created by the classifier.
   */
  TagInfo *getTagInfo();

  Support *getSupport();

  void setPairSupport(Support *pair_support);

  vector<SamLine*> *getSamLines();

private:


  /**
   * Return a substring of the tag of length <length> and which
   * ends at the position corresponding to the end of the break.
   * The position retrieved can be shifted using value <shift>.
   * length > params->max_bases_retrieved => result == NULL
   */
 char *getTagAtEndBreak(SupportBreak *suppoBreak, char *tag,
			uint length, uint shift=0);

  /**
   * Try to find a SNP even is some information is lacking (eg. the break
   * is not complete).
   * This is done by naively retrieving the modified nucleotide on the genome
   * and to check if we succeed to locate (a) k-mer(s) on the genome when
   * correcting the faulty nucleotide.
   * @post if a (or several) SNP(s) is/are found, the TagInfo is updated
   * accordingly.
   */
  void perform_deep_snp_search(SupportBreak *, uint break_num);

  /**
   * Apply several algorithmic filters for each chimera event.
   * This function returns of score. We consider a bioundetermined event 
   * when the score < 0.5 (ie one or more tests have failed).
   *
   * Tests are:  
   *         1 the chimera has a break_length too small
   *         2 the chimera comes from a break with too many merges
   *         3 the chimera is in a repeated zone
   *         4 the chimera has a high support variation
   *         5 the chimera is ambiguous
   *
   * @return a score between 0 and 1
   */
  pair<float,string> getChimeraScore(SupportBreak *sb);
};

class PairedEndReadClassifier : public Classifier {
private:
  Read *r1,*r2;
  ReadIndex *tags;
  LocateOnGenome *genome;
  Parameters *parameters;
  SingleReadClassifier *classifier1;
  SingleReadClassifier *classifier2;
  // 0 : nothing
  // 1 : Class 1 chimera
  // 2 : Class 2 chimera
  // 3 : Class 3 chimera
  // 4 : Class 4 chimera
  uint paired_end_chimera;

public:
  /**
   * Constructor of PairedEndReadClassifier
   *
   * @param r1 the first read of the pair to be classified
   * @param r2 the second read of the pair to be classified
   * @param genome the LocateOnGenome object
   * @param tags the ReadIndex object
   * @param parameters the parameters object
   */
  PairedEndReadClassifier(Read *r1, Read *r2, LocateOnGenome *genome, ReadIndex *tags, Parameters *parameters);

  /*
   * Init process :
   * - creating support profil
   * - creating loc profile
   * /!\ Heavy treatment!!!
   *
   * -> Must be called before calling classify()!
   */
  virtual void init();

  /*
   * Destructor
   * - Delete the two classifiers, this will also delete the two reads passed
   *   in parameters of the constructor
   */
  ~PairedEndReadClassifier();

  /**
   * Classify paired-end reads contained in the classifier
   */
  virtual void classify();

  /**
   * Print two SAM lines corresponding to the paired-end reads classified by the
   * Classifier in the output file in parameter
   *
   * @param os the stream where the SAM lines will be written
   */
  virtual ostream &samOutput(ostream &os);

  /**
   * Write the line(s) in SAM file
   */
  virtual int samOutput(samFile *out, const bam_hdr_t *h);

  /**
   * Print informations about the classification of paired-end reads
   * in the related stream
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single);
  /**
   * Update the statistics about classification with the information about paired-end reads classified
   *
   * @param nb_classes the array containing statistics about each type of classification
   * @param nb_explainable the number of read explainable
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable, uint *nb_single, uint *nb_duplication, uint *nb_multiple, uint *nb_none);

  /**
   * If the paired-end tags contains a chimera in the non sequenced region
   * between the tags informations about it in the pairedEndChimera file.
   * Otherwise nothing is done.
   */
  void writePairedEndChimera(ostream *pairedEndChimera);

  /**
   * Return true if the paired-end tags contains a chimera in the non sequenced
   * region between the tags.
   */
  bool hasPairedEndChimera();

  /**
   * @return the tagInfo that has been created by the classifier.
   */
  TagInfo *getFirstTagInfo();
  TagInfo *getSecondTagInfo();

  vector<SamLine*> getSamLines();

private:

  /**
   * Check if a chimera is valid according to the informations
   * about its paired read location.
   *
   * @param chimera the chimera to analyze
   * @param pairedEndClassifier the paired read classfier
   * @return true if chimera is valid, false if not.
   */
  bool chimeraPairedEndCheck(SpliceInterInfo *chimera, SingleReadClassifier *pairedEndClassifier);

  /*
   * Set @line optional fields deticated to store mate informations
   */
  void setPairedEndOptionalFields(SamLine &line, const SamLine &paired_line);

  /**
   * Write a vector of samlines to the output stream.  According to the primary
   * line and the paired-end read primary line the method update some flags to
   * be SAM consistent.
   */
  void postProcessSamLines(vector<SamLine*> &sam_lines, SamLine &primary_line, SamLine &paired_primary_line, bool is_first_taginfo);
};
#endif
