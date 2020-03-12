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

#ifndef SUPPORT_H
#define SUPPORT_H
#include <list>

#include "libSSA/locateOnGenome.h"
#include "Parameters.h"
#include "SupportBreak.h"
#include "ReadIndex.h"
#include "Read.h"
#include "Bitset.h"

using namespace std;

// TODO BLURF THIS DO NOT LOOK GOOD AT ALL!!!
// Is there a really good reason for that?
// Mike? Nico? Eric? (Ok I'm just kidding for the last one...)
class TagInfo;

class Support {
  
  Parameters *parameters;	/* parameters to init constant values*/
  uint *support; 		/* the support by itself */
  Read *read;       /* the read we ara analyzing */
  uint length;			/* support length */
  LocateOnGenome *genome;	/* genome we are matching on */
  ReadIndex *tags;		/* index of the tags */

  Support *pair_support;

  uint nb_pos_located;		/* number of distinct factors that are located 
				   (at most getLength()) */
  
/*   uint min_occ;			/\* minimal number of occurrences for the support *\/ */
/*   uint max_occ;			/\* maximal number of occurrences for the support *\/ */
  uint nb_single;		/* number of factors single located */
  uint nb_multiple;		/* number of factors located many times (>1) */
  uint nb_duplicate;		// number of duplicates 
  bool almostNormal;          

  // The range support[start_pos_repeat ... end_pos_repeat] is such that
  // each value is >= MIN_OCC_REPETITION
  int start_pos_repeat;	/* start position of the repeat in the support */
  int end_pos_repeat;		/* end position of the repeat in the support */
  
  uint position_of_location;    /* Position in the read of the best k-mer locations */

  ChrPosition *best_location;
  
  ChrPosition **locations;     /* Locations of occurrences of the best k-mer 
				  and its strand (if any) */
  
  uint *nb_locs;                /* number of occurrences of each factor are 
				   located on the genome */

  uint nb_locs_max;             /* max value among the number of occurrences 
				   computed in nb_locs */

  uint nb_breaks;		/* number of times consecutive factors are not
				   located on the genome */
     
  SupportBreak **breaks;


  pair <uint, uint> *ranges_forward; /* Range of occurrences, in the FM-index, 
					on the forward strand. */
  pair <uint, uint> *ranges_reverse; /* Range of occurrences, in the FM-index, 
					on the reverse strand. */

  float average_support;        /* The average support of the read */

  
 public:

  Support(Parameters *p, uint *s, Read *r
	  ,LocateOnGenome *g, ReadIndex *i, Support *psupport = NULL);

  Support(const Support &);

  ~Support();

  // REQUESTS //

  /**
   * @pre i < j && i < nb_breaks && j < nb_breaks
   *      && i >= 0 && j >= 0
   * @return true when two breaks can be merged.
   * A break can be merged when we obvisouly have random
   * locations in the middle of a larger break.
   * This can be detected by checking if genome locations
   * of both breaks are in the same region.
   *
   * We also consider that we can merge two breaks
   * when one of them is on a border and (is a chimeric break
   * or does not have a long enough break).
   */
  bool canMergeBreak(uint i, uint j);

  /**
  * @return the average support in the read
  */
  float getAverage();

  /**
   * @return the variance of support in the read
  */
  float getVariance();

  /**
   * @return the standard deviation of support in the read
  */
  float getStandardDeviation();

  /**
  * @return the coefficient variation of support in the read
  */
  float getCoefficientVariation();
  
  /**
   * @param i: the SupportBreak to return
   * @param consider_strand: should we consider the strand when returning
   *                the break number i?
   * @pre 0 <= i < nb_breaks
   * @return the break number i in the read iff 
   *          - we don't consider the strand,
   *          - it is on forward strand,
   *          - it is unknown
   *         the break number getNbBreaks() - i - 1, otherwise.
   */
  SupportBreak *getBreak(uint i, bool consider_strand=false);

  /**
   * @return end_pos_repeat
   */
  int getEndPosRepeat();

  /**
   * @return tags 
   */
  ReadIndex *getIndexTags();

  /**
   * @return genome
   */
  LocateOnGenome *getGenome();
  
  /**
   * @return the support length
   */
  uint getLength();

  /**
   * @return min_occ
   * @post min_occ > 0 iff getNbSingle()+getNbMultiple() == getLength()
   */
  /* uint getMinOcc();	        */

  /**
   * @return max_occ
   * @post max_occ == 0 iff isNone()
   */
/*   uint getMaxOcc(); */

  /**
   * @return nb_breaks
   */
  uint getNbBreaks();

  /**
   * @return nb_duplicate
   * @post getNbDuplicate() <= getNbMultiple()
   */
  uint getNbDuplicate();

  /**
   * @return nb_locs
   */
  uint *getNbLocs();

  /**
   * @return nb_locs at position i
   */
  uint getNbLocs(uint i);

  /**
   * @return the max value in nb_locs
   */
  uint getNbLocsMax();

  /**
   * @return nb_multiple
   */
  uint getNbMultiple();

  /**
   * @return the number of distinct positions where a factor of length
   *         getThreshold() is located.
   */
  uint getNbPositionsLocated();

  /**
   * @return nb_single
   */
  uint getNbSingle();

  /**
   * @return parameters
   */
  Parameters *getParameters();

  /**
   * @return the 'best' location of a factor (of length getThreshold()) in the
   *         read. 
   *         This is the location in the middle of a "run" of 
   *         single-located factors. If no such location is found
   *         take a location in the middle of a run of located factors.
   *         If nothing is found return NULL
   */
  ChrPosition *getBestLocation();
  
  /**
   * getLocation() is the position of the k-mer at position
   * getPositionOfLocation() in the read.
   * Return read->getLength() if getLocation() == NULL
   */
  uint getPositionOfLocation();

  /**
   * @return position_of_location != read->getLength() 
   */ 
  bool hasPositionOfLocation();

  /**
   * @return the ChrPos location at pos i in locations at position_of_location;
   */
  ChrPosition *getLocation(uint i);
  
  /**
   * @return all the ChrPos locations at position_of_location;
   */
  ChrPosition **getLocations();

  /**
   * @post best_location = new ChrPosition(pos, chr, strand);
   */
  void setBestLocation(ChrPosition *best_location);
    
  /** 
   * Return the range of the k-mer at position i of
   * the profile on a specific strand
   *
   * @param i the number of the kmer if the support
   * @param strand the strand on wich we get the range
   * @return the range of the kmer i in the genome index
   */
  pair<uint, uint> getKmerRange(uint i, uint strand);
  
  /**
   * @return the length of the repetitioon or 0 if not
   */
  int getRepeatLength();

  /**
   * @return start_pos_repeat, -1 if no repeat exists
   */
  int getStartPosRepeat();

  /**
   * @return support
   */
  uint *getSupport();

  /**
   * @return support at position i
   */
  uint getSupport(uint i);

  /**
   * @return pair_support or NULL (in case of pair_support is not defined)
   */
  Support *getPairSupport();
  
  /**
   * @return the read object
   */
  Read *getRead();

  /**
   * @return threshold
   */
  uint getThreshold();

  /**
   * A repetition is built by a number of consecutive k-mers with 
   * a number of occurrence >= parameters->min_occ_repetition
   * @return (getRepeatLength() >= percent_min_unique_repetition*getNbPositionsLocated())
   * 
   */
  bool hasRepeat();

  /**
   * @return isAlmostNormal
   */
  bool isAlmostNormal();

  /**
   * return getLength() == getNbMultiple()+getNbSingle()
   */
  bool isContinuous();

  /**
   * @return getNbDuplicate() >= params->percent_min_duplicate*getNbPositionsLocated()
   */
  bool isDuplicate();

  /**
   * @return ! isDuplicate 
   * && (getNbSingle() < parameters->percent_max_unique_multiple*getNbPositionsLocated()
   */
  bool isMultiple();
  
  /**
   * @return (nb_single == 0 && nb_multiple == 0)
   * ||  (nb_breaks == 1 && breaks[0]->getPositionStartBreak() == 0
   *	&& breaks[0]->getPositionEndBreak() == getLength()-1)
   */
  bool isNone();
  
  /**
   * @return !isNone() && !isMultiple() && !isDuplicate()
   */ 
  bool isSingle();
 private:

  /**
   * All locations of a factor (of length getThreshold()) in
   * the read (at @parameter pos) and its number of occurences. If too much
   * locations return only getGenome()->getNbLocations() locations.
   * @post set locations
   * @post set nb_locs_max
   */
  void setLocationsAt(uint pos);

  /**
   * Check when the support drop between two consecutives value in order to 
   * define a chunck
   */
  bool isFallenSupport(uint s1, uint s2);

  /**
   * Check is the single chunck is not an artefact
   */
  bool isSingleConsistent(uint start, uint length, uint windows);

  /**
   * According to the PE protocol (fr,rf,ff) we check if the pos chosen for the 
   * mates are consistent
   */
  bool isPairedEndOrientationConsistent(ChrPosition *pos1, ChrPosition *pos2);
     
  /**
   * Check if locations oscillate inside a chunck
   * We use min, max, average, variance, standard_deviation and standard_error
   */
  bool isOscillateLoc(uint start, uint length);

  /**
   * Check if the break just before/after the chunck is long enough 
   * and greater than (k-1) by adding the chunck length
   */
  bool isGoodBreak(uint start_chunck, uint chunck_length, uint start_break, uint break_length);
  
  /**
   * @return the median pos of the chunck of representatives locations number
   */
  uint getMedianRun(uint start, uint length, bool accept_single=false);

  /**
   * Choose the position of the k-mer anchor for the mapping process
   */  
  void computeBestLocation(bool only_single, bool accept_single=false);
  
  /**
   * Will merge breaks using the following algorithm:
   * If a break is at the beginning of the read and the next one
   * ! isNiceBreak()
   *   -> merge them (if there are less than MIN_BASES_BEFORE_BREAK nt between them) 
   * If a break is at the end of the read and the previous one
   * ! isNiceBreak()
   *   -> merge them (if there are less than MIN_BASES_BEFORE_BREAK nt between them) 
   * For each short chimeric break:
   *   if isNiceMerge() is true with the next or previous break and is not too far
   *     -> merge
   *   else if the next or previous break is chimeric
   *      merge with the closest one if it is not too far
   *      (< 2 * getThreshold)
   *   else if the next or previous break is at the start or at the 
   *        end of the read
   *     -> merge
   * For each chimeric break (not yet merged):
   *   if isNiceMerge() is true with the next or previous break and is not too far
   *     or if the next break is chimeric and not too far (< 2 * getThreshold())
   *        from the current break
   *     or if the next or previous break is at the start or at the 
   *        end of the read
   *     -> merge
   * For each break:
   *    Merge with the next one iff the resulting break isVeryNiceBreak(),
   *    the current break or the next one are not isVeryNiceBreak()
   *    and if the resulting break lengthis <= getThreshold(). 
   *    We also merge in the case two consecutive breaks are overlapping
   */
  void tryToMergeBreaks();

  /**
   * Will merge locations of the representative k-mer of both paired reads
   * in order to reduce false positives k-mer located. If the intersection 
   * of locations is empty, we do nothing.
   */
  void checkPairedConnection();

  /**
   * Erase locations according to a Bitset of false_positive positions
   */
  void removeLocations(Bitset *false_positives);

  };

#endif
