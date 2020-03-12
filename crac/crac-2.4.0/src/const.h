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

#ifndef CONST_H
#define CONST_H

#define MAX_NB_LOCATED_OCCS 300    /* Maximal number of occurrences we
                                    * get the actual locations */

#define MAX_MAPPING_COMPARISON 1000 /* Maximal number of comparison between two paired reads 
				      in order to choose the best k-mer location */

#define MAX_GOOD_PAIRED_READS 5 /* Number of paired reads which are located on same chr and same strand which are considered goods */

#define PERCENT_VERY_GOOD_PAIRED 0.75 /* Minimum ratio for checking very good candidats (splicing) among the good paired candidats */

#define MAX_SPLICE_LENGTH 300000 /* Maximal genomic insertion length between the locations
                                  *  before and after a break  */
#define MAX_BIO_INS_DEL 15      /* Maximal length of short biological indel */ 
#define MAX_LOCALISATION_DUPLICATION 9 /* Maximal number of occurrences for a factor
					* to be considered as duplicated 
					*/
#define MIN_LOCALISATION_DUPLICATION 2 /* Minimal number of occurrences for a factor
					* to be considered as duplicated 
					*/
#define PERCENT_MIN_UNIQUE_REPETITION 0.05 /* For classifying as repetition, minimal percentage of factors located once */


#define PERCENT_MIN_DUPLICATE 0.15  /* To be classified as duplicate, a support must 
				     * have at least this percentage of factors 
				     * located at most MAX_LOCALISATION_DUPLICATION on the genome 
				     */

#define PERCENT_MIN_SINGLE 0.15  /* To be classified as single, a support
                                  * must have at least this percent of
                                  * single-located factors 
				  */

#define PERCENT_MIN_CHUNCK 0.05  /* To be considered as a chunck, ploc
                                  * must have at least this percent of
                                  * located factors 
				  */

#define PERCENT_MIN_MULTIPLE 0.50  /* To be classified as single, a support
				    * must have at least this percent of 
				    * multiple-located factors > MAX_LOCALISATION_DUPLICATION 
				    */


#define MIN_OCC_REPETITION 20	/* Minimal number of occurrences to a located k-mer to consider a repetition */

#define PERCENT_SUPPORT_VARIATION_ALMOST_NORMAL 0.40 /* 33% */
#define PERCENT_MAX_LOC_VARIATION 0.10 /* 10% */

#define P_VALUE_VARIATION_BIOLOGICAL 0. /* P-Value which is the relation of two average
					  * to determine whether we have a sequence error 
					  * or a biological explanation. */

#define TREAT_MULTIPLE false      /* treat multiple alignments
				  (>max-duplication) rather than considering a no-alignment in the
				  SAM file */

#define NO_AMBIGUITY false      /* discard biological events (splice, snv, indel, chimera) 
				   which have several matches on the reference index */

#define SCORE_INTRA_AMBIGUOUS 0.9 /* Maximal score for intra biological causes
				     under which, the support must fall
				     to be an error. */

#define SCORE_INTER_AMBIGUOUS 0.7 /* Maximal score for inter biological causes
                                   * under which, the support must fall
                                   * to be an error. */

#define MAX_BASES_RANDOMLY_MATCHED 10 /* Maximal number of bases we may match, before the
				       * expected position (e.g. for genome ins/del,
				       * there is 0.25 probability to match a factor of
				       * the tag one position, before the one expected)*/

#define MIN_BASES_BEFORE_BREAK 3 /* Minimal number of bases we must have to consider a
				    new break */
 
#define MAX_BASES_RETRIEVED 15	/* is the number of nucleotides to display in outputfile in case of
				   insertion. */

#define MAX_EXTENSION_LENGTH 10 /* A break is extended right and left by, at most,
				   MAX_EXTENSION_LENGTH bases. */

#define SUPPORT_SCORE_WINDOW_LENGTH 16


#define MIN_SUPPORT_NO_COVER 1.3 /* Value below which the read is considered to 
				    be not supported enough for drawing conclusions.*/

#define MIN_BREAK_LENGTH 0.5 	/* a minimal break length is required for a biological cause */


#define MAX_SUPPORT_OUT_NO_COVER 3.0 /* value below which the outside score break is considered to
					be not enough for drawing conclusions. */


#define MAX_AMBIGUOUS_AVERAGE_HIGH  1.6 /* Value below which high inside average
					   of a support is considered as
					   possibly ambiguous and needs
					   extra test to ensure that there is
					   (not) an error. */

#define MAX_VERYNICE_MERGE 5000 /*  Maximal insertion length to consider a very nice 
				    merge. For example if we have two canditats: one 
				    with an insertion of 10 Kbp and an other of 40 bp,
				    so the first is nice and the second is very nice. 
				 */

#define NB_BASES_SUBSTRACTED_LOW_AVERAGE 5 /* Number of positions that are 
                                             substracted from the threshold
                                             to know how many positions we have
                                             to consider when computing the 
                                             low average inside the support */

#define NB_POSITIONS_CHECK_ONES 5 /* Designate the number of positions for which
                                    we have
                                    to count the number of ones at a border
                                    inside the support (the same number is 
                                    taken outside the support).
                                   */

#define MIN_PERCENTAGE_ONES_INSIDE 0.5 /* The minimal percentage of ones
                                         that must be in the support inside the
                                         break, near its border
                                         (0.5 is 50%)
                                        */

#define MIN_RATIO_SUPPORT_FALL 1.5 /* Minimal ratio we must have between
                                     the percentage of ones inside and
                                     outside the support.
                                     A ratio of 1.5 would mean
                                     that the percentage of ones inside must be
                                     at least 1.5 times the outside one.
                                    */

#define NB_TAGS_INFO_STORED 1000  /* Number of TagInfo stored at a given moment,
                                    per thread,
                                    for accelerating the process
                                    (for multi-threading)
                                   */

#define NB_THREADS 1 /* Number of threads by default */

#define MAX_NB_OVERLAPPING_BREAKS 2   /* Number of overlapping breaks 
                                        'allowed'. This is for merging
                                        breaks and is used to know
                                        what is the maximal distance we allow.
                                        The distance allowed will be 
                                        MAX_NB_OVERLAPPING_BREAKS * threshold 
                                       */

#define DEEP_SNP_SEARCH false /* If true, will search SNPs in the breaks that
                                are not complete at the start or at the end of
                                the read. */
#define NUMBER_NUCLEOTIDES_SNP_COMPARISON 8 /* Only used if DEEP_SNP_SEARCH is
                                              true.
                                              Corresponds to the number of extra
                                              nucleotides that are compared
                                              after (or brefore) the SNP
                                              to make sure that there is only
                                              one substitution (and that
                                              it is not a more complicated 
                                              cause) 
					     */
#define DETAILED_SAM false      /* Output detailed information for each
                                  break in the SAM output 
				 */
#define INPUT_FIFO "classify.fifo" /* Default name for the input fifo
                                     for the main server 
				    */
#define OUTPUT_FIFO "classify.out.fifo" /* Default name for the output 
                                          fifo for the main server 
					 */

#define EPSILON 0.01            /* A small value (to make the difference when
                                  comparing floats for example) 
				 */

#define DELAY_BEFORE_READ_IS_PROCESSED 1000 /* Time in microseconds to wait
                                               before a read has been
                                               processed so that we can output
                                               it */

#define NB_MASKS 19
#define MASK_NOTHING 0
#define MASK_MULTIPLE 1		/* the tag is located multiple times (and is not a DUPLICATION) */
#define MASK_CLASSIFIED 2
#define MASK_DUPLICATION 4
#define MASK_ALMOST_NORMAL 8	/* as normal but something goes down at the end, for example */
#define MASK_NORMAL 16
#define MASK_SNP 32
#define MASK_SEQ_ERR 64
#define MASK_REPETITION 128	/* a factor of the tag is repeated (but not the whole tag, apart from DUPLICATIONs) */
#define MASK_SPLICE 256
#define MASK_INTRA_TRANSPLICING 512
#define MASK_INTER_TRANSPLICING 1024
#define MASK_NONE 2048		/* not located */
#define MASK_BIOLOGICAL_UNDETERMINED 4096   /* Maybe splice or SNP
					       but we cannot conclude 
					       using the data*/
#define MASK_BIOLOGICAL_TAG_INDEL 8192 
#define MASK_UNDETERMINED_ERROR 16384  /* There is an error but we cannot
				        tell what it is really due to
				       */
#define MASK_SINGLE 32768
#define MASK_SPLICE_NO_COVER 65536 
#define MASK_PAIRED_END_CHIMERA 131072

#define INFO_NB_TYPES 11
#define INFO_SNP 0
#define INFO_SEQ_ERR 1
#define INFO_SPLICE 2
#define INFO_INTRASPLICE 3
#define INFO_INTERSPLICE 4
#define INFO_BIOLOGICAL_UNDETERMINED 5
#define INFO_REPETITION 6
#define INFO_DUPLICATION 7
#define INFO_BIOLOGICAL_TAG_INDEL 8
#define INFO_UNDETERMINED_ERROR 9
#define INFO_SPLICE_NO_COVER 10


#define NB_ELEMENTS_PER_TYPE 3  /* Internal constant for the number of cells
				   to allocate in tagInfo for info_each_type */
#define NB_INFO_STORED 6 /* Internal constant for the number of cells to allocate
			    in tagInfo for types_in_order */
#define NB_INFO_BREAKS 4        /* Internal constant for the number of cells
				   allocated in tagInfo for storing
				   the number of causes per break. */
#define NB_BREAKS_ALLOCATION_IN_SUPPORT 10
#define NB_CIGAR_ELEM 6

#define STRINGENT_CHIMERA false
#define STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES 4 /* Maximum number of merges to construct the break where the stringent chimera
                                                    has been found */

#define MIN_CHIMERA_SCORE 0.7 /* Default score to consider a chimera */
#define MIN_STRINGENT_CHIMERA_SCORE 0.9 /* Min score to consider a stringent chimera */

#define MAX_SIZE_MESSAGE_UNDETERMINED_ERROR 256

#define READS_INDEX "JELLYFISH"

#define USE_X_IN_CIGAR 0

#define PAIRED_END_ORIENTATION FORWARD_REVERSE /* 0 => fr, 1 => rf, 2 => ff */

/* Exceptions */
#define ILLEGAL_STATE_EXCEPTION 1

#endif
