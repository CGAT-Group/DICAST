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

#ifndef PARAMETERS_H
#define PARAMETERS_H 
#include <cstdlib>
#include <iostream>
#include "types.h"
#include "const.h"

class Parameters {

 public:
  uint threshold;
  uint max_nb_located_occs;
  uint max_splice_length;
  uint max_bio_ins_del;
  uint max_localisation_duplication;
  uint min_localisation_duplication;
  float percent_min_unique_repetition;
  float percent_min_duplicate;
  uint min_occ_repetition;
  float percent_support_variation_almost_normal;
  float p_value_variation_biological;
  float score_intra_ambiguous;
  float score_inter_ambiguous;
  uint max_bases_randomly_matched;
  uint min_bases_before_break;
  uint max_bases_retrieved;
  uint min_break_length;
  int max_extension_length;
  uint support_score_window_length;
  uint nb_bases_substracted_low_average;
  uint nb_positions_check_ones;
  uint nb_tags_info_stored;
  uint nb_threads;
  float min_support_no_cover;
  float max_support_out_no_cover;
  float max_ambiguous_average_high;
  float min_perc_ones_inside;
  float min_ratio_support_fall;
  float max_nb_overlapping_breaks;
  float percent_min_single;
  float percent_min_multiple;
  bool deep_snp_search;
  uint number_nucleotides_snp_comparison;
  bool detailed_sam;
  uint max_verynice_merge;
  bool stringent_chimera;
  uint stringent_chimera_max_number_of_merges;
  float min_stringent_chimera_score;
  bool no_ambiguity;
  uint treat_multiple;
  bool show_progressbar;
  uint reads_length;
  uint max_mapping_comparison; 
  float percent_max_loc_variation;
  uint max_good_paired_reads;
  float percent_very_good_paired;
  bool use_x_in_cigar;
  float percent_min_chunck;
  paired_end_orientation_t paired_end_orientation; 
  float min_chimera_score;
  
  Parameters(uint threshold=~0,
	     uint max_nb_located_occs=MAX_NB_LOCATED_OCCS,
	     uint max_splice_length=MAX_SPLICE_LENGTH,
	     uint max_bio_ins_del=MAX_BIO_INS_DEL,
	     uint max_localisation_duplication=MAX_LOCALISATION_DUPLICATION,
	     uint min_localisation_duplication=MIN_LOCALISATION_DUPLICATION,
	     float percent_min_unique_repetition=PERCENT_MIN_UNIQUE_REPETITION,
	     float percent_min_duplicate=PERCENT_MIN_DUPLICATE,
	     uint min_occ_repetition=MIN_OCC_REPETITION,
	     float percent_support_variation_almost_normal=PERCENT_SUPPORT_VARIATION_ALMOST_NORMAL,
	     float p_value_variation_biological=P_VALUE_VARIATION_BIOLOGICAL,
             float score_intra_ambiguous=SCORE_INTRA_AMBIGUOUS,
             float score_inter_ambiguous=SCORE_INTER_AMBIGUOUS,
	     uint max_bases_randomly_matched=MAX_BASES_RANDOMLY_MATCHED,
	     uint min_bases_before_break=MIN_BASES_BEFORE_BREAK,
	     uint max_bases_retrieved = MAX_BASES_RETRIEVED,
	     int max_extension_length = MAX_EXTENSION_LENGTH,
	     uint support_score_window_length = SUPPORT_SCORE_WINDOW_LENGTH,
             uint nb_bases_substracted_low_average = NB_BASES_SUBSTRACTED_LOW_AVERAGE,
             uint nb_positions_check_ones = NB_POSITIONS_CHECK_ONES,
             uint nb_tags_info_stored = NB_TAGS_INFO_STORED,
             uint nb_threads = NB_THREADS,
	     float min_support_no_cover = MIN_SUPPORT_NO_COVER,
	     float max_support_out_no_cover = MAX_SUPPORT_OUT_NO_COVER,
	     float max_ambiguous_average_high = MAX_AMBIGUOUS_AVERAGE_HIGH,
             float min_perc_ones_inside = MIN_PERCENTAGE_ONES_INSIDE,
             float min_ratio_support_fall = MIN_RATIO_SUPPORT_FALL,
             float max_nb_overlapping_breaks = MAX_NB_OVERLAPPING_BREAKS,
	     float percent_min_single = PERCENT_MIN_SINGLE,
	     float percent_min_multiple = PERCENT_MIN_MULTIPLE,
             bool deep_snp_search = DEEP_SNP_SEARCH,
             uint number_nucleotides_snp_comparison = NUMBER_NUCLEOTIDES_SNP_COMPARISON,
             bool detailed_sam = DETAILED_SAM,
	     uint max_verynice_merge = MAX_VERYNICE_MERGE,
	     bool stringent_chimera = STRINGENT_CHIMERA,
	     uint stringent_chimera_max_number_of_merges = STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES,
	     float min_stringent_chimera_score = MIN_STRINGENT_CHIMERA_SCORE,
	     bool no_ambiguity=NO_AMBIGUITY,
	     uint treat_multiple=0,
	     bool show_progressbar=true,
	     uint reads_length=0,
	     uint max_mapping_comparison=MAX_MAPPING_COMPARISON,
	     float percent_max_loc_variation=PERCENT_MAX_LOC_VARIATION,
	     uint max_good_paired_reads=MAX_GOOD_PAIRED_READS,
	     float percent_very_good_paired=PERCENT_VERY_GOOD_PAIRED,
	     bool use_x_in_cigar=USE_X_IN_CIGAR,
	     float percent_min_chunck=PERCENT_MIN_CHUNCK,
	     paired_end_orientation_t paired_end_orientation=PAIRED_END_ORIENTATION,
	     float min_chimera_score=MIN_CHIMERA_SCORE
	     );

  void init();
  
};

#endif
