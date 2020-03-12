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

#include "Parameters.h"
#include "const.h"

using namespace std;

Parameters::Parameters(uint threshold,
		       uint max_nb_located_occs,
		       uint max_splice_length,
		       uint max_bio_ins_del,
		       uint max_localisation_duplication,
		       uint min_localisation_duplication,
		       float percent_min_unique_repetition,
		       float percent_min_duplicate,
		       uint min_occ_repetition,
		       float percent_support_variation_almost_normal,
		       float p_value_variation_biological,
                       float score_intra_ambiguous,
                       float score_inter_ambiguous,
		       uint max_bases_randomly_matched,
		       uint min_bases_before_break,
		       uint max_bases_retrieved,
		       int max_extension_length,
		       uint support_score_window_length,
                       uint nb_bases_substracted_low_average,
                       uint nb_positions_check_ones,
                       uint nb_tags_info_stored,
                       uint nb_threads,
		       float min_support_no_cover,
		       float max_support_out_no_cover,
		       float max_ambiguous_average_high,
                       float min_perc_ones_inside,
                       float min_ratio_support_fall,
                       float max_nb_overlapping_breaks,
                       float percent_min_single,
                       float percent_min_multiple,
                       bool deep_snp_search,
                       uint number_nucleotides_snp_comparison,
                       bool detailed_sam,
		       uint max_verynice_merge,
		       bool stringent_chimera,
		       uint stringent_chimera_max_number_of_merges,
		       float min_stringent_chimera_score,
		       bool no_ambiguity,
		       uint treat_multiple,
		       bool show_progressbar,
		       uint reads_length,
		       uint max_mapping_comparison,
		       float percent_max_loc_variation,
		       uint max_good_paired_reads,
		       float percent_very_good_paired,
		       bool use_x_in_cigar,
		       float percent_min_chunck,
		       paired_end_orientation_t paired_end_orientation,
		       float min_chimera_score
		       ):
  threshold(threshold)
  ,max_nb_located_occs(max_nb_located_occs)
  ,max_splice_length(max_splice_length)
  ,max_bio_ins_del(max_bio_ins_del)
  ,max_localisation_duplication(max_localisation_duplication)
  ,min_localisation_duplication(min_localisation_duplication)
  ,percent_min_unique_repetition(percent_min_unique_repetition)
  ,percent_min_duplicate(percent_min_duplicate)
  ,min_occ_repetition(min_occ_repetition)
  ,percent_support_variation_almost_normal(percent_support_variation_almost_normal)
  ,p_value_variation_biological(p_value_variation_biological)
  ,score_intra_ambiguous(score_intra_ambiguous)
  ,score_inter_ambiguous(score_inter_ambiguous)
  ,max_bases_randomly_matched(max_bases_randomly_matched)
  ,min_bases_before_break(min_bases_before_break)
  ,max_bases_retrieved(max_bases_retrieved)
  ,max_extension_length(max_extension_length)
  ,support_score_window_length(support_score_window_length)
  ,nb_bases_substracted_low_average(nb_bases_substracted_low_average)
  ,nb_positions_check_ones(nb_positions_check_ones)
  ,nb_tags_info_stored(nb_tags_info_stored)
  ,nb_threads(nb_threads)
  ,min_support_no_cover(min_support_no_cover)
  ,max_support_out_no_cover(max_support_out_no_cover)
  ,max_ambiguous_average_high(max_ambiguous_average_high)
  ,min_perc_ones_inside(min_perc_ones_inside)
  ,min_ratio_support_fall(min_ratio_support_fall)
  ,max_nb_overlapping_breaks(max_nb_overlapping_breaks)
  ,percent_min_single(percent_min_single)
  ,percent_min_multiple(percent_min_multiple)
  ,deep_snp_search(deep_snp_search)
  ,number_nucleotides_snp_comparison(number_nucleotides_snp_comparison)
  ,detailed_sam(detailed_sam)
  ,max_verynice_merge(max_verynice_merge)
  ,stringent_chimera(stringent_chimera)
  ,stringent_chimera_max_number_of_merges(stringent_chimera_max_number_of_merges)
  ,min_stringent_chimera_score(min_stringent_chimera_score)
  ,no_ambiguity(no_ambiguity)
  ,treat_multiple(treat_multiple)
  ,show_progressbar(show_progressbar)
  ,reads_length(reads_length)
  ,max_mapping_comparison(max_mapping_comparison)
  ,percent_max_loc_variation(percent_max_loc_variation)
  ,max_good_paired_reads(max_good_paired_reads)
  ,percent_very_good_paired(percent_very_good_paired)
  ,use_x_in_cigar(use_x_in_cigar)
  ,percent_min_chunck(percent_min_chunck)
  ,paired_end_orientation(paired_end_orientation)
  ,min_chimera_score(min_chimera_score)
{
  if (max_nb_located_occs < min_occ_repetition){
    cerr << "max_nb_located_occs must be greater than min_loc_repetition" << endl << endl;
    exit(2);
  }
  if (reads_length > 0 && reads_length < threshold) {
    cerr << "Minimum read length (-m option) should be greater than k length (-k option)." << endl << endl;
    exit(2);
  }
}

void Parameters::init() {
  min_break_length = (uint) (MIN_BREAK_LENGTH*threshold);
}
