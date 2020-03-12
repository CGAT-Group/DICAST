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
*                            Technologie).	                                  *
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

#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>

#include "types.h"
#include "utils.h"
#include "classifyTags.h"

#include <gzstream.h>
#include <gkArrays.h>
#include "htslib/sam.h"

//#include <gkArraysTypes.h>
#ifdef HAVE_LIBGKARRAYS
#include "GkReadIndex.h"
#endif //HAVE_LIBGKARRAYS
#ifdef HAVE_LIBJELLYFISH
#include "JellyReadIndex.h"
#endif //HAVE_LIBJELLYFISH

#include "libSSA/locateOnGenome.h"
#include "libReadsInfo/samHeader.h"
#include "locateTags.h"
#include "libSSA/cracIndex.h"
#include "Parameters.h"
/* TODO 
 * for now locateServer has been disable with the abstraction of the ReadIndex
 * structure This should be resolved soon...
 */
#include "locateServer.h" 
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

using namespace std;

#define SNP 8
#define BIOTAG 9
#define SEQERR 10
#define SPLICE 11
#define SPLICE_NOCOVER 12
#define CHIMERA 13
#define UNDETERMINED 14
#define REPETITION 15
#define DUPLICATION 16
#define NOTHING 17
#define NORMAL 18
#define ALMOSTNORMAL 19
#define MULTIPLE 20
#define NONE 21
#define BIOUNDETERMINED 22
#define SINGLE 23

#define MAX_SPLICE 24
#define MAX_BIO_INDEL 25
// #define MAX_SEQERR_INDEL 26 
#define MAX_DUPLICATION 27 /* Maximal number of occurrences for a factor
					  to be considered as duplicated */
#define MIN_DUPLICATION 28 /* Minimal number of occurrences for a factor
					  to be considered as duplicated */
#define PERCENT_MIN_SINGLE_REPETITION 29 /* For classifying as repetition, minimal percentage of factors located once */
#define PERCENT_MIN_MULTIPLE_CST 30 /* To be classified as single, a support
                                  * must have at least PERCENT_MIN_SINGLE%
                                  * single-located factors 
				  */
#define PERCENT_MIN_DUPLICATION 31  /* To be classified as duplicate, a support must 
					 have at least this percentage of factors 
					 located at most MAX_LOCALISATION_DUPLICATION on the genome */
#define MIN_LOC_REPETITION 32	/* Minimal number of occurrences of a factor, for being a repetition */
#define PERCENT_VARIATION_COVERAGE 33 /* 50% */
#define P_VALUE_VARIATION_COVERAGE 34 /* P-Value which is the relation of two average
					  * to determine whether we have a sequence error 
					  * or a biological explanation. */
#define MAX_BIOLOGICAL_DEVIATION 35 /* Maximal deviation allowed for the support
					  * inside the break, to consider it as a biological event*/

#define MAX_RANDOMLY_MATCHED 36 /* Maximal number of bases we may match, before the
				       * expected position (e.g. for genome ins/del,
				       * there is 0.25 probability to match a factor of
				       * the tag one position, before the one expected)*/

#define MAX_EXTENSION 38 /* A break is extended right and left by, at most,
				 * MAX_EXTENSION_LENGTH bases. */
#define MIN_BEFORE_BREAK 39 /* Minimal number of bases we must have to consider a
				    new break */ 
#define MAX_RETRIEVED 40	/* Maximal number of bases retrieve for information
				 *  when discovering  a sequence error. */

#define SUPPORT_SCORE_WINDOW 41

#define MINIMUM_SUPPORT_NO_COVER 42

#define MINIMUM_BREAK_LENGTH 43

#define MAXIMUM_SUPPORT_OUT_NO_COVER 44

#define MAX_AMBIGUOUS_AVG_HIGH 45
#define NB_BASES_SUBSTRACTED_LOW_AVG 46
#define NB_POS_CHECK_ONES 47
#define MIN_PERC_ONES_INSIDE 48
#define MIN_RATIO_SUPPORT_FALL_CST 49
#define NB_TAGS_INFO_CST 50
#define NB_THREADS_CST 51

#define ALL_OUTPUT 52

#define MAX_NB_OVERLAPPING_BREAKS_CST 53
#define PERCENT_MIN_SINGLE_CST 54 /* To be classified as single, a support
                                  * must have at least PERCENT_MIN_SINGLE%
                                  * single-located factors 
				  */
#define SCORE_INTRA_AMBIGUOUS_CST 55
#define SCORE_INTER_AMBIGUOUS_CST 56

#define SAM 57

#define DEEP_SNP_SEARCH_CST 58
#define NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST 59
#define DETAILED_SAM_CST 60
#define SERVER_CST 61
#define INPUT_FIFO_CST 62
//63 for '?'
#define OUTPUT_FIFO_CST 64

#define MAX_VERYNICE_MERGE_CST 65
#define MAX_NB_LOCATED_OCCS_CST 66  // -n

#define READS_LENGTH 67
#define IS_STRANDED_CST 68

#define DGE 69
#define STRINGENT_CHIMERA_CST 70
#define NO_AMBIGUITY_CST 71
#define PAIRED_END_CHIMERA 72
#define GZ_OUTPUT 73
#define SUMMARY 74
#define PROGRESSBAR 75
#define TREAT_MULTIPLE_CST 76
#define STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES_CST 77
#define MIN_STRINGENT_CHIMERA_SCORE_CST 78
#define READS_INDEX_CST 79
#define MAX_MAPPING_COMPARISON_CST 80
#define USE_X_IN_CIGAR_CST 81
#define BAM 82
#define FORWARD_REVERSE_ORIENTATION_CST 83
#define REVERSE_FORWARD_ORIENTATION_CST 84
#define FORWARD_FORWARD_ORIENTATION_CST 85
#define MIN_CHIMERA_SCORE_CST 86

static const char* paired_end_orientation_options[NB_ENUM_PAIRED_END_ORIENTATION] = {"--fr","--rf","--ff"};

void print_help() {
  printf(PACKAGE_STRING);
  printf("\tCompiled on %s.\n\n",__DATE__);

  printf("   -h, --help           <none>          print this help and exit\n");
  printf("   -f, --full-help      <none>          print a complete help and exit\n");
  printf("   -v                   <none>          print version and exit\n\n");
 
  printf("Mandatory arguments\n");  
  printf("   -i                   <FILE>          set genome index file (without the extension filename)\n");
  printf("   -r                   <FILE> [FILE2]  set read file. Specify FILE2 in case of paired-end reads\n");
  printf("   -k                   <INT>           set k-mer length\n");
  printf("   -o, --sam            <FILE>          set SAM output filename or print on STDOUT with \"-o -\" argument\n\n");
  

  printf("Optional arguments\n");
  printf("  * Protocol\n");
  printf("   --fr/--rf/--ff       <none>          set the mates alignement orientation (DEFAULT %s)\n",paired_end_orientation_options[PAIRED_END_ORIENTATION]);
  printf("   --stranded           <none>          set the read mapping with for a strand specific library (DEFAULT non-strand specific)\n\n");
  printf("  * Efficiency\n");
  printf("   --nb-threads         <INT>           set the number of worker threads (DEFAULT %d)\n", NB_THREADS);
  printf("   --reads-length, -m   <INT>           set read length in case of all reads have the same length to optimize\n\
                                        CPU and memory times\n"); 
  printf("   --treat-multiple     <INT>           write alignments with multiple locations (with a fixed limit) in the SAM file rather than only one occurrence\n"); 
  printf("   --max-locs           <INT>           set the maximum number of locations on the reference index (DEFAULT %d)\n\n", MAX_NB_LOCATED_OCCS);
  printf("  * Output\n");
  printf("   --bam                <none>          sam output is encode in binary format(BAM)\n");
  printf("   --summary            <FILE>          set output summary file with some statistics about mapping and classification\n\n");
  printf("  * Accuracy\n");
  printf("   --no-ambiguity       <none>          discard biological events (splice, snv, indel, chimera) which have several matches on the reference index\n\n");

}

void print_fullhelp() {
  print_help();
  printf("\nOptional output arguments\n");
  printf("   --all                              <FILE>     set output base filename for all homemade file formats (see man page for more details)\n\
                                                 in addition to SAM output \n");
  printf("   --gz                               <none>     all output files specified after this argument are gzipped (included SAM output file)\n");
#ifdef HAVE_LIBPROGRESSBAR
  printf("   --show-progressbar                 <none>     show a progress bar for the process times on STDERR\n\n"); 
  printf("   --use-x-in-cigar                   <none>     use X cigar operator when CRAC identifies a mismatch\n\n"); 
#else
  printf("\n"); 
#endif

  // printf("  * Mapping\n");
  // printf("   --single                           <FILE>     set output single file\n"); 
  // printf("   --duplicate                        <FILE>     set output duplication file\n"); 
  // printf("   --multiple                         <FILE>     set output multiple file\n"); 
  // printf("   --none                             <FILE>     set output none file\n"); 
  // printf("   --normal                           <FILE>     set output normal file\n"); 
  // printf("   --almost-normal                    <FILE>     set output almost normal file\n\n"); 

  // printf("  * Biological causes\n");
  // printf("   --snv                              <FILE>     set output SNV file\n"); 
  // printf("   --indel                            <FILE>     set output short indel file\n");
  // printf("   --splice                           <FILE>     set output splice junction file\n"); 
  // printf("   --weak-splice                      <FILE>     set output coverless splice junction file\n"); 
  // printf("   --chimera                          <FILE>     set output chimera junction file\n"); 
  // printf("   --paired-end-chimera               <FILE>     set output for paired-end chimera file\n");
  // printf("   --biological                       <FILE>     set output bio-undetermined file\n\n");  

  // printf("  * Sequence errors\n");
  // printf("   --errors                           <FILE>     set output sequence errors file\n\n"); 

  // printf("  * Repetition\n");
  // printf("   --repeat                           <FILE>     set output repetition file\n\n"); 
   
  // printf("  * Other causes\n");
  // printf("   --undetermined                     <FILE>     set output undetermined file\n"); 
  // printf("   --nothing                          <FILE>     set output nothing file\n\n"); 

  printf("Optional process for specific research\n");
  printf("   --deep-snv                         <none>     will search hard to find SNPs\n");
  printf("   --stringent-chimera                <none>     will search chimeras with more accuracy (but less sensitivity)\n\n");

  printf("Optional process launcher (once must be selected)\n");
  printf("  * Exact matching tool\n");
  printf("   --emt                              <none>     launch CRAC-emt for exact mapping of short reads\n\n");
  printf("  * Server tool (for debugging) \n");
  printf("   --server                           <none>     launch CRAC server,the output arguments will\n\
                                                 not be taken into account\n");
  printf("   --input-name-server                <STRING>   DEFAULT %s\n", INPUT_FIFO);
  printf("   --output-name-server               <STRING>   DEFAULT %s\n\n", OUTPUT_FIFO);
  printf("Additional settings for users\n");
  printf("  * Sam output file\n");  
  printf("   --detailed-sam                     <none>     more informations are added in SAM output file\n\n"); 
 
  printf("  * Mapping\n");
  printf("   --min-percent-single-loc           <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_SINGLE);
  printf("   --min-duplication                  <INT>      DEFAULT %d\n", MIN_LOCALISATION_DUPLICATION);
  printf("   --max-duplication                  <INT>      DEFAULT %d\n", MAX_LOCALISATION_DUPLICATION);
  printf("   --min-percent-duplication-loc      <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_DUPLICATE);
  printf("   --min-percent-multiple-loc         <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_MULTIPLE);
  printf("   --min-repetition                   <INT>      DEFAULT %d\n", MIN_OCC_REPETITION);  
  printf("   --min-percent-repetition-loc       <FLOAT>    DEFAULT %.2f\n\n", PERCENT_MIN_UNIQUE_REPETITION);
  
 
  printf("  * Biological causes\n");
  printf("   --max-splice-length                <INT>      DEFAULT %d\n", MAX_SPLICE_LENGTH);
  printf("   --max-bio-indel                    <INT>      DEFAULT %d\n", MAX_BIO_INS_DEL);
  printf("   --min-chimera-score                <FLOAT>    DEFAULT %.2f\n", MIN_CHIMERA_SCORE);
  printf("   --max-bases-retrieved              <INT>      DEFAULT %d\n\n", MAX_BASES_RETRIEVED);

  printf("  * Undetermined\n");
  printf("   --min-support-no-cover             <FLOAT>    DEFAULT %.2f\n\n", MIN_SUPPORT_NO_COVER);


  printf("Additional settings for advanced users\n");
  printf("  * Break verification and fusion (merging mirage breaks)\n");
  printf("   --min-break-length                 <FLOAT>    DEFAULT %.2f\n", MIN_BREAK_LENGTH);
  printf("   --max-bases-randomly-matched       <INT>      DEFAULT %d\n", MAX_BASES_RANDOMLY_MATCHED);
  printf("   --max-extension-length             <INT>      DEFAULT %d\n\n", MAX_EXTENSION_LENGTH);

  printf("  * Threading\n");
  printf("   --nb-tags-info-stored              <INT>      DEFAULT %d\n\n", NB_TAGS_INFO_STORED);

  printf("  * Data structures\n");
  printf("   --reads-index                      <STRING>   DEFAULT %s\n\n", READS_INDEX);

  printf("  * Deep SNV search option\n");
  printf("   --nb-nucleotides-snv-comparison    <INT>      DEFAULT %d\n\n", NUMBER_NUCLEOTIDES_SNP_COMPARISON);

  printf("  * Stringent Chimera search option\n");
  printf("   --max-number-of-merges             <INT>      DEFAULT %d\n", STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES);
  printf("   --min-score-chimera-stringent      <FLOAT>    DEFAULT %.2f\n", MIN_STRINGENT_CHIMERA_SCORE);
  printf("   --max_mapping_comparison           <INT>      DEFAULT %d\n\n", MAX_MAPPING_COMPARISON);
  // printf("   --percent-variation-coverage       <FLOAT> DEFAULT %.2f\n", PERCENT_SUPPORT_VARIATION_ALMOST_NORMAL);
  // printf("   --p-value-variation-coverage       <FLOAT> DEFAULT %.2f\n", P_VALUE_VARIATION_BIOLOGICAL);
  // printf("   --score-intra-ambiguous            <FLOAT> DEFAULT %.2f\n", SCORE_INTRA_AMBIGUOUS);
  // printf("   --score-inter-ambiguous            <FLOAT> DEFAULT %.2f\n", SCORE_INTER_AMBIGUOUS);
  // printf("   --min-bases-before-break           <INT>   DEFAULT %d\n", MIN_BASES_BEFORE_BREAK);
  // printf("   --support-score-window-length      <INT>   DEFAULT %d\n", SUPPORT_SCORE_WINDOW_LENGTH);
  // printf("   --max-support-out-no-cover         <FLOAT> DEFAULT %.2f\n", MAX_SUPPORT_OUT_NO_COVER);
  // printf("   --max-ambiguous-high-average       <FLOAT> DEFAULT %.2f\n", MAX_AMBIGUOUS_AVERAGE_HIGH);
  // printf("   --max-verynice-merge               <INT>   DEFAULT %d\n", MAX_VERYNICE_MERGE);
  // printf("   --nb-bases-substracted-low         <INT>   DEFAULT %d\n", NB_BASES_SUBSTRACTED_LOW_AVERAGE);
  // printf("   --nb-positions-check-ones          <INT>   DEFAULT %d\n", NB_POSITIONS_CHECK_ONES);
  // printf("   --min-perc-ones-inside             <FLOAT> DEFAULT %.2f\n", MIN_PERCENTAGE_ONES_INSIDE);
  // printf("   --min-ratio-support-fall           <FLOAT> DEFAULT %.2f\n", MIN_RATIO_SUPPORT_FALL);
  // printf("   --max-nb-overlapping-breaks        <INT>   DEFAULT %d\n\n", MAX_NB_OVERLAPPING_BREAKS);

}

void print_version() {
  printf("CRAC version " PACKAGE_VERSION "\n\n");
  printf("/******************************************************************************\n");
  printf("*  Copyright © 2009-2016 -- LIRMM/CNRS                                        *\n");
  printf("*                           (Laboratoire d'Informatique, de Robotique et de   *\n");
  printf("*                           Microélectronique de Montpellier /                *\n");
  printf("*                           Centre National de la Recherche Scientifique)     *\n");
  printf("*                           LIFL/INRIA                                        *\n");
  printf("*                           (Laboratoire d'Informatique Fondamentale de       *\n");
  printf("*                           Lille / Institut National de Recherche en         *\n");
  printf("*                           Informatique et Automatique)                      *\n");
  printf("*                           LITIS                                             *\n");
  printf("*                           (Laboratoire d'Informatique, du Traitement de     *\n");
  printf("*                           l'Information et des Systèmes).                   *\n");
  printf("*                                                                             *\n");
  printf("*  Copyright © 2011-2016 -- IRB/INSERM                                        *\n");
  printf("*                           (Institut de Recherches en Biothérapie /          *\n");
  printf("*                           Institut National de la Santé et de la Recherche  *\n");
  printf("*                           Médicale).                                        *\n");
  printf("*                                                                             *\n");
  printf("*  Copyright © 2015-2016 -- AxLR/SATT                                         *\n");
  printf("*                           (Lanquedoc Roussilon /                            *\n");
  printf("*                            Societe d'Acceleration de Transfert de           *\n");
  printf("*                            Technologie).	                                *\n");
  printf("*                                                                             *\n");
  printf("*  Programmeurs/Progammers:                                                   *\n");
  printf("*                    Nicolas PHILIPPE <nphilippe.resear@gmail.com>            *\n");
  printf("*                    Mikaël SALSON    <mikael.salson@lifl.fr>                 *\n");
  printf("*                    Jérôme Audoux    <jerome.audoux@gmail.com>               *\n");
  printf("*   with additional contribution for the packaging of:	                *\n");
  printf("*                    Alban MANCHERON  <alban.mancheron@lirmm.fr>              *\n");
  printf("*                                                                             *\n");
  printf("*   Contact:         CRAC list   <crac-bugs@lists.gforge.inria.fr>            *\n");
  printf("*   Paper:           CRAC: An integrated RNA-Seq read analysis                *\n");
  printf("*                    Philippe N., Salson M., Commes T., Rivals E.             *\n");
  printf("*                    Genome Biology 2013; 14:R30.                             *\n");
  printf("*                                                                             *\n");
  printf("*  -------------------------------------------------------------------------  *\n");
  printf("*                                                                             *\n");
  printf("*   This File is part of the CRAC program.                                    *\n");
  printf("*                                                                             *\n");
  printf("*   This program is free software: you can redistribute it and/or modify      *\n");
  printf("*   it under the terms of the GNU General Public License as published by      *\n");
  printf("*   the Free Software Foundation, either version 3 of the License, or (at     *\n");
  printf("*   your option) any later version.  This program is distributed in the       *\n");
  printf("*   hope that it will be useful, but WITHOUT ANY WARRANTY; without even       *\n");
  printf("*   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR       *\n");
  printf("*   PURPOSE.  See the GNU General Public License for more details.  You       *\n");
  printf("*   should have received a copy of the GNU General Public License along       *\n");
  printf("*   with this program.  If not, see <http://www.gnu.org/licenses/>.           *\n");
  printf("*                                                                             *\n");
  printf("******************************************************************************/\n");
}

// Write SamHeader in the specified samFile and return a
// baù_hdr_t structure
bam_hdr_t* writeSamHeader(LocateOnGenome *genome
          , int argc
          , char **argv
          , samFile* sam_file) {

  SamHeader sam_header;
  // Display information on sequence names and lengths
  for (uint i=0; i < genome->getGenomeInfo()->getNbChr(); i++) {
    sam_header.addReference((char *)genome->getGenomeInfo()->getChrName(i),
        genome->getGenomeInfo()->getChrLength(i));
  }
  sam_header.setVersionNumber(PACKAGE_VERSION);
  for (int i=0; i < argc; i++)
    sam_header.addArgument(argv[i]);

  return sam_header.writeBamHeader(sam_file);
}

int main(int argc, char **argv) {
  fstream file;
  int nb_streams;
  ostream *streams[18] = {NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL, NULL
                          ,NULL,NULL,NULL};
  const char *extensions[18] = {
    "snp", "indel", "error", "splice", 
    "coverlessSplice", "chimera", "undetermined", "repetition", 
    "duplication", "nothing", "normal", "almostNormal", "multiple",
    "none", "bioUndetermined", "single", "pairedEndChimera", 
    "summary"
  };
  nb_streams = 18;

  bool gzOutput = false;
  readsReader *reads_reader = NULL;
  gkArrays *gk = NULL;
  ReadIndex *indexTags = NULL;
  bool is_stranded = false;
  ClassifyTags *ct;
  LocateOnGenome *log = NULL;
  Parameters *params = new Parameters();
  params->show_progressbar = false;
  uint threshold = ~0;
  char *reads = 0;
  char *reads2 = 0;
  char *indexGenome = 0;
  char *fileBasename = 0;
  char *filename = 0;
  uint max_length = 0;
  int nb_streams_to_init = 0;
  bool use_server = false;
  bool use_emt = false;
  bool use_paired_end_chimera = false;
  bool use_all_outputs = false;
  bool find_output_file = false;
  char *fifo_input_name = NULL;
  char *fifo_output_name = NULL;
  char *reads_index = NULL;
  string sam_filename;
  ostream *summary = NULL;
  samFile *sam_f = NULL;
  bool binary_sam = false;

  int option_index = 0;
  int opt;
  
  static const char *short_options =  "hi:r:o:k:m:fv";
  
  static struct option long_options[] = {
    { "help",    no_argument,       NULL, 'h' },
    { "full-help",    no_argument,       NULL, 'f' },
    { "version",    no_argument,       NULL, 'v' },
    { "all", required_argument, NULL, ALL_OUTPUT},
    { "reads-length", required_argument, NULL, READS_LENGTH},
    { "stranded", no_argument, NULL, IS_STRANDED_CST},
    { "sam",    required_argument, NULL, SAM },
    { "snv",    required_argument, NULL, SNP },
    { "indel",    required_argument, NULL, BIOTAG },
    { "errors",    required_argument, NULL, SEQERR },
    { "splice",    required_argument, NULL, SPLICE },
    { "weak-splice",    required_argument, NULL, SPLICE_NOCOVER },
    { "chimera",    required_argument, NULL, CHIMERA },
    { "undetermined",    required_argument, NULL, UNDETERMINED },
    { "repeat",    required_argument, NULL, REPETITION },
    { "duplicate",    required_argument, NULL, DUPLICATION },
    { "nothing",    required_argument, NULL, NOTHING },
    { "normal",    required_argument, NULL, NORMAL },
    { "almost-normal",    required_argument, NULL, ALMOSTNORMAL },
    { "multiple",    required_argument, NULL, MULTIPLE},
    { "none",    required_argument, NULL, NONE },
    { "biological",    required_argument, NULL, BIOUNDETERMINED },
    { "single",    required_argument, NULL, SINGLE },
    { "max-locs", required_argument, NULL, MAX_NB_LOCATED_OCCS_CST},
    { "max-splice-length",    required_argument, NULL, MAX_SPLICE },
    { "max-bio-indel",    required_argument, NULL, MAX_BIO_INDEL },
    { "min-duplication",    required_argument, NULL, MIN_DUPLICATION },
    { "max-duplication",    required_argument, NULL, MAX_DUPLICATION },
    { "min-percent-multiple-loc",    required_argument, NULL, PERCENT_MIN_MULTIPLE_CST },
    { "min-percent-duplication-loc",    required_argument, NULL, PERCENT_MIN_DUPLICATION },
    { "min-percent-single-loc", required_argument, NULL, PERCENT_MIN_SINGLE_CST},
    { "min-percent-repetition-loc",    required_argument, NULL, PERCENT_MIN_SINGLE_REPETITION },
    { "max-verynice-merge",    required_argument, NULL, MAX_VERYNICE_MERGE_CST },
    { "min-repetition",    required_argument, NULL, MIN_LOC_REPETITION },
    { "percent-variation-coverage",    required_argument, NULL, PERCENT_VARIATION_COVERAGE },
    { "p-value-variation-coverage",    required_argument, NULL, P_VALUE_VARIATION_COVERAGE },
    { "score-intra-ambiguous", required_argument, NULL, SCORE_INTRA_AMBIGUOUS_CST},
    { "score-inter-ambiguous", required_argument, NULL, SCORE_INTER_AMBIGUOUS_CST},
    { "max-bio-deviation",    required_argument, NULL, MAX_BIOLOGICAL_DEVIATION },
    { "max-bases-randomly-matched",    required_argument, NULL, MAX_RANDOMLY_MATCHED },
    { "min-bases-before-break",    required_argument, NULL, MIN_BEFORE_BREAK },
    { "max-bases-retrieved",    required_argument, NULL, MAX_RETRIEVED },
    { "max-extension-length",    required_argument, NULL, MAX_EXTENSION },
    { "support-score-window-length",    required_argument, NULL, SUPPORT_SCORE_WINDOW },
    { "min-support-no-cover",    required_argument, NULL, MINIMUM_SUPPORT_NO_COVER },
    { "min-break-length",    required_argument, NULL, MINIMUM_BREAK_LENGTH },
    { "max-support-out-no-cover",    required_argument, NULL, MAXIMUM_SUPPORT_OUT_NO_COVER },
    { "max-ambiguous-high-average",    required_argument, NULL, MAX_AMBIGUOUS_AVG_HIGH },
    { "nb-bases-substracted-low",    required_argument, NULL, NB_BASES_SUBSTRACTED_LOW_AVG },    
    { "nb-positions-check-ones",    required_argument, NULL, NB_POS_CHECK_ONES },
    { "min-perc-ones-inside",    required_argument, NULL, MIN_PERC_ONES_INSIDE },    
    { "min-ratio-support-fall",    required_argument, NULL, MIN_RATIO_SUPPORT_FALL_CST },
    { "nb-tags-info-stored", required_argument, NULL, NB_TAGS_INFO_CST},
    { "nb-threads", required_argument, NULL, NB_THREADS_CST},
    { "max-nb-overlapping-breaks", required_argument, NULL, 
      MAX_NB_OVERLAPPING_BREAKS_CST},
    { "deep-snv", no_argument, NULL, DEEP_SNP_SEARCH_CST},
    { "nb-nucleotides-snv-comparison", required_argument, NULL, 
      NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST},
    { "detailed-sam", no_argument, NULL, DETAILED_SAM_CST},
    { "server", no_argument, NULL, SERVER_CST},
    { "input-name-server", required_argument, NULL, INPUT_FIFO_CST},
    { "output-name-server", required_argument, NULL, OUTPUT_FIFO_CST},
    { "emt", no_argument, NULL, DGE },
    { "stringent-chimera", no_argument, NULL, STRINGENT_CHIMERA_CST },
    { "min-chimera-score", required_argument, NULL, MIN_CHIMERA_SCORE_CST },
    { "no-ambiguity", no_argument, NULL, NO_AMBIGUITY_CST },
    { "paired-end-chimera", required_argument, NULL, PAIRED_END_CHIMERA },
    { "gz", no_argument, NULL, GZ_OUTPUT},
    { "summary",    required_argument, NULL, SUMMARY },
#ifdef HAVE_LIBPROGRESSBAR
    { "show-progressbar", no_argument, NULL, PROGRESSBAR },
#endif
    { "treat-multiple", required_argument, NULL, TREAT_MULTIPLE_CST },
    { "max-number-of-merges", required_argument, NULL, STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES_CST },
    { "min-score-chimera-stringent", required_argument, NULL, MIN_STRINGENT_CHIMERA_SCORE_CST },
    { "reads-index", required_argument, NULL, READS_INDEX_CST},
    { "max-mapping-comparison", required_argument, NULL, MAX_MAPPING_COMPARISON_CST},
    { "use-x-in-cigar", no_argument, NULL, USE_X_IN_CIGAR_CST},
    { "fr", no_argument, NULL, FORWARD_REVERSE_ORIENTATION_CST},
    { "rf", no_argument, NULL, REVERSE_FORWARD_ORIENTATION_CST},
    { "ff", no_argument, NULL, FORWARD_FORWARD_ORIENTATION_CST},
    { "bam", no_argument, NULL, BAM},
    { 0, 0, 0, 0 } //terminator
  };
  
  while ((opt = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1){
    switch(opt) {
      
    case 'h':
      print_help();
      exit(0);
    case 'f' : 
      print_fullhelp();
      exit(0);
    case 'v' : 
      print_version();
      exit(0);
    case 'i':
      indexGenome = optarg;
      break;
    case 'r':
      reads = optarg;   
      break;
    case 'k':
      threshold = atoi(optarg); 
      params->threshold = threshold;
      break;
    case 'o':
      sam_filename = (const char *)optarg;
      //if (strcmp(sam_filename.c_str(),"-") != 0){ 
	//streams[16] = create_stream(gzOutput, optarg); 
  //    }else{
	//streams[16] = &cout;
  //    }
      break;
    case 'm':
      params->reads_length = atol(optarg);
      break;
    case ALL_OUTPUT:
      fileBasename = optarg;
      use_all_outputs = true;
      break;
    case READS_LENGTH:
      params->reads_length = atol(optarg);
      break;
    case IS_STRANDED_CST:
      is_stranded = true;
      break;
    case MAX_NB_LOCATED_OCCS_CST:
      params->max_nb_located_occs = atol(optarg);
      break;
    case SNP:
      streams[0] = create_stream(gzOutput, optarg);
      break;
    case BIOTAG:
      streams[1] = create_stream(gzOutput, optarg);
      break;
    case SEQERR:
      streams[2] = create_stream(gzOutput, optarg);
      break;
    case SPLICE:
      streams[3] = create_stream(gzOutput, optarg);
      break;
    case SPLICE_NOCOVER:
      streams[4] = create_stream(gzOutput, optarg);
      break;
    case CHIMERA:
      streams[5] = create_stream(gzOutput, optarg);
      break;
    case UNDETERMINED:
      streams[6] = create_stream(gzOutput, optarg);
      break;
    case REPETITION:
      streams[7] = create_stream(gzOutput, optarg);
      break;
    case DUPLICATION:
      streams[8] = create_stream(gzOutput, optarg);
      break;
    case NOTHING:
      streams[9] = create_stream(gzOutput, optarg);
      break;
    case NORMAL:
      streams[10] = create_stream(gzOutput, optarg);
      break;
    case ALMOSTNORMAL:
      streams[11] = create_stream(gzOutput, optarg);
      break;
    case MULTIPLE:
      streams[12] = create_stream(gzOutput, optarg);
      break;
    case NONE:
      streams[13] = create_stream(gzOutput, optarg);
      break;
    case BIOUNDETERMINED:
      streams[14] = create_stream(gzOutput, optarg);
      break;
    case SINGLE:
      streams[15] = create_stream(gzOutput, optarg);
      break;
    case SAM:
      //streams[16] = create_stream(gzOutput, optarg);
      sam_filename = (const char *)optarg;
      break;
    case PAIRED_END_CHIMERA:
      streams[16] = create_stream(gzOutput, optarg);
      use_paired_end_chimera = true;
      break;
    case SUMMARY:
      streams[17] = create_stream(gzOutput, optarg);
      summary = streams[17]; 
      break;
    case MAX_SPLICE:
      params->max_splice_length = atoi(optarg);
      break;
    case MAX_BIO_INDEL:
      params->max_bio_ins_del = atoi(optarg);
      break;
    case MAX_DUPLICATION:
      params->max_localisation_duplication = atoi(optarg);
      break;
    case MIN_DUPLICATION:
      params->min_localisation_duplication = atoi(optarg);
      break;
    case PERCENT_MIN_SINGLE_REPETITION:
      params->percent_min_unique_repetition = atof(optarg);
      break;
    case PERCENT_MIN_SINGLE_CST: 
      params->percent_min_single = atof(optarg);
      break;
    case PERCENT_MIN_MULTIPLE_CST:
      params->percent_min_multiple = atof(optarg);
      break;
    case PERCENT_MIN_DUPLICATION:
      params->percent_min_duplicate = atof(optarg);
      break;
    case MIN_LOC_REPETITION:
      params->min_occ_repetition = atoi(optarg);
      break;
    case PERCENT_VARIATION_COVERAGE:
      params->percent_support_variation_almost_normal = atof(optarg);
      break;
    case P_VALUE_VARIATION_COVERAGE:
      params->p_value_variation_biological = atof(optarg);
      break;
    case SCORE_INTRA_AMBIGUOUS_CST:
      params->score_intra_ambiguous = atof(optarg);
      break;
    case SCORE_INTER_AMBIGUOUS_CST:
      params->score_inter_ambiguous = atof(optarg);
      break;
    case MAX_RANDOMLY_MATCHED:
      params->max_bases_randomly_matched = atoi(optarg);
      break;
    case MAX_VERYNICE_MERGE_CST:
      params->max_verynice_merge = atoi(optarg);
      break;
    case MIN_BEFORE_BREAK:
      params->min_bases_before_break = atoi(optarg);
      break;
    case MAX_RETRIEVED:
      params->max_bases_retrieved = atoi(optarg);
      break;
    case MAX_EXTENSION:
      params->max_extension_length = atoi(optarg);
      break;
    case SUPPORT_SCORE_WINDOW:
      params->support_score_window_length = atoi(optarg);
      break;
    case MINIMUM_SUPPORT_NO_COVER:
      params->min_support_no_cover = atof(optarg);
      break;
    case MINIMUM_BREAK_LENGTH:
      params->min_break_length = atoi(optarg);
      break;
    case MAXIMUM_SUPPORT_OUT_NO_COVER:
      params->max_support_out_no_cover = atof(optarg);
      break;
    case MAX_AMBIGUOUS_AVG_HIGH:
      params->max_ambiguous_average_high = atof(optarg);
      break;
    case NB_BASES_SUBSTRACTED_LOW_AVG:
      params->nb_bases_substracted_low_average = atoi(optarg);
      break;
    case NB_POS_CHECK_ONES:
      params->nb_positions_check_ones = atoi(optarg);
      break;
    case MIN_PERC_ONES_INSIDE:
      params->min_perc_ones_inside = atof(optarg);
      break;
    case MIN_RATIO_SUPPORT_FALL_CST:
      params->min_ratio_support_fall = atof(optarg);
      break;
    case NB_TAGS_INFO_CST:
      params->nb_tags_info_stored = atoi(optarg);
      break;
    case NB_THREADS_CST:
      params->nb_threads = atoi(optarg);
      break;
    case MAX_NB_OVERLAPPING_BREAKS_CST:
      params->max_nb_overlapping_breaks = atof(optarg);
      break;
    case DEEP_SNP_SEARCH_CST:
      params->deep_snp_search = true;
      break;
    case NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST:
      params->number_nucleotides_snp_comparison = atoi(optarg);
      break;
    case DETAILED_SAM_CST:
      params->detailed_sam = true;
      break;
    case SERVER_CST:
      use_server=true;
      break;
    case STRINGENT_CHIMERA_CST:
      params->stringent_chimera = true;
      break;
    case INPUT_FIFO_CST:
      fifo_input_name = new char[strlen(optarg)+1];
      strcpy(fifo_input_name, optarg);
      break;
    case OUTPUT_FIFO_CST:
      fifo_output_name = new char[strlen(optarg)+1];
      strcpy(fifo_output_name, optarg);
      break;
    case DGE:
      use_emt=true;
      break;
    case GZ_OUTPUT:
      gzOutput = true;
      break;
    case NO_AMBIGUITY_CST:
      params->no_ambiguity=true;
      break;
    case PROGRESSBAR:
      params->show_progressbar = true;
      break;
    case TREAT_MULTIPLE_CST:
      params->treat_multiple= atoi(optarg);
      break;
    case STRINGENT_CHIMERA_MAX_NUMBER_OF_MERGES_CST:
      params->stringent_chimera_max_number_of_merges = atoi(optarg);
      break;
    case MIN_STRINGENT_CHIMERA_SCORE_CST:
      params->min_stringent_chimera_score = atof(optarg);
      break;
    case READS_INDEX_CST:
      reads_index = new char[strlen(optarg)+1];
      strcpy(reads_index, optarg);
      break;
    case MAX_MAPPING_COMPARISON_CST:
      params->max_mapping_comparison = atoi(optarg);
      break;
    case USE_X_IN_CIGAR_CST:
      params->use_x_in_cigar = true;
      break;
    case BAM:
      binary_sam = true;
      break;
    case FORWARD_REVERSE_ORIENTATION_CST:
      params->paired_end_orientation = FORWARD_REVERSE;
      break;
    case REVERSE_FORWARD_ORIENTATION_CST:
      params->paired_end_orientation = REVERSE_FORWARD;
      break;
    case FORWARD_FORWARD_ORIENTATION_CST:
      params->paired_end_orientation = FORWARD_FORWARD;
      break;
    case MIN_CHIMERA_SCORE_CST:
      params->min_chimera_score = atof(optarg);
      break;
    case '?':
      print_help();
      exit(1);
    default:
      print_help();
      exit(1);
    }
  }

  // Init parameters with values we have fixed with options
  params->init();

  // little hack to get the second file argument of -r option when doing paired-end
  if (optind < argc)
  {
      reads2 = argv[optind++];
  }

  // Init output files if we use --all option
  if(use_all_outputs) {
    nb_streams_to_init = nb_streams;
    // If we are not using paired-end reads we dont create the last stream which is
    // dedicated for paired-end reads only
    if(reads2 == NULL) {
      nb_streams_to_init--;
    }
    max_length = 0;
    for (int i = 0; i < nb_streams_to_init; i++) {
      uint current_length = strlen(extensions[i]);
      if (current_length > max_length)
        max_length = current_length;
    }
    filename = new char[strlen(fileBasename)+max_length+2];
    for (int i=0; i < nb_streams_to_init; i++) {
      sprintf(filename, "%s.%s", fileBasename, extensions[i]);
      streams[i] = create_stream(gzOutput, filename);
//      if (i==16)
//	sam_filename = filename;
    }
    delete [] filename;
  }

  if (indexGenome == NULL  || reads == 0 || 
      threshold == (uint)~0 || (sam_filename.empty() && !use_server)){
    print_help();
    cerr << endl << "Missing at least a mandatory argument among: -i <index> -r <reads> -k <k_length> -o <output_SAM>" << endl << endl;    
    exit(1);
  }else{
    for (int i=0; i < nb_streams && !find_output_file ; i++) {
      if (streams[i] != NULL){
	find_output_file = true;
      }
    }
    if(!sam_filename.empty())
      find_output_file = true;
    if (!find_output_file && !use_server){
      print_help();
      cerr << endl << "Obviously, you must specify at least one output argument" << endl << endl;
      exit(1);
    }
    // check if the there is paired-end input file when using paired-end options
    if (reads2 == NULL && use_paired_end_chimera) {
       print_help();
       cerr << endl << "Please specify a second read input file if you want to use a specific option for paired-end tags." << endl;
    }
  }

  // We open sam_file with HTS lib
  if(!sam_filename.empty()) {
    string mode = "w";
    if(binary_sam)
      mode += "b";
      //mode += "c";
    else if(gzOutput)
      mode += "g";
    if((sam_f = sam_open(sam_filename.c_str(),mode.c_str())) == 0) {
      cerr << "Failed to open " << sam_filename << "for writing" << endl;
      exit(1);
    }
    //hts_set_fai_filename(sam_f, "/data/genomes/GRCh38/GRCh38.fa.fai");
  }

  // Checking fifos for server
  if (use_server) {
    struct stat buf;
    // If no name was specified, use default ones.
    if (! fifo_input_name)
      fifo_input_name = (char *)INPUT_FIFO;
    if (! fifo_output_name)
      fifo_output_name = (char *)OUTPUT_FIFO;

    if (stat(fifo_input_name,&buf) == 0 || stat(fifo_output_name,&buf) == 0) {
      cerr << fifo_input_name << " or " << fifo_output_name << " already exist."
	   << endl << "Please specity another name or delete the existing files." << endl;
      exit(4);
    }
  }

  // Launching exact-matching for reads mapping
  if (use_emt){
    if (reads2 != NULL)
      cerr << "Warning: CRAC--emt version does not support paired-ends protocol, only the first input file "<<reads <<" is considered" << endl;
    if (summary != NULL)
      *summary << "Running CRAC-emt version " <<PACKAGE_VERSION <<" optimized for exact mapping reads " << endl << endl;
    // Indexation of the genome
    if (threshold > 0)
      cerr << "Warning: Running CRAC--emt with k=" << threshold << " which means that it only maps the first located kmer of length "<< threshold <<" for each read" << endl;
    else
      cerr << "Warning: Running CRAC--emt with k=0 to map the entire reads" << endl;
    log = new LocateOnGenome((uchar *)indexGenome); 
    log->setNbLocations(params->max_nb_located_occs);
    log->startChrono();
    LocateTags *lt; 
    lt = new LocateTags(log, params, reads, threshold);
    // put headers in the sam file
    bam_hdr_t* sam_h = writeSamHeader(log,argc,argv,sam_f);
    // launch the locate process
    lt->locate(sam_f,sam_h);
    
    if (summary != NULL) {
      *summary << "Time to locate all tags: "<< log->getElapsedTime() << " s" << endl;       
      *summary << endl << "----------------------------------" << endl;
      *summary << "           Some STATISTICS            " << endl;
      *summary << endl << "----------------------------------" << endl;
      *summary << "Total number of reads analyzed: " << lt->getNbTags() << endl << endl;
      *summary << "Single: " << lt->getNbSingle() << " (" << lt->getNbSingle()*100./lt->getNbTags() << "%)" << endl;
      *summary << "Multiple: " << lt->getNbMultiple() << " (" << lt->getNbMultiple()*100./lt->getNbTags() << "%)" << endl;
      *summary << "None: " << lt->getNbNone() << " (" << lt->getNbNone()*100./lt->getNbTags() << "%)" << endl;
    }
    // Close sam file
    sam_close(sam_f);
    if(sam_h) {
      // FIXME
      //free(sam_h->text);
      //bam_hdr_destroy(sam_h);
    }
    delete lt;
  }
  // Launching CRAC for RNA-Seq mapping 
  else{
    if (! reads_index)
      reads_index = (char *)READS_INDEX;
    // Indexation of reads
    float start = getChrono();

    if(reads2 != 0) {
      reads_reader = new readsReader(reads, reads2, threshold,  params->reads_length);
    } else {
      reads_reader = new readsReader(reads, threshold,  params->reads_length);
    }

    if(strcmp(reads_index,"JELLYFISH") == 0) {
#ifdef HAVE_LIBJELLYFISH
      // Stranded in CRAC means not_canonical in Jellyfish 
      if(reads2 != 0) {
        indexTags = new JellyReadIndex(threshold, 10000, !is_stranded, params->nb_threads, reads, reads2);
      } else {
        indexTags = new JellyReadIndex(threshold, 10000, !is_stranded, params->nb_threads, reads);
      }
#else
      cerr << "JELLYFISH READS INDEX NOT FOUND" << endl;
      exit(1);
#endif
    } else if(strcmp(reads_index,"GKARRAYS") == 0) {
#ifdef HAVE_LIBGKARRAYS
      try {
        if(reads2 != 0) {
          gk = new gkArrays(reads, reads2, threshold, false, params->reads_length, is_stranded, params->nb_threads, params->show_progressbar);
        } else {
          gk = new gkArrays(reads, threshold, false, params->reads_length, is_stranded, params->nb_threads, params->show_progressbar);
        }
      indexTags = new GkReadIndex(gk);
      } catch (bad_alloc e) {
        cerr << "You do not have enough memory to store all the reads. Aborting" << endl;
        exit(42);
      }
#else
      cerr << "GKARRAYS READS INDEX NOT FOUND" << endl;
      exit(1);
#endif
    } else {
      cerr << reads_index << " does not correspond to a read index (please see CRAC help)." << endl;
      exit(1);
    }

    if (summary != NULL){
      if (strcmp(reads_index,"GKARRAYS") == 0){
	if (gk->isLarge()){
	  *summary << "Running CRAC version " <<PACKAGE_VERSION <<" optimized for large data set ";
	  if (params->reads_length){
	    *summary << "and fixed read length " ;
	  }else{
	    *summary << "and variable read length ";
	  }
	}else{
	  *summary << "Running CRAC version " <<PACKAGE_VERSION <<" optimized for small data set < 4Gbp (eg. less than 50M reads of 75 bp) ";
	  if (params->reads_length){
	    *summary << "and fixed read length " ;
	  }else{
	    *summary << "and variable read length " ;
	  }
	}
      }
      if (is_stranded)
	*summary << "for a strand specific library." << endl;
      else
	*summary << "for a non-strand specific library." << endl;
      
      *summary << endl <<  "Time to index all reads: " << (getChrono()-start)/1000000. << " s" << endl;
    }      
    
    // Indexation of the genome
    log = new LocateOnGenome((uchar *)indexGenome); 
    log->setNbLocations(params->max_nb_located_occs);

    // classify the reads
    ct = new ClassifyTags(log, indexTags, reads_reader, params);
    if (use_server) {
      locateServer(ct, fifo_input_name, fifo_output_name);
    } else {
      log->startChrono();
      // put headers in the sam file
      //writeSamHeader(streams[16],log,argc,argv,sam_filename);
      bam_hdr_t* sam_h = writeSamHeader(log,argc,argv,sam_f);
      //*streams[16] << sam_h.c_str();
      //streams[16] << ss.str();
      // put headers with feature for each files
      ct->getHeaders(streams[0], streams[1], streams[2], streams[3],
		     streams[4], streams[5], streams[6], streams[7],
		     streams[8], streams[9], streams[10], streams[11],
		     streams[12], streams[13], streams[14], streams[15],
		     streams[16]);
      // and classify
      ct->classify(streams[0], streams[1], streams[2], streams[3],
		   streams[4], streams[5], streams[6], streams[7],
		   streams[8], streams[9], streams[10], streams[11],
		   streams[12], streams[13], streams[14], streams[15],
		   //streams[16], streams[17], bam, sam_h);
		   streams[16], sam_f, sam_h);

      sam_close(sam_f);
      if(sam_h) {
        // FIXME
        //free(sam_h->text);
        bam_hdr_destroy(sam_h);
      }


      if (summary != NULL){
	uint nb_reads = ct->getNbReads();
	if (reads2 != NULL)
	  nb_reads*=2;
	
	*summary << "Time to classify all reads: "<< log->getElapsedTime() << " s" << endl; 
	
	// making statistics on STDOUT
	*summary << endl << "----------------------------------" << endl;
	*summary << "           Some STATISTICS            " << endl;
	*summary << endl << "----------------------------------" << endl;
	*summary << "Total number of reads analyzed: " << nb_reads << endl << endl;
	*summary << "Single: " << ct->getNbSingle() << " (" << ct->getNbSingle()*100./nb_reads << "%)" << endl;
	*summary << "Multiple: " << ct->getNbMultiple() << " (" << ct->getNbMultiple()*100./nb_reads << "%)" << endl;
	*summary << "None: " << ct->getNbNone() << " (" << ct->getNbNone()*100./nb_reads << "%)" << endl;
	*summary << "Duplication: " << ct->getNbDuplication() << " (" << ct->getNbDuplication()*100./nb_reads << "%)" << endl;
	
	// *summary << endl << "Warning: the sum of the four previous categories may not be equal to 100%." << endl << " This is normal: reads are considered by chunks." << endl << " In a given read, a chunk may appear multiple times, while another just appears once." << endl ; 
	*summary << endl << "----------------------------------" << endl;
	
	*summary << "Explainable: " << ct->getNbExplainable() << " (" << ct->getNbExplainable()*100./nb_reads << "%)" << endl << endl;
	
	*summary << "Repetition: " << ct->getNbRepetition() << " (" << ct->getNbRepetition()*100./nb_reads << "%)" << endl ;
	*summary << "Normal: " << ct->getNbNormal() << " (" << ct->getNbNormal()*100./nb_reads << "%)" << endl ;
	*summary << "Almost-Normal: " << ct->getNbAlmostNormal() << " (" << ct->getNbAlmostNormal()*100./nb_reads << "%)" << endl ;
	*summary << "Sequence-Errors: " << ct->getNbSeqErr() << " (" << ct->getNbSeqErr()*100./nb_reads << "%)" << endl ;
	*summary << "SNV: " << ct->getNbSNP() << " (" << ct->getNbSNP()*100./nb_reads << "%)" << endl ;
	*summary << "Short-Indel: " << ct->getNbBioTagIndel() << " (" << ct->getNbBioTagIndel()*100./nb_reads << "%)" << endl ;
	*summary << "Splice: " << ct->getNbSplice() << " (" << ct->getNbSplice()*100./nb_reads << "%)" << endl ;
	*summary << "Weak-Splice: " << ct->getNbSpliceNoCover() << " (" << ct->getNbSpliceNoCover()*100./nb_reads << "%)" << endl ;
	*summary << "Chimera: " << ct->getNbSpliceInter() << " (" << ct->getNbSpliceInter()*100./nb_reads << "%)" << endl ;
	if(reads2 != NULL) {
	  *summary << "Paired-end Chimera: " << ct->getNbPairedEndChimera() << " (" << ct->getNbPairedEndChimera()*100./nb_reads << "%)" << endl ;
	}
	*summary << "Bio-Undetermined: " << ct->getNbBioUndetermined() << " (" << ct->getNbBioUndetermined()*100./nb_reads << "%)" << endl ;
	*summary << "Undetermined: " << ct->getNbUndetermined() << " (" << ct->getNbUndetermined()*100./nb_reads << "%)" << endl ;
	// *summary << "Nothing: " << ct->getNbNothing() << " (" << ct->getNbNothing()*100./nb_reads << "%)" << endl ;
      }
    }  
    delete ct;
  }
  
  for (int i = 0; i < nb_streams; i++) {
    if (streams[i] != NULL && streams[i] != &cout)
      delete streams[i]; 
  }
  delete gk;
  delete indexTags;
  delete params;
  delete log;
  delete reads_reader;
  if (reads_index && reads_index != (char*) READS_INDEX)
   delete [] reads_index;
  if(fifo_input_name)
   delete fifo_input_name;
  if(fifo_output_name)
   delete fifo_output_name;
  return 0;
}
