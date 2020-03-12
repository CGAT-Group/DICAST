/******************************************************************************
*                                                                             *
*  Copyright © 2010-2013 -- IRB/INSERM                                        *
*                           (Institut de Recherches en Biothérapie /          *
*                           Institut National de la Santé et de la Recherche  *
*                           Médicale)                                         *
*                           LIFL/INRIA                                        *
*                           (Laboratoire d'Informatique Fondamentale de       *
*                           Lille / Institut National de Recherche en         *
*                           Informatique et Automatique)                      *
*                           LIRMM/CNRS                                        *
*                           (Laboratoire d'Informatique, de Robotique et de   *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique)     *
*                           LITIS                                             *
*                           (Laboratoire d'Informatique, du Traitement de     *
*                           l'Information et des Systèmes).                   *
*                                                                             *
*                                                                             *
*  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*                   Thierry LECROQ   <thierry.lecroq@univ-rouen.fr>           *
*                   Martine LÉONARD  <Martine.Leonard@univ-rouen.fr>          *
*                   Éric RIVALS      <eric.rivals@lirmm.fr>                   *
*                                                                             *
*  Programmeurs                                                               *
*      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*                   Jérôme AUDOUX    <jerome.audoux@etud.univ-montp2.fr>      *
*  with additional contribution for the packaging of:	                      *
*                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               *
*                                                                             *
*  Contact:         Gk-Arrays list   <crac-gkarrays@lists.gforge.inria.fr>    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce fichier fait partie de la librairie Gk-arrays.                          *
*                                                                             *
*  La librairie Gk-arrays  a  pour objectif d'indexer de grands ensembles de  *
*  lectures de séquences issues du séquençage haut-débit.                     *
*                                                                             *
*  Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  *
*  développement  et à la reproduction du  logiciel par  l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  *
*  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  *
*  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  *
*  pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  *
*  termes.                                                                    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This File is part of the Gk-arrays library.                                *
*                                                                             *
*  The Gk-arrays library aims at indexing k-factors from a huge set of        *
*  sequencing reads.                                                          *
*                                                                             *
*  This software is governed by the CeCILL-C license under French law and     *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL-C   *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL-C license and that you accept its terms.           *
*                                                                             *
******************************************************************************/
 
#ifndef GKARRAYS_H
#define GKARRAYS_H

#include <utility>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include "gkArraysTypes.h"
#include "solArray.h"

#include "sux/rank9b.h"
#include "readsReader.h"

namespace gkarrays {

  extern "C" {
    /* One can use this function to test the library availability */
    char *libGkArraysVersion();
  }

  class gkArrays {
  private:
    uchar *alltags; 
    uintSA nb_tags;
    SolArray *gkSA;
    SolArray *gkISA;
    SolArray *gkCFPS;
 
    uint threshold;
    bool use_bitvector;
    uint tag_length;
    uint64_t *read_lengths;
    rank9b *read_lengths_select;
    rank9b *read_lengths_rank;

    readsReader *reads;

    bool is_large;
    bool is_stranded;
    uint nb_threads;

    /**
     * A static variable to control whether progress bar should be displayed or not...
     */
    static bool show_progressbar;

    /**
     * A static variable to store the number of columns for displkay the progress bar...
     */
    static unsigned int nb_cols;

  public:
    /**
     * Construct the read index
     *
     * @param tags_file Name of the file containg the reads
     * @param threshold length of k-mers we have to use
     * @param use_bitvector: true iff we must store the array using a bit vector
     *        (slower but more space efficient)
     * @param tag_length length of the reads. If a shorter read is found, 
     *             it raises
     *             an error. If a longer read is found, only the prefix
     *             of tag_length characters is kept.
     *             If tag_length == 0 (default), just gess what the read length
     *             is.
     * @param stranded: true iff we know which strand has been sequenced and,
     *        therefore, (for instance) AACG must not be considered as equal
     *        to its revcomp (CGTT).
     * @param nb_threads allows to build GkSA on a multi-thread architecture
     * @param show_progressbar display a progress bar when building the index
     */
    gkArrays(char *tags_file, uint threshold, bool use_bitvector=false, uint tag_length=0, bool stranded=false,uint nb_threads=1, bool show_progressbar = false);

    /**
     * Alternative to construct the read index with paired-end reads
     *
     * @param tags_file1 Name of the file containing the reads of the first pair
     * @param tags_file2 Name of the file containing the reads of the second pair
     * @param threshold length of k-mers we have to use
     * @param use_bitvector: true iff we must store the array using a bit vector
     *        (slower but more space efficient)
     * @param tag_length length of the reads. If a shorter read is found, 
     *             it raises
     *             an error. If a longer read is found, only the prefix
     *             of tag_length characters is kept.
     *             If tag_length == 0 (default), just gess what the read length
     *             is.
     * @param stranded: true iff we know which strand has been sequenced and,
     *        therefore, (for instance) AACG must not be considered as equal
     *        to its revcomp (CGTT).
     * @param nb_threads allows to build GkSA on a multi-thread architecture
     * @param show_progressbar display a progress bar when building the index
     */
    gkArrays(char *tags_file1, char *tags_file2, uint threshold, bool use_bitvector=false, uint tag_length=0, bool stranded=false,uint nb_threads=1, bool show_progressbar = false);

    ~gkArrays();

    /**
     * @return true iff the read is not suitable ie. if it is shorter than
     *         the specified length (if any) or shorter than the specified
     *         k-mer length.
     */
    static bool isDiscarded(uint actual_length, uint theoretical_length=0, uint k=0);

    /**
     * Convert a position from P-position to Q-position
     * (if you do not understand this, please read our article!).
     * That converts a position as in the concatenation of reads
     * to the position in GkIFA (for example).
     * In the article, values of GkSA are also renumbered to Q-position
     * but we do not renumber them in practice (it is quite useless).
     * @param i: a P-position
     * @return a Q-position
     */
    uintSA convertPposToQpos(uintSA i);

    /**
     * Gives the end position of a given read in the concatenation of reads.
     * @param tag_num: tag number
     * @return the end position of the read #tag_num in C_R (the concatenation
     *         of reads)
     */
    uintSA getEndPosOfTagNum(uint tag_num);

    /**
     * @param i the index position in the array (starting at 0).
     * @return the value of GkCFA at the given index ie. the number of k-factors
     * of rank i, where i is the requested index.
     */
    uintSA getGkCFA(uintSA i);

    /**
     * @return the number of elements in the GkCFA array. In other terms it
     * corresponds to the number of distinct k-mers in the input.
     */
    uintSA getGkCFALength();

    /**
     * @param i the index position in the array (starting at 0).
     * @return the value of GkISA at the given index ie. the rank of the k-factor
     * at position P-position i.
     */
    uintSA getGkISA(uintSA i);
    /**
     * @param i the index position in the array (starting at 0).
     * @return the value of GkSA at the given index ie. the P-position of the 
     * k-factor whose rank is i
     */
    uintSA getGkSA(uintSA i);

    /**
     * @return the number of entries in gkSA (ie. the number of P-positions)
     */
    uintSA getGkSALength();

    /**
     * @return the number of P-positions in Cr from a number of reads (fixed length or not)
     * This function is available before the construction of gkSA.
     */
    uintSA getNbPposition(uintSA nb_reads);

    /**
     * @return the number of tags (or reads) indexed in the Gk Arrays
     */
    uint getNbTags();

    /**
     * @param tag_num The number of the tag in the index
     * @param pos_factor Position of the factor in the tag
     * @param multiplicity Counts (if false) only once a tag
     *                      that contains the factor many times
     * @return Return the number of tags sharing the factor
     *         starting at position pos_factor in the tag
     *         tag_num.
     *         This is the number of elements returned
     *         by the function getTagsWithFactor(.)
     */
    uint getNbTagsWithFactor(uint tag_num, uint pos_factor, bool multiplicity=0);

    /**
     * @return the number of threads the GkArrays have been told to use.
     *         The threads can be used for the construction.
     */
    uint getNbThreads();

    /**
     * @param i The number of the tag in the index
     * @return  the tag number of the paired-end read associated with i 
     *          or -1 if reads are not paired-end.
     */
    uint getPair(uint i);

    /**
     * @return the rank of the P-k factor starting at position pos_factor in the 
     * read number tag_num.
     */
    uintSA getPosInCommon(uint tag_num,uint pos_factor); 
    
    /**
     * @return the object that allows to get a readIterator
     */
    readsReader *getReads();
  
    /**
     * Gives the start position of a given read in the concatenation of reads.
     * @param tag_num: tag number
     * @return the start position of the read #tag_num in C_R (the concatenation
     *         of reads)
     */
    uintSA getStartPosOfTagNum(uint tag_num);

    /**
     * Gives the start Q-position of a given read in the ISA array
     * @param tag_num: tag number
     * @return the start Q-position of the read #tag_num in GkISA.
     */
    uintSA getStartQPosOfTagNum(uint tag_num);
  
    /**
     * @param i Tag number
     * @return an array whose length is getSupportLength(i)
     *         and where the value at position k is the number of occurrences
     *         of the k-factor starting at position k in the reads
     *         among all the Pk-factors.
     */
    uint *getSupport(uint i);

    /**
     * Return the length of the support.
     * @return getTagLength(i) - getThreshold()+1
     */
    uint getSupportLength(uint i=0);

    /**
     * @param i the read number to be retrieved
     * @return the read number i.
     */
    char *getTag(uint i);

    /**
     * @param i Tag number (if the length is not constant)
     * @return the length of the read.
     */
    uint getTagLength(uint i=0);

    /**
     * @param i The number of the tag in the index
     * @param p Position of the factor in the tag
     * @param l The length of the factor
     * @return the factor at the position p in the tag number i
     */
    char *getTagFactor(uint i, uint p, uint l);

    /**
     * Gives the number of a read
     * @param pos a position in SA or in the concatenated
     *             sequence of reads
     * @return the read number where this position lies
     */
    uint getTagNum(uintSA pos);

    /**
     * Return the number of tag and the relative position in that tag
     * corresponding to a given position in the concatenation of reads
     * @param pos position in the concatenation of reads
     * @return a pair whose fist element is the tag number and the second 
     *         element is the position in the tag.
     */
    std::pair<uint, uint> getTagNumAndPosFromAbsolutePos(uintSA pos);

    /**
     * @param tag_num The number of the tag in the index
     * @param pos_factor Position of the factor in the tag
     * @return Return an array that contains each tag number
     *         where the factors matches.
     * @post The array is sorted
     */
    uint *getTagNumWithFactor(uint tag_num, uint pos_factor);

    /**
     * @param tag_num The number of the tag in the index
     * @param pos_factor Position of the factor in the tag
     * @return Return an array composed of pairs 
     *         (tag, pos) corresponding to all the Pk-factors
     *         equal to the Pk-factor 
     *         starting at position pos_factor in the tag
     *         tag_num.
     * @post The array is sorted according to read number and read position
     */
    std::pair <uint, uint> *getTagsWithFactor(uint tag_num, uint pos_factor);

    /**
     * @param factor the pattern to be searched.
     * @param factor_length the length of the factor, should be <= getThreshold()
     * @param nb_fact nb_fact is used to give the number of occurrences
     *                 in the array.
     * @return Return an array composed of pairs 
     *         (tag, pos) corresponding to all the Pk-factors
     *         equal to the k-factor factor
     */
    std::pair <uint, uint> *getTagsWithFactor(char *factor, uint factor_length,
					      uint &nb_fact);

    /**
     * @param pos The position from where we want to retrieve a text subtring.
     * The position must be given in the original text (not the filtered one).
     * @param length the length of the substring to be retrieved.
     * @return text factor at position pos of length length.
     *         The returned string is NULL-terminated.
     */
    char *getTextFactor(uintSA pos, uint length);

    /**
     * @return return the length of the k-factors (ie. return k).
     */
    uint getThreshold();
  
    /**
     *@return the array type used for building GkSA and GkISA (either
     *        SMALL_ARRAY, LARGE_ARRAY or OPTIMAL_ARRAY).
     */
    array_type getType();

    /**
     *@return true if the nbPposition > 2^32 
     */
    bool isLarge();

    /**
     * @return true iff the position does not lie in the threshold - 1 
     * last characters of a read, ie. if it is a P-position.
     */
    bool isPposition(uintSA pos);

    /**
     * @return true iff the GkArrays have been built as a strand-dependant
     *         index. Therefore a k-mer and its revcomp won't be considered
     *         as equal.
     */
    bool isStranded();

    /**
     * @param i the number of the tag in the index
     * @return true if the tag is the first member of is pair
     *         in case of paired-end files. False either
     */
    bool isTheFirstMemberOfPair(uint i);

    /**
     * This static method return the value of the static show_progressbar boolean value.
     */
    static bool ShowProgressBar();

    /**
     * This static method return the value of the static nb_cols.
     */
    static unsigned int GetNbColumns();

    private:

    /**
     * The function init is called at the end of each constructor
     * Init() build the gKarrays and related structures
     */
    void init();

  };

}
#endif

// Local Variables:
// mode:c++
// End:
