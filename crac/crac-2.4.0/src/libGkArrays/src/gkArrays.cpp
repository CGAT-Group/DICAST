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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include "basic.h"
#include "utils.h"
#include "gkArrays.h"
#include "gkArraysBuilder.h"
#include "const.h"
#include <sys/ioctl.h>
#ifdef HAVE_LIBPROGRESSBAR
#include <libProgressBar/progressBar.h>
#include <unistd.h>
#endif

// extern "C" {
// #include <divsufsort.h>
// #include <divsufsort64.h>
// }


using namespace std;

namespace gkarrays {

  extern "C" {
    char *libGkArraysVersion() {
      return (char*)PACKAGE_VERSION;
    }
  }

  ////////////////////////////////////////////////////////////

  bool gkArrays::show_progressbar = false;
  unsigned int gkArrays::GetNbColumns() {
    return nb_cols;
  }

  unsigned int gkArrays::nb_cols = 80;
  bool gkArrays::ShowProgressBar() {
    return show_progressbar;
  }

  ////////////////////////////////////////////////////////////

  bool gkArrays::isDiscarded(uint actual_length, uint theoretical_length, uint k) {
    return (actual_length < theoretical_length || actual_length < k);
  }

  gkArrays::gkArrays(char *tags_file, uint threshold, bool use_bitvector,
                     uint tag_length, bool stranded, uint nb_threads,
		     bool show_progressbar)
    :threshold(threshold), use_bitvector(use_bitvector), tag_length(tag_length), 
     is_stranded(stranded), nb_threads(nb_threads)
  {
    gkArrays::show_progressbar = show_progressbar;
    reads = new readsReader(tags_file,threshold,tag_length);
    init();  
  }

  gkArrays::gkArrays(char *tags_file1, char *tags_file2, uint threshold, bool use_bitvector,
                     uint tag_length, bool stranded, uint nb_threads,
		     bool show_progressbar)
    :threshold(threshold), use_bitvector(use_bitvector), tag_length(tag_length), 
     is_stranded(stranded), nb_threads(nb_threads)
  {
    gkArrays::show_progressbar = show_progressbar;
    reads = new readsReader(tags_file1,tags_file2,threshold,tag_length);
    init();  
  }

  void gkArrays::init() {
    if (threshold > tag_length && tag_length){
      cerr << "k-mer length is greater than read length" << endl;
      exit(9);
    }
    uintSA totalLength = 0;
    is_large = false;
    nb_tags = 0;
    
    // Reading reads to know the total length we have.
    readIterator *it = reads->begin(true);
    while (!it->isFinished()) {
      totalLength += (tag_length > 0 ? tag_length : it->getLength());
      nb_tags++;
      ++(*it);
    }

    delete it;

    this->tag_length = tag_length;
  
    // Initialise the mask which stores the read lengths
    if (!tag_length){
      read_lengths = (uint64_t*)Calloc<uint64_t>(nb_digits(totalLength+2, 1));
      // Add an extra 1 at the beginning
      bitset(read_lengths,0);
      // And we also have an extra 0 at the end.
    }    

    if (nb_tags == 0) {
      cerr << "There is no read in the input file!" << endl;
      exit(12);
    }

    alltags = Calloc<uchar>(totalLength + 1);
    uintSA currentPos = 0;
    it = reads->begin();

#ifdef HAVE_LIBPROGRESSBAR
    DoccY::ProgressBar *PB = NULL;
    unsigned int pb_step = 0;
    unsigned int pb_refresh = 0;
    if (show_progressbar) {
#ifdef TIOCGSIZE
      struct ttysize ts;
      nb_cols = ioctl(STDIN_FILENO, TIOCGSIZE, &ts) ? 80 : ts.ts_cols;    
#elif defined(TIOCGWINSZ)
      struct winsize ts;
      nb_cols = ioctl(STDIN_FILENO, TIOCGWINSZ, &ts) ? 80 : ts.ws_col;
#endif /* TIOCGSIZE */

      PB = new DoccY::ProgressBar ("Scanning Reads  ", nb_tags, nb_cols, cerr, false);
      PB->ShowPercent();
      PB->ShowTime();
      PB->update();
      pb_refresh = nb_tags / (nb_cols ? nb_cols : 100);
      pb_refresh = pb_refresh > 0 ? pb_refresh : 1;
    }
#endif
    // strore each tags in the char array alltags[]
    while (!it->isFinished()) {
      if (tag_length) {
	if (it->getLength() >= tag_length) {
	  strncpy((char *) &(alltags[currentPos]), it->getSequence(),
		  tag_length);
	  currentPos += tag_length;
	}
      } else {
	if (threshold <= it->getLength()) {
	  strncpy((char *) &(alltags[currentPos]), it->getSequence(),
		  it->getLength());
	  currentPos += it->getLength();
	  bitset(read_lengths, currentPos);
	  // read_lengths[currentPos] = 1;
	}
      }
      ++(*it);

#ifdef HAVE_LIBPROGRESSBAR
      if (show_progressbar) {
	PB->Step(!(++pb_step % pb_refresh));
      }
#endif
    }
#ifdef HAVE_LIBPROGRESSBAR
    if (show_progressbar) {
      PB->update(false);
      delete PB;
    }
#endif

    delete it;
    dnaFiltration((char *)alltags, totalLength);

    if (!tag_length){
      // Create the rank and select structure on the bit vector
      read_lengths_rank = new rank9b(read_lengths, totalLength+2);
      read_lengths_select = read_lengths_rank;
    }


    // CR with DNA 2 bits coding 
    for (uintSA i=0; i < totalLength/4; i+=1) {
      alltags[i] = (uchar)DNAtoInt((char *)&alltags[i*4],4);
    }
    
    alltags[totalLength/4] = (uchar)(DNAtoInt((char *)&alltags[(totalLength/4)*4],
					      (totalLength)%4) 
				     << (8 - 2*(totalLength%4)));
    
    alltags = Realloc<uchar>(alltags, totalLength/4+1);

    // Determine what type of arrays we are going to store
    is_large = (totalLength >= (uint) numeric_limits<uint>::max());
    array_type type_to_build = (is_large ? LARGE_ARRAY : SMALL_ARRAY);
    if (use_bitvector)
      type_to_build = OPTIMAL_ARRAY;

    // Creates the structure
    GkArraysBuilder gkBuild(alltags,totalLength,threshold,this,type_to_build,nb_threads);
    
    gkSA = gkBuild.getGkSA();
    gkISA = gkBuild.getGkIFA();
    gkCFPS = gkBuild.getGkCFPS();
  }
  
  gkArrays::~gkArrays() {
    free(alltags);
    if (!tag_length){
      delete read_lengths_select;
      // delete read_lengths_rank;
      free(read_lengths);
    }
    delete gkSA;
    delete gkISA;
    delete gkCFPS;
    delete reads;
  }

  uintSA gkArrays::getEndPosOfTagNum(uint tag_num) {
    return getStartPosOfTagNum(tag_num) + getTagLength(tag_num) - 1;
  }
  
  uintSA gkArrays::getGkCFA(uintSA i) {
    return (*gkCFPS)[i+1] - (*gkCFPS)[i];
  }

  uintSA gkArrays::getGkCFALength() {
    return gkCFPS->length()-1;
  }
  
  uintSA gkArrays::getGkISA(uintSA i) {
    return (*gkISA)[i];
  }

  uintSA gkArrays::getGkSA(uintSA i) {
    return (*gkSA)[i];
  }

  uintSA gkArrays::getGkSALength() {
    return gkSA->length();
  }
  
  uintSA gkArrays::getNbPposition(uintSA nb_reads){
    return getStartQPosOfTagNum(nb_reads-1)+getTagLength(nb_reads-1)-getThreshold()+1;
  }
  
  uint gkArrays::getNbTags(){
    return nb_tags;
  }

  uint gkArrays::getNbTagsWithFactor(uint tag_num, uint pos_factor, 
				     bool multiplicity){
    uintSA pos_in_common = getPosInCommon(tag_num,pos_factor);
    uint nb;
    nb =(*gkCFPS)[pos_in_common+1] - (*gkCFPS)[pos_in_common];
    if (multiplicity || nb==1) {
      return nb;
    } else {
      nb = 1;
      for (uintSA i=(*gkCFPS)[pos_in_common]+1;
           i < (*gkCFPS)[pos_in_common+1]; i++) {
        // Do the two positions belong to the same read?
        if (getTagNum((*gkSA)[i]) != getTagNum((*gkSA)[i-1]))
          nb++;
      }
    }
    return nb;
  }

  uint gkArrays::getPair(uint i){
    if(reads->isPairedEnd()){
      if(i%2 == 0) return i+1;
      else return i-1;
    } else {
      return -1;
    }
  }

  bool gkArrays::isTheFirstMemberOfPair(uint i){
    if(reads->isPairedEnd()){
      if(i%2 == 0) return true;
      else return false;
    } else {
      cerr << "Ask for a paired-end function using non-paired-end reads" << endl;
      return false;
    }
  }

  readsReader *gkArrays::getReads() {
    return reads;
  }

  uintSA gkArrays::getStartPosOfTagNum(uint tag_num) {
    if (tag_length){
      return (uintSA) tag_num * tag_length;
    }else{
      // +1 when starting at 1
      return read_lengths_select->select(tag_num);
    }  
  }
  
  uintSA gkArrays::getStartQPosOfTagNum(uint tag_num) {
    if (tag_length){ 
      return (uintSA)tag_num * tag_length - tag_num * (threshold - 1);
    }else{
      // tag_num+1 when starting at 1
      return read_lengths_select->select(tag_num) - tag_num * (threshold - 1);
    }
  }

  char *gkArrays::getTag(uint i) {
    uint this_tag_length = getTagLength(i);
    uint length_char = this_tag_length+1;
    char dna[5];
    char *final_dna = new char[length_char];
    uint previous_length;
    uintSA readPos = getStartPosOfTagNum(i);
    uintSA nextReadPos = readPos + getTagLength(i);
    dna[4]=0;
    final_dna[length_char-1] = 0;
    intToDNA(alltags[readPos/4],4,dna);
    previous_length = 4-(readPos%4);
    strncpy(final_dna, &dna[readPos%4], previous_length);

    uint j;
    for (j=1 ; previous_length+4 <= this_tag_length ; j++) {
      intToDNA(alltags[readPos/4+j],4,dna);
      strncpy(&final_dna[previous_length], dna, 4);
      previous_length+=4;
    }
    if (nextReadPos%4 > 0) {
      intToDNA(alltags[nextReadPos/4], 4, dna);
      strncpy(&final_dna[previous_length], dna, nextReadPos%4);
    }
    return final_dna;
  }

  char *gkArrays::getTagFactor(uint i, uint p, uint fact_length) {
    uint length_char = fact_length+1;
    char dna[5];
    char *final_dna = new char[length_char];
    uint previous_length = 0;
    uintSA readPos = getStartPosOfTagNum(i);
    dna[4]=0;
    final_dna[length_char-1] = 0;

    // Retrieving text factors that end a compacted int
    intToDNA(alltags[(p+readPos)/4],4,dna);
    // Either we go to the end of the int or we stop before if we want to retrieve a short factor
    previous_length = (fact_length < 4-((p+readPos)%4) ? fact_length : 4-((p+readPos)%4));
    strncpy(final_dna, &dna[(p+readPos)%4], previous_length);

    // Retrieving text factors that are fully contained in a compacted int
    uint j = 0;
    for (j=1 ; previous_length+4 <= fact_length ; j++) {
      intToDNA(alltags[(p+readPos)/4+j],4,dna);
      strncpy(&final_dna[previous_length], dna, 4);
      previous_length+=4;
    }
    // Retrieving factors that are prefix in a compacted int
    if ((p+fact_length+readPos)%4 > 0 && previous_length < fact_length) {
      intToDNA(alltags[(p+readPos)/4+j], 4, dna);
      strncpy(&final_dna[previous_length], dna, (p+fact_length+readPos)%4);
    }
    return final_dna;
  }

  uint gkArrays::getTagNum(uintSA pos) {
    if (tag_length){
      return pos/tag_length;
    }else{
      return read_lengths_rank->rank(pos+1) - 1;
    }
  }

  std::pair<uint, uint> gkArrays::getTagNumAndPosFromAbsolutePos(uintSA pos) {
    std::pair<uint, uint> result;
    result.first = getTagNum(pos);
    result.second = pos - getStartPosOfTagNum(result.first);
    return result;
  }
  
  uint *gkArrays::getTagNumWithFactor(uint tag_num, uint pos_factor) {
    uintSA pos_in_common = getPosInCommon(tag_num,pos_factor);
    uintSA nb_init;
    uint *result;
    uint nb = 1;

    nb_init = (*gkCFPS)[pos_in_common+1] - (*gkCFPS)[pos_in_common];
    result = Calloc<uint>(nb_init);
    result[0] = getTagNum((*gkSA)[(*gkCFPS)[pos_in_common]]);
    for (uintSA i=(*gkCFPS)[pos_in_common]+1;
         i < (*gkCFPS)[pos_in_common+1]; i++) {
      uint tag_num = getTagNum((*gkSA)[i]);
      if (tag_num != result[nb-1]) {
        result[nb] = tag_num;
        nb++;
      }
    }
    if (nb != nb_init) {
      return Realloc(result, nb);
    }
    return result;    
  }

  pair <uint, uint> *gkArrays::getTagsWithFactor(uint tag_num, uint pos_factor){
    uintSA pos_in_common = getPosInCommon(tag_num,pos_factor);
    uint nb_fact = getNbTagsWithFactor(tag_num,pos_factor, true);
    pair <uint, uint> *tagsAndPos = new pair<uint,uint>[nb_fact];
    uintSA nbc = (*gkCFPS)[pos_in_common];
    for (uint i = 0 ; i < nb_fact ; i++){
      if (tag_length){
        tagsAndPos[i] = pair<uint,uint>((*gkSA)[nbc+i]/tag_length,
                                        (*gkSA)[nbc+i]%tag_length);
      }else{
        tagsAndPos[i] = getTagNumAndPosFromAbsolutePos((*gkSA)[nbc+i]);
      }
    }
    return tagsAndPos;
  }
  
  pair <uint, uint> *gkArrays::getTagsWithFactor(char *factor, uint factor_length
						 , uint &nb_fact) {
    // Binary search over the suffix array to find the occurrences
    // Code from Gonzalez and Makinen
    uintSA left=0;
    uintSA right = gkSA->length()-1;
    uintSA currentIndex=(left+right)/2;
    uintSA min=0;
    uintSA max=0;
    char *text_factor;
    /* Search for the left boundary */
    while (left < right-1) {
      text_factor = getTextFactor((*gkSA)[currentIndex], factor_length);
      if (strcmp(text_factor, factor) >= 0) {
	right = currentIndex;
      } else {
	left = currentIndex;
      }
      currentIndex = (right+left)/2;
      delete [] text_factor;
    }
    min=right;

    /* Search for the right boundary */
    left=min-1;
    right=gkSA->length()-1;
    currentIndex=(left+right)/2;
    while (left < right-1) {
      text_factor = getTextFactor((*gkSA)[currentIndex], factor_length);
      if (strcmp(text_factor, factor) <= 0) {
	left = currentIndex;
      } else {
	right = currentIndex;
      }
      currentIndex = (right+left)/2;
      delete [] text_factor;
    }
    max=left;

    nb_fact = max-min+1;
    pair<uint, uint> *occurrences = new pair<uint, uint>[nb_fact];
    for (uint i=0; i < nb_fact; i++) {
      occurrences[i] = getTagNumAndPosFromAbsolutePos((*gkSA)[min+i]);
    }

    return occurrences;
  }

  uint *gkArrays::getSupport(uint i) {
    uint this_tag_length = getTagLength(i);
    uint posMax_tag = this_tag_length-threshold+1;
    uint *support = new uint[posMax_tag];
    uintSA startPos_tag = getStartQPosOfTagNum(i);
    for (uint j=0; j < posMax_tag; j++) {
      support[j] = (*gkCFPS)[(*gkISA)[j+startPos_tag]+1] 
        - (*gkCFPS)[(*gkISA)[j+startPos_tag]];
    }

    return support;
  }

  uint gkArrays::getTagLength(uint i){
    if (tag_length){  
      return tag_length;
    }else{
      // i+2 and i+1 when starting at 1
      return read_lengths_select->select(i+1) 
	- read_lengths_select->select(i);
    }  
  }

  char *gkArrays::getTextFactor(uintSA pos, uint length) {
    std::pair<uint,uint> infos = getTagNumAndPosFromAbsolutePos(pos);
    return getTagFactor(infos.first,infos.second, length);
  }

  uint gkArrays::getThreshold(){
    return threshold;
  }

  uint gkArrays::getSupportLength(uint i) {
    return getTagLength(i) - threshold +1;
  }

  array_type gkArrays::getType(){
    return gkSA->getType();
  }

  bool gkArrays::isLarge(){
    return is_large;
  }

  bool gkArrays::isPposition(uintSA pos) {
    if (tag_length){  
      return pos % tag_length <= tag_length - threshold;
    }else{
      uint tag_num = getTagNum(pos);
      // tag_num+2 when starting at 1
      uintSA end_read = read_lengths_select->select(tag_num + 1);
      return end_read  - pos >= threshold;
    }
  }
  
  bool gkArrays::isStranded(){
    return is_stranded;
  }

  uintSA gkArrays::convertPposToQpos(uintSA i) {
    return i - getTagNum(i) * (threshold - 1);
  }

  uintSA gkArrays::getPosInCommon(uint tag_num, uint pos_factor){
    uintSA start_pos = getStartQPosOfTagNum(tag_num);
    return (*gkISA)[start_pos + pos_factor];
  }
}

