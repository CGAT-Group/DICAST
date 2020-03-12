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

#include <algorithm>
#include<iostream>
#include <stdlib.h>
#include <assert.h>
#include "locateOnGenome.h"
#include "utils.h"
using namespace std;

LocateOnGenome::LocateOnGenome(uchar *file_base): out_uniq(&cout), out_multiple(&cout), out_none(&cout), out_sam(&cout) {
  // opening index
  int err = load_index((char *)file_base,(void **)&index);
  if (err) {
    cerr << error_index(err) << endl;
    exit(2);
  }

  // opening conf file
  genomeInfo = new GenomeInfo((char *)file_base);

  nb_tags_treated = 0;
  nb_displayed = 0;
  display_nb = false;
}


LocateOnGenome::~LocateOnGenome(){
  free_index(index);
  delete genomeInfo;
}


ulong LocateOnGenome::getAbsolutePosition(ChrPosition *chrPos) {
  return genomeInfo->getGenomePosition((uchar *)chrPos->getChrPosition(), 
                                       chrPos->getRelativePosition());
}

const uchar *LocateOnGenome::getChromosome(ulong position) {
  return genomeInfo->getChrNameFromPosition(position).first;
}

ulong LocateOnGenome::getIdChromosomeFromName(char *chr_name){
  return genomeInfo->getNumFromName((uchar *)chr_name);
}

ulong LocateOnGenome::getIdChromosome(ulong position) {
  return genomeInfo->getNumFromPosition(position).first;
}

ChrPosition *LocateOnGenome::getChrPos(ulong pos, int strand) {
  if (pos == (ulong)~0)
    return NULL;
  //pair<const uchar *, const ulong> infos = genomeInfo->getChrNameFromPosition(pos);
  //return new ChrPosition(infos.second, (char *)infos.first, strand);
  ChrPosition *chr_pos = genomeInfo->getChrPositionFromPosition(pos);
  chr_pos->setStrand(strand);
  return chr_pos;
}

float LocateOnGenome::getElapsedTime() {
  struct timeval t;
  
  if (gettimeofday(&t,NULL) == -1)
    return 0;
  return (t.tv_sec - chrono_sec) + (t.tv_usec-chrono_usec)/1000000.0;
}

ulong LocateOnGenome::getNbLocations() {
  return nb_displayed;
}

pair<uint, uint>  LocateOnGenome::getFMIndexRange(uchar *tag, ulong tag_length) {
  ulong sp, ep;
  index->getRange(tag, tag_length, sp, ep);
  return pair<uint, uint>(sp, ep);
}

pair<uint, uint> LocateOnGenome::getFMIndexReverseRange(uchar *tag, ulong tag_length) {
  ulong sp, ep;
  index->getReverseCompRange(tag, tag_length, sp, ep);
  return pair<uint, uint>(sp, ep);
}

void LocateOnGenome::getOccurrencesFMIndexRange(ulong sp, ulong ep, ulong **occ, ulong m, ulong nb_disp) {
  if (nb_disp > 0)
    index->returnLocate(sp, ep, nb_disp, occ, m);
  else
    index->returnLocate(sp, ep, nb_displayed, occ, m);
}

GenomeInfo * LocateOnGenome::getGenomeInfo() {
  return genomeInfo;;
}

uchar *LocateOnGenome::getGenomeSubstring(ulong pos, ulong length, int strand, 
                                          char *chr_name) {
  uchar *result;
  uint num_chr;

  // At this point, both length_chr[num_chr - 1] and total_length[nb_chr] are defined (or should be)
  if (chr_name) {
    num_chr = getIdChromosomeFromName(chr_name)+1;
    if (num_chr > genomeInfo->getNbChr()) { // Discard wrong queries
      pos = 0;
      length = 0;
      num_chr = 0;
    }
    if (pos >= genomeInfo->getChrLength(num_chr - 1)) { // Discard wrong queries
      pos = 0;
      length = 0;
    }
    // At this point, we know that either pos < length_chr[num_chr - 1] or pos = 0
    if (pos + length > genomeInfo->getChrLength(num_chr - 1)) { // Adjust erroneous queries
      length = genomeInfo->getChrLength(num_chr - 1) - pos; // No possible unsigned int overflow
    }
    pos = genomeInfo->getGenomePosition(num_chr-1, pos);
  } else {
    if (! genomeInfo->isValidPosition(pos)) { // Discard wrong queries
      pos = 0;
      length = 0;
    }
    // At this point, we know that either pos < total_length[nb_chr] or pos = 0
    if (! genomeInfo->isValidPosition(pos + length - 1)) { // Adjust erroneous queries
      length = genomeInfo->getGenomeLength() - pos; // No possible unsigned int overflow
    }
  }
  if (!length) { // zero length query is an empty sequence
    result = (uchar *)malloc(sizeof(uchar));
    result[0] = 0;
  } else {
    // At this point, we know that 1/ pos + length - 1 doesn't produce unsigned int overflow
    // and that 2/ the subsequence of length 'length' starting at position 'pos' exists.
    index->extract(pos, pos+length-1, &result);
    if (strand == -1) {
      uchar *result2 = (uchar *)malloc((length+1)*sizeof(uchar));
      for (uint i=0; i < length; i++) {
	result2[i] = complementDNA(result[length-i-1]);
      }
      result2[length]=0;
      free(result);
      return result2;
    }
//   else
//     index->extract(pos-length+1, pos, &result);
  }
    return result;
}

void LocateOnGenome::getLocations(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens, ulong **poccs_sens, ulong **poccs_antisens) {
  uchar *up_tag = new uchar[tag_length+1];
  up_tag[tag_length] = 0;
  for (uint i=0; i < tag_length; i++)
    up_tag[i] = toupper(tag[i]);
  nb_sens = index->locate(up_tag, tag_length, poccs_sens,nb_displayed);
  nb_antisens = 0;
  if (nb_sens < nb_displayed) {
    nb_antisens = index->locateReverseComp(up_tag, tag_length, poccs_antisens,nb_displayed-nb_sens);
  }
  delete [] up_tag;
}

pair<ChrPosition **,uint>LocateOnGenome::getLocations(uchar *tag, uint klength, uint tagLength, uint &position_of_locations) {
  ChrPosition **locations=NULL;
  uint nb_locs_multiple=0;
  
  ulong *occs_fwd;
  uint nb_occs_fwd=0;
  ulong *occs_rev;
  uint nb_occs_rev=0;
  pair<uint,uint> ranges_forward;
  pair<uint,uint> ranges_reverse;

  uint i = 0;
  uint nb_occs = 0;
  position_of_locations=0;  
  // searching the first k-mer located
  while (i < (tagLength-klength+1) && nb_occs == 0){
   // multiple is initialized to nb_max_diplayed (--max-locs)
    nb_locs_multiple=getNbLocations(); 
    // get ranges in FM-index  
    ranges_forward = getFMIndexRange((uchar *)&tag[i],klength);
    ranges_reverse = getFMIndexReverseRange((uchar *)&tag[i],klength);

    // First, on forward strand
    if (ranges_forward.first <= ranges_forward.second) {
      getOccurrencesFMIndexRange(ranges_forward.first, ranges_forward.second,&occs_fwd);
      nb_occs_fwd = ranges_forward.second - ranges_forward.first + 1;
    }
    // We check that the number of occurrences is smaller than nb_locs_multiple before the reverse strand treatment
    if (nb_occs_fwd < nb_locs_multiple) {
    if (ranges_reverse.first <= ranges_reverse.second) {
      getOccurrencesFMIndexRange(ranges_reverse.first,
					 ranges_reverse.second,
				 &occs_rev);
      nb_occs_rev = ranges_reverse.second - ranges_reverse.first + 1;
    }
    
    // Again, we check that the total number of occurrences is smaller than nb_locs_multiple
    // and we adjust the nb_locs_multiple if the limit is not reached
    if ((nb_occs_fwd+nb_occs_rev) < nb_locs_multiple)
      nb_locs_multiple=nb_occs_rev+nb_occs_fwd;
    }
    nb_occs = nb_locs_multiple;
    position_of_locations = i;
    i++;
  }
    
  // We save locations
  if (nb_locs_multiple>0){
    locations = (ChrPosition **)malloc(sizeof (ChrPosition *)*nb_locs_multiple);
    
    uint nb_displayed = min(nb_occs_fwd,nb_locs_multiple);
    for (uint i=0 ; i<nb_displayed ; i++){
      locations[i] = getChrPos(occs_fwd[i],1);
    }
    for (uint i=nb_displayed ; i<nb_locs_multiple ; i++){
      locations[i] = getChrPos(occs_rev[i-nb_displayed],-1);
    }
  }
  
  // We delete tempory tables
  if (nb_occs_fwd)
    free(occs_fwd);
  if (nb_occs_rev)
    free(occs_rev);
  
  return pair<ChrPosition **,uint>(locations,nb_locs_multiple);
}

ulong LocateOnGenome::getNbOccsLastLocatedTag() {
  return nbOccs_lastTag;
}

void LocateOnGenome::getNbOccurrences(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens) {
  uchar *up_tag = new uchar[tag_length+1];
  up_tag[tag_length] = 0;
  for (uint i=0; i < tag_length; i++)
    up_tag[i] = toupper(tag[i]);
  nb_sens = index->count(up_tag, tag_length);
  nb_antisens = index->countReverseComp(up_tag, tag_length);
  delete [] up_tag;
}

ulong LocateOnGenome::getNbTagsTreated() {
  return nb_tags_treated;
}


void LocateOnGenome::locateTag(uchar *tag, ulong tag_length) {
  ulong *occ_sens;
  ulong *occ_antisens;
  ulong n_sens, n_antisens;

  getLocations(tag, tag_length, n_sens, n_antisens, &occ_sens, &occ_antisens);

  nbOccs_lastTag = n_sens + n_antisens;
  if (nbOccs_lastTag == 1) {
    displaySingleLocation(tag,tag_length,n_sens, n_antisens, occ_sens, occ_antisens);
  } else if (nbOccs_lastTag > 1) {
    displayMultipleLocations(tag,tag_length,n_sens, n_antisens, occ_sens, occ_antisens);
  } else {
    displayNoLocation(tag, tag_length);
  }
  nb_tags_treated++;
}

void LocateOnGenome::setDisplayNumbers(bool disp) {
  display_nb = disp;
}

void LocateOnGenome::setNbLocations(ulong nb_displayed) {
   this->nb_displayed=nb_displayed;
}

void LocateOnGenome::setOutputStreams (ostream *uniq, ostream *multiple, ostream *none) {
  if (uniq != NULL)
    out_uniq = uniq;
  else
    out_uniq = &cout;
  if (multiple != NULL)
    out_multiple = multiple;
  else
    out_multiple = &cout;
  if (none != NULL)
    out_none = none;
  else
    out_none = &cout;
}
 
void LocateOnGenome::setOutputUniq(ostream *s) {
  out_uniq = s;
}
void LocateOnGenome::setOutputMultiple(ostream *s) {
  out_multiple = s;
}
void LocateOnGenome::setOutputNone(ostream *s) {
  out_none = s;
}

void LocateOnGenome::setOutputSAM(ostream *s) {
  out_sam = s;
}

bool LocateOnGenome::startChrono() {
  struct timeval time;
  
  if (gettimeofday(&time,NULL) == -1)
    return false;
  chrono_sec = time.tv_sec;
  chrono_usec = time.tv_usec;
  return true;
}

void LocateOnGenome::displaySingleLocation(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens) {
  tag[tag_length]=0;
  if (display_nb)
    *out_uniq << getNbTagsTreated() << "|";
  *out_uniq << tag;
  if (n_sens > 0) {
    display_occ(out_uniq, occ_sens, n_sens, tag_length, true);
    free(occ_sens);
  } else {
    display_occ(out_uniq, occ_antisens, n_antisens, tag_length, false);
    free(occ_antisens);
  }
  *out_uniq << endl;
}

void LocateOnGenome::displayMultipleLocations(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens) {
  tag[tag_length]=0;
  if (display_nb)
    *out_multiple << getNbTagsTreated() << "|";
  *out_multiple << tag;
  if (n_sens > 0) {
    display_occ(out_multiple, occ_sens, n_sens, tag_length, true);
    free(occ_sens);
  } 
  if (n_antisens > 0) {
    display_occ(out_multiple, occ_antisens, n_antisens, tag_length, false);
    free(occ_antisens);
  }
  *out_multiple << endl;
}

void LocateOnGenome::displayNoLocation(uchar *tag, ulong tag_length) {
  tag[tag_length]=0;
  if (display_nb)
    *out_none << getNbTagsTreated() << "|";
  *out_none << tag << endl;
}

void LocateOnGenome::display_occ(ostream *stream, ulong *occ, ulong n, ulong tag_length, bool sens) {
  ulong num_chr = 0;
  ulong old_num_chr = genomeInfo->getNbChr() + 1;
  int strand = (sens) ? 1 : -1;
  sort(occ,occ+n);
  for (ulong i=0; i < n; i++) {
    pair<const ulong, const ulong> info_pos = genomeInfo->getNumFromPosition(occ[i]);
    num_chr = info_pos.first;
    if (num_chr != old_num_chr) {
      old_num_chr = num_chr;
      *stream << "|"<< genomeInfo->getChrName(num_chr) << ","<< strand <<"\t";
    }else{
      *stream << "," ;
    }
    // Positioning according to the chromosomes
    *stream << info_pos.second ;
  }
}
