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

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sys/types.h>
#include <sys/times.h>
#include <cassert>
#include <cstring>
#include <stdint.h>
#include "SSA.h"
#include "utils.h"
//---------------------------------------------------------------------------

/* Three function to variables to manage parameters */
static bool is_delimeter(char *delimiters, char c) {
  int i=0,len_delimiters=strlen(delimiters);
  bool is=false;
  for (i=0;i<len_delimiters;i++)
    if (c == delimiters[i]) is=true;
  return is;
}

static void parse_parameters(char *options, int *num_parameters, char ***parameters, char *delimiters) {
  int i=0,j=0,temp=0,num=0, len_options=strlen(options);
  char *options_temp;
  while  (i<len_options) {
    while ((i<len_options) && is_delimeter(delimiters,options[i])) i++;
    temp=i;
    while ((i<len_options) && !is_delimeter(delimiters,options[i])) i++;
    if (i!=temp) num++;
  }
  (*parameters) = (char **) malloc(num*sizeof(char * ));
  i=0;
  while  (i<len_options) {
    while ((i<len_options) && is_delimeter(delimiters,options[i])) i++;
    temp=i;
    while ((i<len_options) && !is_delimeter(delimiters,options[i])) i++;
    if (i!=temp) {
      (*parameters)[j]=(char *) malloc((i-temp+1)*sizeof(char));
      options_temp = options+temp;
      strncpy((*parameters)[j], options_temp, i-temp);
      ((*parameters)[j])[i-temp] = '\0';
      j++;
    }
  }
  *num_parameters = num;
}

static void free_parameters(int num_parameters,char ***parameters) {
  int i=0;
  for (i=0; i<num_parameters;i++)
    free((*parameters)[i]);
  free((*parameters));
}

void swap(ulong *pos, ulong i, ulong j) {
  ulong tmp = pos[i];
  pos[i] = pos[j];
  pos[j] = tmp;
}

int compare (const void * a, const void * b)
{
  if (*(ulong *)a < *(ulong *)b) 
    return -1;
  else if (*(ulong *)a == *(ulong *)b)
    return 0;
  else
    return 1;
}


void sortPositions(ulong *pos, ulong *suffixes, ulong length, ulong mult) {
  memcpy(suffixes, pos, length*sizeof(ulong));
  qsort(suffixes, length, sizeof(ulong), compare);

  ulong min, max;
  ulong index;
  bool found;
  uint *ranks = new uint[length];
  for (ulong i =0; i < length ;  i++) {
    min = 0; 
    max = length-1;
    found = false;
    //    cerr << i << " (searching "<<pos[i]<<")"<< endl;
    while (! found) {
      index = (min+max)/2;
      if (pos[i] < suffixes[index]) {
        max = index-1;
      } else if (pos[i] > suffixes[index]) {
        min = index+1;
      } else {
        found = true;
      }
    }
    ranks[index] = i;
  }

  for (ulong i = 0; i < length; i++) {
    suffixes[i] = ranks[i]*mult;
  }
  delete [] ranks;


}

/**********************************************************************************************************
 * succinct full-text index                                                                               *
 **********************************************************************************************************/


TFMindex::TFMindex(uchar *bwt, ulong length, ulong samplerate, bool free_text,ulong factor, int *error, uchar terminator) { 
  this->n = length;
  this->samplerate = samplerate;
  //        uchar *bwt = BWT(text,free_text,factor);
       
  if (bwt == NULL) *error = 1;
  uint terminator_pos = 0;
  map0 = terminator;
  // caller can delete text now
  ulong i;
  for (i=0;i<size_uchar+1;i++)
    C[i]=0;
  for (i=0;i<n;++i) {
    if (bwt[i] == terminator) {
      terminator_pos = i;
    } 
    C[bwt[i]]++;
  }
  ulong prev=C[0], temp;
  C[0]=0;
  for (i=1;i<size_uchar+1;i++) {
    temp = C[i];
    C[i]=C[i-1]+prev;
    prev = temp;
  }
  codetable = makecodetable(bwt,n);

       
  //        for (uint i = 0; i < size_uchar; i++) {
  //          cout << (uchar)i << " " << codetable[i].code << " " << codetable[i].bits << endl;
  //        }
  alphabetrank = new THuffAlphabetRank(bwt,n,codetable,0,factor);

  // Computing the position in F where the terminator is.
  terminator_pos = C[alphabetrank->charAtPos(terminator_pos)] + alphabetrank->rank(alphabetrank->charAtPos(terminator_pos), terminator_pos)-1;
  
  ulong *sampledpositions = (ulong *)calloc(this->n/W+1, sizeof(ulong));
  this->suffixes = new ulong[(this->n-1)/samplerate+1];
  this->positions = new ulong[(this->n-1)/samplerate+2];
  for (i=0;i<this->n/W+1;i++)
    sampledpositions[i]=0u;
  ulong sa_value = this->n-1;
  ulong pos = terminator_pos;

  for (i=0;i<this->n;i++) {
    if (sa_value % samplerate == 0) {
      SetField(sampledpositions,1,pos,1);
      //           this->suffixes[pos] = sa_value;
      this->positions[sa_value/samplerate] = pos;
    }
    pos = C[alphabetrank->charAtPos(pos)] + alphabetrank->rank(alphabetrank->charAtPos(pos), pos)-1;
    sa_value--;
  }

  sortPositions(positions,suffixes,(this->n-1)/samplerate+1, samplerate);
  positions[(this->n-1)/samplerate+1] = positions[0];

  this->sampled = new BitRankW32Int(sampledpositions,this->n,true,factor);
  //        delete [] bwt;
  *error=0;
}

ulong TFMindex::count(uchar *pattern, ulong m) {
  // use the FM-search replacing function Occ(c,1,i) with alphabetrank->rank(c,i)
  ulong sp, ep;
  getRange(pattern, m, sp, ep);
  if (sp<=ep) {
    return ep-sp+1;
  } else {
    return 0;
  }
}

ulong TFMindex::countReverseComp(uchar *pattern, ulong m) {
  ulong sp = 1;
  ulong ep = get_length();
  for (ulong i = 0; i < m && sp <= ep ; i++) {
    iterCountReverseComp(pattern[i], sp, ep);
  }
  if (sp <= ep)
    return ep - sp + 1;
  return 0;
}

void TFMindex::getRange(uchar *pattern, ulong m, ulong &sp, ulong &ep) {
  uchar *pat=pattern;
  if (map0 != 0 )
    pat = pattern0(pattern,m);
  uchar c = pat[m-1]; ulong i=m-1;
  sp = C[c];
  ep = C[c+1]-1;
  while (sp<=ep && i>=1) {
    c = pat[--i];
    // cout << "sp = " << sp <<", ep = " << ep << endl;
    iterCount(c, sp, ep);
  }
  if (map0 != 0) free(pat);
}

void TFMindex::getReverseCompRange(uchar *pattern, ulong m, ulong &sp, ulong &ep) {
  uchar c = complementDNA(pattern[0]);
  sp = C[c];
  ep = C[c+1]-1;
  for (ulong i = 1; i < m && sp <= ep ; i++) {
    iterCountReverseComp(pattern[i], sp, ep);
  }
}

void TFMindex::iterCount(uchar c, ulong &sp, ulong &ep) {
  sp = C[c]+ (sp == 0 ? 0 : alphabetrank->rank(c,sp-1));
  ep = C[c]+alphabetrank->rank(c,ep);
  // Before decrementing ep, make sure it is non-zero
  if (ep > 0)
    ep--;
  else {
    // ep should be negative (but we can't), meaning that we have an empty
    // range. We modify sp so that we have an empty range (buurk!)
    sp = 1;
  }
}

void TFMindex::iterCountReverseComp(uchar c, ulong &sp, ulong &ep) {
  uchar c2 = complementDNA(c);
  iterCount(c2, sp, ep);
}

ulong TFMindex::locateReverseComp(uchar *pattern, ulong m, ulong **occ,ulong nb_displayed) {
  ulong sp, ep;
  getReverseCompRange(pattern, m, sp, ep);
  if (sp <= ep)
    return returnLocate(sp, ep, nb_displayed, occ);
  else
    return 0;
}

ulong TFMindex::returnLocate(ulong sp, ulong ep, ulong nb_displayed, ulong **occ, long m) {
  ulong locate=0;
  uchar c;
  ulong i=sp, j;
  ulong dist,rank_tmp;
  *occ = (ulong *)malloc(sizeof(ulong)*(ep-sp+1));
  while (i<=ep && locate < nb_displayed) {
    j=i,dist=0;
    while (!sampled->IsBitSet(j)) {
      c = alphabetrank->charAtPos2(j,&rank_tmp);
      j = C[c]+rank_tmp; // LF-mapping The rank_tmp is already decremented by 1
      ++dist;
    }
    (*occ)[locate]=suffixes[sampled->rank(j)-1]+dist+m;
    locate ++ ;
    ++i;
  }
  return locate;
}
  
ulong TFMindex::locate(uchar *pattern, ulong m, ulong **occ,ulong nb_displayed) {
  // use the FM-search replacing function Occ(c,1,i) with alphabetrank->rank(c,i)
  ulong sp, ep, i;
  uchar c;
  getRange(pattern, m, sp, ep);
  if (nb_displayed == 0)
    nb_displayed = ep-sp+1;

  if (sp<=ep) {
    ulong matches = ep-sp+1;
    ulong locate=0;
    *occ = (ulong *) malloc(matches*sizeof(ulong));
    i=sp;
    ulong j,dist,rank_tmp;
    while (i<=ep && locate < nb_displayed) {
      j=i,dist=0;
      while (!sampled->IsBitSet(j)) {
        c = alphabetrank->charAtPos2(j,&rank_tmp);
        j = C[c]+rank_tmp; // LF-mapping The rank_tmp is already decremented by 1
        ++dist;
      }
      (*occ)[locate]=suffixes[sampled->rank(j)-1]+dist;
      locate ++ ;
      ++i;
    }
    return locate;
  } else {
    *occ = NULL;
    return 0;
  }
}

int TFMindex::display(uchar *pattern, ulong m, ulong numc, ulong *numocc, uchar **snippet_text, ulong **snippet_lengths, ulong **occ) {
  // use the FM-search replacing function Occ(c,1,i) with alphabetrank->rank(c,i)
  ulong sp, ep;
  getRange(pattern, m, sp, ep);
  uchar c; ulong i;
  if (sp<=ep) {
    (*numocc) = ep-sp+1;
    ulong locate=0;
    *occ = (ulong *) malloc((*numocc)*sizeof(ulong));
    *snippet_lengths = (ulong *) malloc((*numocc)*sizeof(ulong));
    *snippet_text = (uchar *) malloc((*numocc)*(m+2*numc)*sizeof(uchar));
    uchar *text_aux=*snippet_text;

    i=sp;
    ulong j,dist,x,rank_tmp;
    while (i<=ep) {
      j=i,dist=0;
      while (!sampled->IsBitSet(j)) {
        c = alphabetrank->charAtPos2(j,&rank_tmp);
        j = C[c]+rank_tmp; // LF-mapping The rank_tmp is already decremented by 1
        ++dist;
      }
      x = suffixes[sampled->rank(j)-1]+dist;
      (*occ)[locate]= x;
      ulong from;
      if (x>numc) from = x-numc;
      else from=0;
      ulong to=min(x+m+numc-1,n-2); // n-1 is '\0'
      ulong len =to-from+1;
      ulong skip;
      ulong j = positions[to/samplerate+1];
      if ((to/samplerate+1) == ((n-1)/samplerate+1))
        skip = n-1 - to;
      else
        skip = samplerate-to%samplerate-1;
      for (ulong dist=0;dist<skip+len;dist++) {
        c = alphabetrank->charAtPos2(j,&rank_tmp);
        j = C[c]+rank_tmp; // LF-mapping The rank_tmp is already decremented by 1
        if (dist>= skip) {
          if (c == map0)
            text_aux[len-dist-1+skip]='\0';
          else 
            text_aux[len-dist-1+skip]=c;
        }
      }
      (*snippet_lengths)[locate]=len;
      text_aux+=m+2*numc;
             
      locate ++ ;
      i++;
    }
  } else {
    *snippet_lengths = NULL;
    *snippet_text = NULL;
    *occ = NULL;
    (*numocc) = 0;
  }
  return 0;
}

ulong TFMindex::get_length () {
  return n-1;
}

TFMindex::~TFMindex() {
  delete alphabetrank;
  delete sampled;
  delete [] codetable;
  delete [] suffixes;
  delete [] positions;
}

ulong TFMindex::SpaceRequirement(){
  return sizeof(TFMindex)+alphabetrank->SpaceRequirementInBits()/8+
    sampled->SpaceRequirementInBits()/8+sizeof(codetable)*size_uchar+
    ((this->n-1)/samplerate+1)*2*sizeof(ulong);
}

ulong TFMindex::SpaceRequirement_locate(){
  return sizeof(TFMindex)+alphabetrank->SpaceRequirementInBits()/8+
    sampled->SpaceRequirementInBits()/8+sizeof(codetable)*size_uchar+
    ((this->n-1)/samplerate+1)*sizeof(ulong);
}

ulong TFMindex::SpaceRequirement_count(){
  return sizeof(TFMindex)+alphabetrank->SpaceRequirementInBits()/8+
    sizeof(codetable)*size_uchar;
}

int TFMindex::extract(ulong from, ulong to, uchar **snippet) {
  if (to> n) to=n-2; // n-1 is '\0'
  if (from > to) {
    *snippet = NULL;
    return 0;
  }
  ulong len =to-from+1;
  *snippet = (uchar *) malloc((len+1)*sizeof(uchar));
  (*snippet)[len] = 0;
  uchar c;
  ulong skip,rank_tmp;
  ulong j = positions[to/samplerate+1];
  if ((to/samplerate+1) == ((n-1)/samplerate+1)) 
    skip = n-1 - to;
  else 
    skip = samplerate-to%samplerate-1;
  for (ulong dist=0;dist<skip+len;dist++) {
    c = alphabetrank->charAtPos2(j,&rank_tmp);
    j = C[c]+rank_tmp; // LF-mapping The rank_tmp is already decremented by 1
    if (dist>= skip) {
      if (c == map0)
        (*snippet)[len-dist-1+skip]='\0';
      else 
        (*snippet)[len-dist-1+skip]=c;
    }
  }
  return len;
}

int TFMindex::save(char *filename) {
  char fnamext[1024];        
  FILE *f;
  sprintf (fnamext,"%s.ssa",filename);
  f = fopen (fnamext,"w");
  if (f == NULL) return 20;
  if (fwrite (&n,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&samplerate,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&map0,sizeof(uchar),1,f) != 1) return 21;
  if (fwrite (C,sizeof(ulong),size_uchar+1,f) != size_uchar+1) return 21;
  if (save_codetable(f,codetable) !=0) return 21;
  if (alphabetrank->save(f) !=0) return 21;
  if (sampled->save(f) !=0) return 21;
  if (fwrite (suffixes,sizeof(ulong),(this->n-1)/samplerate+1,f) != (this->n-1)/samplerate+1) return 21;
  if (fwrite (positions,sizeof(ulong),(this->n-1)/samplerate+2,f) != (this->n-1)/samplerate+2) return 21;
  fclose(f);
  return 0;
}

int TFMindex::load(char *filename) {
  char fnamext[1024];        
  FILE *f;
  int error;
  sprintf (fnamext,"%s.ssa",filename);
  f = fopen (fnamext,"r");
  if (f == NULL) return 23;
  if (fread (&n,sizeof(ulong),1,f) != 1) return 25;
  if (fread (&samplerate,sizeof(ulong),1,f) != 1) return 25;
  if (fread (&map0,sizeof(uchar),1,f) != 1) return 25;
  if (fread (C,sizeof(ulong),size_uchar+1,f) != size_uchar+1) return 25;
  if (load_codetable(f,&codetable) !=0) return 25;
  alphabetrank = new THuffAlphabetRank(f,codetable,&error);
  if (error !=0) return error;
  sampled = new BitRankW32Int(f,&error);
  if (error !=0) return error;
  this->suffixes = new ulong[(this->n-1)/samplerate+1];
  if (!this->suffixes) return 1;
  if (fread (suffixes,sizeof(ulong),(this->n-1)/samplerate+1,f) != (this->n-1)/samplerate+1) return 25;
  this->positions = new ulong[(this->n-1)/samplerate+2];
  if (!this->positions) return 1;
  if (fread (positions,sizeof(ulong),(this->n-1)/samplerate+2,f) != (this->n-1)/samplerate+2) return 25;
  fclose(f);
  return 0;
}

TFMindex::TFMindex(char *filename, int *error) {
  *error = load(filename);
}

////////////////////////
////Building the Index//
////////////////////////
int build_index(uchar *text, ulong length, char *build_options, void **index){
  ulong samplerate=64;
  ulong factor=4; 
  int error;
  char delimiters[] = " =;";
  int j,num_parameters;
  char ** parameters;
  int free_text=false; /* don't free text by default */
  uchar terminator = 0;

  if (build_options != NULL) {
    parse_parameters(build_options,&num_parameters, &parameters, delimiters);
    for (j=0; j<num_parameters;j++) {
      if ((strcmp(parameters[j], "samplerate") == 0 ) && (j < num_parameters-1) ) {
        samplerate=atoi(parameters[j+1]);
        j++;
      } else if ((strcmp(parameters[j], "factor") == 0 ) && (j < num_parameters-1) ) {
        factor=atoi(parameters[j+1]);
        j++;
      } else if (strcmp(parameters[j], "free_text") == 0 )
        free_text=true;
      else if (strcmp(parameters[j], "terminator") == 0 && (j < num_parameters - 1))
        terminator = parameters[j+1][0];
    }
    free_parameters(num_parameters, &parameters);
  }
  TFMindex *FMindex = new TFMindex(text,length,samplerate,free_text,factor,&error, terminator);
  (*index) = FMindex;
  if (error != 0) return error;
  return 0;
}
int save_index(void *index, char *filename) {
  TFMindex *_index=(TFMindex *) index;
  return _index->save(filename);
}
int load_index(char *filename, void **index){
  int error;
  TFMindex *FMindex = new TFMindex(filename, &error);
  (*index) = FMindex;
  return error;
}

int free_index(void *index){
  TFMindex *_index=(TFMindex *) index;
  delete _index;
  return 0;
}

int index_size(void *index, ulong *size) {
  TFMindex *_index=(TFMindex *) index;
  (*size) = _index->SpaceRequirement();
  return 0;
}
int index_size_count(void *index, ulong *size) {
  TFMindex *_index=(TFMindex *) index;
  (*size) = _index->SpaceRequirement_count();
  return 0;
}
int index_size_locate(void *index, ulong *size) {
  TFMindex *_index=(TFMindex *) index;
  (*size) = _index->SpaceRequirement_locate();
  return 0;
}
////////////////////////
////Querying the Index//
////////////////////////
int count(void *index, uchar *pattern, ulong length, ulong *numocc){
  TFMindex *_index=(TFMindex *) index;
  (*numocc)= _index->count(pattern,length);
  return 0;
}
int locate(void *index, uchar *pattern, ulong length, ulong **occ, ulong *numocc, ulong nb_displayed){
  TFMindex *_index=(TFMindex *) index;
  (*numocc)= _index->locate(pattern,length,occ, nb_displayed);
  return 0;
}
/////////////////////////
//Accessing the indexed//
/////////////////////////
int extract(void *index, ulong from, ulong to, uchar **snippet, ulong *snippet_length){
  TFMindex *_index=(TFMindex *) index;
  (*snippet_length)= _index->extract(from,to, snippet);
  return 0 ;
}
int display(void *index, uchar *pattern, ulong length, ulong numc, ulong *numocc, uchar **snippet_text, ulong **snippet_lengths) {
  TFMindex *_index=(TFMindex *) index;
  ulong *occ;
  int aux = _index->display(pattern, length, numc, numocc, snippet_text, snippet_lengths, &occ);
  free(occ);
  return aux;
}
int get_length(void *index, ulong *length){
  TFMindex *_index=(TFMindex *) index;
  (*length)=(_index)->get_length();
  return 0;
}
////////////////////
////Error handling//
////////////////////
const char *error_index(int e){
  switch(e) {
    case 0:  return "No error"; 
    case 1:  return "Out of memory"; 
    case 2:  return "The text must end with a \\0"; 
    case 5:  return "You can't free the text if you don't copy it"; 
    case 20: return "Cannot create files"; 
    case 21: return "Error writing the index"; 
    case 22: return "Error writing the index"; 
    case 23: return "Cannot open index";
    case 24: return "Cannot open text";
    case 25: return "Error reading the index";
    case 26: return "Error reading the index";
    case 27: return "Error reading the text"; 
    case 28: return "Error reading the text"; 
    case 99: return "Not implemented"; 
    default: return "Unknown error";
  }
}

