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

#include "utils.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <gzstream.h>
#include <cstring>

using namespace std;


uintSA DNAtoInt(char *dna, uintSA dna_length){
  uintSA dna_int = 0;
  for (uintSA i=0 ; i< dna_length ; i++){
    dna_int <<= 2;      
    dna_int |= convNuc(dna[i]);
  }
  return dna_int;
}

uint convNuc(char nuc){
   switch (nuc){
   case 'a' : case 'A' : return 0 ;
   case 'c' : case 'C' : return 1 ;
   case 'g' : case 'G' : return 2 ;
   case 't' : case 'T' : return 3 ;
   default : cerr << "invalid nucleotide "<< nuc <<endl; exit(4);
   }  
   return 0;
}

char majNuc(char nuc) {
   switch (nuc){
   case 'a' : return 'A';
   case 'c' : return 'C';
   case 'g' : return 'G';
   case 't' : return 'T';
   default : return nuc;
   }  
   return 0;
}

char *majNucNCpy(char *dest, const char *source, size_t n) {
    char *start = dest;

      while (n && (*dest++ = majNuc(*source++))) n--;
        if (n) while (--n) *dest++ = '\0';
          return start;
}

char intToNuc(uint c) {
  switch(c) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  default: cerr << "Invalid number (0-3): " << c <<endl; exit(5);
  }
}

void intToDNA(uint64 code, uint dna_length, char *dna) {
  uint64 mask = 3;
  for (uint i=0; i < dna_length; i++) {
    dna[dna_length-i-1] = intToNuc(code & mask);
    code >>=2;
  }
}

uint64 factorsToInt(uintSA pos, uchar *dna, uint length) {
  uintSA posFirstFactor = 4 - (pos % 4);
  uintSA numFirstFactor = pos / 4;
  uintSA posLastFactor = (length - posFirstFactor) % 4;
  uint64 factor = dna[numFirstFactor] & ((1 << (posFirstFactor*2)) - 1);

  for(uintSA i = 1 ; i <= (length-posFirstFactor)/4 ; i++){
    factor <<= 8;
    factor |= dna[numFirstFactor+i];
  }
  factor <<= posLastFactor*2;

  factor |= dna[numFirstFactor+(length-posFirstFactor)/4+1] >> (8-2*posLastFactor) & ((1 << (posLastFactor*2)) - 1);
  
  return factor;
}

int comparUint(const void *a1, const void* a2){
  if ((* (uint *) a1) < (* (uint *) a2)){
    return -1;
  }else{
    return ((* (uint *) a1) == (* (uint *) a2)) ? 0 : 1; 
  }
}

bool pairCompare(const std::pair<uint, uint>& firstElem, const std::pair<uint, uint>& secondElem) {
    return firstElem.first < secondElem.first;
}    

uint posBitInMask(uint mask) {
  uint i = 1;
  uint pos = 0;
  if (mask == 0)
    return sizeof(uint)*8;
  while ((mask & i) == 0) {
    i <<= 1;
    pos++;
  }
  return pos;
}

float getChrono() {
	struct timeval time;
	time_t sec;
	suseconds_t usec;

	if (gettimeofday(&time,NULL) != -1) {
		sec = time.tv_sec;
		usec = time.tv_usec;
		return (float)sec*1000000+usec;
	}
	return 0;
}

string intToString(uint i) {
  ostringstream oss;
  oss << i;
  return oss.str();
}

ostream *create_stream(bool gzipped, char *name) {
  if (gzipped) {
    uint len = strlen(name);
    char *gzName = new char[len+4];
    ostream *result;
    sprintf(gzName, "%s.gz", name);
    result = new ogzstream((const char *)gzName, ofstream::out);
    delete [] gzName;
    return result;
  }
  return new ofstream((const char *)name, ofstream::out);
}

std::string concatStringAndInt(const char*str, int i) {
  char buffer[10];
  sprintf(buffer, "%d", i);
  std::string output = str;
  output += buffer;
  return output;
}

std::string concatStringAndFloat(const char*str, float i) {
  char buffer[100];
  sprintf(buffer, "%f", i);
  std::string output = str;
  output += buffer;
  return output;
}
