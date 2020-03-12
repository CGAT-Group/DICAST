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
#include <iostream>
#include "readsReader.h"

using namespace std;


void readIterator::setPrintWarnings(bool isVisible)
{
  this->warningsVisible = isVisible;
}

bool readIterator::printWarnings() {
  return warningsVisible;
}

singleReadIterator::singleReadIterator(const char* filename, uint k, 
                                       uint length, bool printWarnings, 
                                       bool autoDiscard, 
                                       bool autoFirstIteration):
k(k),length(length),autoDiscard(autoDiscard) {
  readsFile = gzopen(filename, "r");
  seq = kseq_init(readsFile);
  finished = false;
  filenameConstructor = true;
  tagNumber = 0;
  setPrintWarnings(printWarnings);
  if (autoFirstIteration) ++(*this);
}

singleReadIterator::singleReadIterator(kseq_t *seq, uint k, uint length):k(k),
  length(length), 
  seq(seq),
  filenameConstructor(false) {
  autoDiscard = true;
  tagNumber = 0;
  setPrintWarnings(false);
}

singleReadIterator::singleReadIterator(const singleReadIterator& mit):k(mit.k),
  length(mit.length),
  seq(mit.seq),
  finished(mit.finished),
  filenameConstructor(false),
  autoDiscard(mit.autoDiscard),
  tagNumber(mit.tagNumber)
{}

singleReadIterator::~singleReadIterator() {
  if (seq && filenameConstructor) {
    kseq_destroy(seq);
    gzclose(readsFile);
  }
}

readIterator& singleReadIterator::operator++() {
  
  if (finished || kseq_read(seq) <= 0) {
    finished = true;
  } else {
    tagNumber++;
    if (isDiscarded() && autoDiscard) {
      if(printWarnings()) {
        cerr << "Warning: read number " << tagNumber 
            << " is ignored";
        if (k > getLength()) {
          cerr << " because its length is lower than the k-mer length"
              << endl;
        } else if (getLength() < length) {
          cerr << " because its length is shorter than expected"
              << endl;
        }
      }
      // Read should be ignored (its length isn't good), go to the next one
      ++(*this);
    }
  }
  return *this;
}

singleReadIterator singleReadIterator::operator++(int)  {
  singleReadIterator orig = *this;
  ++(*this);
  return orig;
}

bool singleReadIterator::operator==(const singleReadIterator& rhs)  {
  return ((seq == NULL || finished) && (rhs.seq == NULL || rhs.finished))
    || (seq == rhs.seq);
}
bool singleReadIterator::operator!=(const singleReadIterator& rhs)  {
  return ! (*this == rhs);
}

kseq_t &singleReadIterator::operator*()  {
  return *seq;
}

size_t singleReadIterator::getLength() {
  return seq->seq.l;
}

char *singleReadIterator::getName()  {
  return seq->name.s;
}

char *singleReadIterator::getQuality()  {
  return seq->qual.s;
}

uint singleReadIterator::getReadNumber() {
  // minus 1 because we want the numbering to start from 0
  return tagNumber-1;
}

char *singleReadIterator::getSequence()  {
  return seq->seq.s;
}

bool singleReadIterator::isFinished()  {
  return finished;
}

bool singleReadIterator::isDiscarded() {
  return (getLength() < length  || getLength() < k);
}

void singleReadIterator::setAutoDiscard(bool value){
  autoDiscard = value;
} 

//////////////////////////////////////////////////////////// 

pairedEndReadIterator::pairedEndReadIterator(const char* filename1, 
                                             const char* filename2, 
                                             uint k, uint length, 
                                             bool printWarnings,
                                             bool removePEIdentifiers):
k(k),length(length)
{
  this->readIterator1 = new singleReadIterator(filename1,k,length,printWarnings,false,false);
  this->readIterator2 = new singleReadIterator(filename2,k,length,printWarnings,false,false);
  this->is_it1CurrentIterator = false;
  this->setPrintWarnings(printWarnings);
  this->remove_paired_end_identifiers = removePEIdentifiers;
  ++(*this);
}

pairedEndReadIterator::pairedEndReadIterator(const pairedEndReadIterator& mit){
  this->readIterator1 = mit.readIterator1;
  this->readIterator2 = mit.readIterator2;
  this->is_it1CurrentIterator = mit.is_it1CurrentIterator;
}

pairedEndReadIterator::~pairedEndReadIterator()
{
  delete readIterator1;
  delete readIterator2;
}

readIterator& pairedEndReadIterator::operator++() 
{
  bool isReadPairValid = false;

  if(!isFinished()) {
    if(!is_it1CurrentIterator) {
      while(!isReadPairValid && !isFinished()) {
        ++(*readIterator1);
        ++(*readIterator2);
        if(isFinished()) break;
        if(!readIterator1->isDiscarded() && !readIterator2->isDiscarded()){
          isReadPairValid = true;
        } else if(printWarnings()) {
          cerr << "Warning: both reads of the pair number " << readIterator1->getReadNumber()
              << " are ignored";
          if (k > readIterator1->getLength()) {
            cerr << " because the first read length is lower than the k-mer length"
                << endl;
          } else if (k > readIterator2->getLength()){
            cerr << " because the second read length is lower than the k-mer length"
                << endl;
          } else if (readIterator1->getLength() < length) {
            cerr << " because the first read length is shorter than expected"
                << endl;
          } else if (readIterator2->getLength() < length) {
            cerr << " because the second read length is shorter than expected"
                << endl;
          }
        }
      }
      is_it1CurrentIterator = true;
    } else {
      is_it1CurrentIterator = false;
    }
  }
  return *this;
}

pairedEndReadIterator pairedEndReadIterator::operator++(int) 
{
  pairedEndReadIterator orig = *this;
  ++(*this);
  return orig;
}

bool pairedEndReadIterator::operator==(const pairedEndReadIterator& rhs)  {
  if(is_it1CurrentIterator){
    return readIterator1 == rhs.readIterator1;
  } else {
    return readIterator2 == rhs.readIterator2;
  }
}

bool pairedEndReadIterator::operator!=(const pairedEndReadIterator& rhs) { 
  return !(*this == rhs);
}

kseq_t &pairedEndReadIterator::operator*() {
  if(is_it1CurrentIterator){
    return *(*readIterator1);
  } else {
    return *(*readIterator2);
  }
}  

size_t pairedEndReadIterator::getLength() {
  return getCurrentReadIterator()->getLength();
}

char *pairedEndReadIterator::getName() {
  char *name = getCurrentReadIterator()->getName();
  if(name) {
    // We remove /1 and /2 paired-end identifiers, if there is any
    int l = strlen(name);
    if(remove_paired_end_identifiers && l > 2 && name[l-2] == '/' && 
       (name[l-1] == '1' || name[l-1] == '2')) {
      name[l-2] = '\0';
    }
  }
  return name;
}

readIterator& pairedEndReadIterator::getPair() {
  if(is_it1CurrentIterator) {
    is_it1CurrentIterator = false;
  } else {
    is_it1CurrentIterator = true;
  }
  return *this;
}

char *pairedEndReadIterator::getQuality() {
  return getCurrentReadIterator()->getQuality();
}

uint pairedEndReadIterator::getReadNumber() {
  if(is_it1CurrentIterator) {
    return readIterator1->getReadNumber() + readIterator2->getReadNumber();
  } else {
    return readIterator1->getReadNumber() + readIterator2->getReadNumber() + 1;
  }
}

char *pairedEndReadIterator::getSequence() {
  return getCurrentReadIterator()->getSequence();
}

bool pairedEndReadIterator::isFinished() {
  return readIterator1->isFinished() || readIterator2->isFinished();
}

bool pairedEndReadIterator::isTheFirstMemberOfPair() {
  return is_it1CurrentIterator;
}

singleReadIterator *pairedEndReadIterator::getCurrentReadIterator() {
  if(is_it1CurrentIterator) {
    return readIterator1;
  } else {
    return readIterator2;
  }
}


readsReader::readsReader(char *filename, uint k, uint length)
{
  this->filename1 = filename;
  this->k = k;
  this->length = length;
  this->is_pairedEnd = false;
}

readsReader::readsReader(char *filename1, char *filename2, uint k, uint length)
{
  this->filename1 = filename1;
  this->filename2 = filename2;
  this->k = k;
  this->length = length;
  this->is_pairedEnd = true;
}

readIterator *readsReader::begin(bool printWarnings) {
  if(is_pairedEnd){
    return new pairedEndReadIterator(filename1,filename2,k,length,printWarnings);
  } else {
    return new singleReadIterator(filename1,k,length,printWarnings);
  }
}

bool readsReader::isPairedEnd() {
  return is_pairedEnd;
}
