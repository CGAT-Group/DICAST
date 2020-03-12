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
#include <iostream>
#include "readsReader.h"

using namespace std;

namespace gkarrays {

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
                                               bool printWarnings):
  k(k),length(length)
  {
    this->readIterator1 = new singleReadIterator(filename1,k,length,printWarnings,false,false);
    this->readIterator2 = new singleReadIterator(filename2,k,length,printWarnings,false,false);
    this->is_it1CurrentIterator = false;
    this->setPrintWarnings(printWarnings);
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
    return getCurrentReadIterator()->getName();
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

}
