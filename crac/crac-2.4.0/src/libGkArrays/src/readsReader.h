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

#ifndef READSREADER_H
#define READSREADER_H

#include <iterator>
#include <zlib.h>
#include "kseq/kseq.h"

KSEQ_INIT(gzFile, gzread)

namespace gkarrays {
  /**
   * readIterator allows to retrieve information on each read by traversing
   * them in order, one after one.
   * If the k-mer length or the read length are given, only reads that are
   * long enough to fulfill those values are considered.
   */
  class readIterator: public std::iterator<std::forward_iterator_tag, kseq_t> {
  bool warningsVisible;
  public:
  /**
   * Virtual destructor, necessary in c++
   */
  virtual ~readIterator(){}

  /**
   * Go to next sequence
   */
  virtual readIterator& operator++() = 0;

  /**
   * @return the kseq_t related to the current read
   */
  virtual kseq_t& operator*() = 0;

  /**
   * @return the length of the read
   */
  virtual size_t getLength() = 0;

  /**
   * @return the name of the read
   */
  virtual char *getName() = 0;

  /**
   * @return the quality of the read
   */
  virtual char *getQuality() = 0;

  /**
   * @return the number of the read in the file
   */
  virtual uint getReadNumber() = 0;

  /**
   * @return the sequence of the read
   */
  virtual char *getSequence() = 0;

  /**
   * @return true iff we have read the whole file
   */
  virtual bool isFinished() = 0;

  /**
   * @return true if warnings are visible
   *         false either.
   */
  bool printWarnings();

  /**
   * @param isVisible Printwarnins visible value.
   */
  void setPrintWarnings(bool isVisible);
  };

  /**
  * singleReadIterator allows to retrieve information on each read by traversing
  * them in order, one after one.
  * If the k-mer length or the read length are given, only reads that are
  * long enough to fulfill those values are considered.
  */
  class singleReadIterator: public readIterator {
  uint k;
  uint length;
  kseq_t *seq;
  bool finished;
  bool filenameConstructor;
  gzFile readsFile;
  bool autoDiscard; // if autoDiscard is true, the iterator will skip bad reads
  uint tagNumber;
  public:
  /**
   * Constructor. Reads are retrieved from the given file.
   * @param filename: file we must iterate on
   * @param k: k-mer length used (0 if unknown or not applicable)
   * @param length: Read length (0 for variable)
   * @param printWarnings: if true print warnings when skipping a read
   * @param autoDiscard: if true skip automatically bad reads
   * @param autoFirstIteration: automatically perform the first iteration 
   *        (meaning that it initializes the data structure)
   */
  singleReadIterator(const char* filename, uint k=0, uint length=0, 
                     bool printWarnings=false, bool autoDiscard=true, 
                     bool autoFirstIteration=true);

  /**
   * Constructor. singleReadIterator starts at the given sequence.
   * @param seq: the sequence to start at.
   * @param k: k-mer length used (0 if unknown or not applicable)
   * @param length: Read length (0 for variable)
   */
  singleReadIterator(kseq_t* seq, uint k=0, uint length=0);

  /**
   * Copy constructor
   */
  singleReadIterator(const singleReadIterator& mit);

  /**
   * Destructor
   */
  ~singleReadIterator();

  /**
   * Go to next sequence prefix
   */
  virtual readIterator& operator++();
  /**
   * Go to next sequence postfix
   */
  virtual singleReadIterator operator++(int);

  virtual bool operator==(const singleReadIterator& rhs);
  virtual bool operator!=(const singleReadIterator& rhs);
  virtual kseq_t& operator*();

  /**
   * @return the name of the read
   */
  virtual char *getName();

  /**
   * @return the quality of the read
   */
  virtual char *getQuality();

  /**
   * @return the sequence of the read
   */
  virtual char *getSequence();

  /**
   * @return the number of the read in the file
   */
  virtual uint getReadNumber();

  /**
   * @return the length of the read
   */
  virtual size_t getLength();

  /**
   * @return true iff we have read the whole file
   */
  virtual bool isFinished();

  /**
   * @return true if real tag length is inferior to k or to wanted length
   */
  bool isDiscarded();

  /**
   * @param value Value to set AutoDiscard mode
   */
  void setAutoDiscard(bool value);

  };

  /**
  * pairedEndReadIterator allows to retrieve information on each read by traversing
  * them in order, one after one in the case of paired-end reads.
  * If the k-mer length or the read length are given, only reads that are
  * long enough to fulfill those values are considered.
  */
  class pairedEndReadIterator: public readIterator {
  singleReadIterator *readIterator1;
  singleReadIterator *readIterator2;
  bool is_it1CurrentIterator;
  uint k;
  uint length;

  public:
  /**
   * Constructor. Reads are retrieved from the given file.
   * @param filename1: file we must iterate on
   * @param filename2: second file (for paired-end) we must iterate on
   * @param k: k-mer length used (0 if unknown or not applicable)
   * @param length: Read length (0 for variable)
   * @param printWarnings: if true print warnings when skipping a read
   */
  pairedEndReadIterator(const char* filename1, const char* filename2, uint k=0, 
                        uint length=0, bool printWarnings=false);

  /**
   * Copy constructor
   */
  pairedEndReadIterator(const pairedEndReadIterator& mit);
  
  /**
   * Destructor
   */
  virtual ~pairedEndReadIterator();

  /**
   * Go to next sequence prefix
   */
  virtual readIterator& operator++();
  /**
   * Go to next sequence postfix
   */
  virtual pairedEndReadIterator operator++(int);
  virtual bool operator==(const pairedEndReadIterator& rhs);
  virtual bool operator!=(const pairedEndReadIterator& rhs);
  virtual kseq_t& operator*();

  /**
   * @return the length of the read
   */
  virtual size_t getLength();

  /**
   * @return the name of the read
   */
  virtual char *getName();

  /**
   * @return the paired iterator with the current one
   */
  readIterator& getPair();

  /**
   * @return the quality of the read
   */
  virtual char *getQuality();

  /**
   * @return the number of the read in the file
   */
  virtual uint getReadNumber();

  /**
   * @return the sequence of the read
   */
  virtual char *getSequence();

  /**
   * @return true iff we have read the whole file
   */
  virtual bool isFinished();
  
  /**
   * @return true if the current read is the first
   *         member of the pair
   */
  bool isTheFirstMemberOfPair();

  private:

  /**
   * @return the current singleReadIterator according to the isIt1CurrentIterator boolean
   */
  singleReadIterator *getCurrentReadIterator();
  };

  /**
   * readsReader is a class that allows you to store informations about
   * an RNA_Seq experiment. It works both for single reads and paired-end reads.
   * You can get an iterator using method begin() in order to go trough all reads.
   */
  class readsReader {

  char *filename1;
  char *filename2; // in case of paired end;
  uint k;
  uint length;
  bool is_pairedEnd;

  public :
  /**
   * Single reads constructor
   * @param filename The name of the file containg reads
   * @param k length of k-mers we have to use (0 for variable-length reads)
   * @param length length of the reads. If a shorter read is found, 
   *        it raises an error. If a longer read is found, only the prefix
   *        of tag_length characters is kept.
   *        If tag_length == 0 (default), just gess what the read length
   *        is.
   */
  readsReader(char *filename, uint k=0, uint length=0);

  /**
   * Single reads constructor
   * @param filename1 The name of the first file containing reads
   * @param filename2 The name of the second file containing reads
   * @param k length of k-mers we have to use (0 for variable-length reads)
   * @param length length of the reads. If a shorter read is found, 
   *        it raises an error. If a longer read is found, only the prefix
   *        of tag_length characters is kept.
   *        If tag_length == 0 (default), just gess what the read length
   *        is.
   */
  readsReader(char *filename1, char *filename2, uint k=0, uint length=0); 

  /**
   * @param printWarnings Value of printWarnings option
   * @return an read iterator that goes through the read file(s)
   */
  readIterator *begin(bool printWarnings=false);

  /**
   * @return true is the reads are paired-end (i.e there is two reads files)
   */
  bool isPairedEnd();

  };
}

#endif
