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

#ifndef READSREADER_H
#define READSREADER_H

#include <iterator>
#include <zlib.h>
#include "kseq/kseq.h"

KSEQ_INIT(gzFile, gzread)

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
  bool remove_paired_end_identifiers;
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
                        uint length=0, bool printWarnings=false, bool removePEIdentifiers=true);

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

#endif
