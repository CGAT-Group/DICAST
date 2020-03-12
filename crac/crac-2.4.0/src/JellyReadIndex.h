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


#include <config.h>

#ifdef HAVE_LIBJELLYFISH
#ifndef JELLYREADINDEX_H
#define JELLYREADINDEX_H

#include "ReadIndex.h"

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>

using jellyfish::mer_dna;
typedef jellyfish::stream_manager<char**> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > read_parser;

class mer_counter : public jellyfish::thread_exec {
  int             nb_threads_;
  mer_hash&       ary_;
  bool            canonical_;
  stream_manager  streams_;
  sequence_parser parser_;

public:
 mer_counter(int nb_threads, mer_hash& ary, bool canonical, char** file_begin, char** file_end) :
    nb_threads_(nb_threads),
    ary_(ary),
    canonical_(canonical),
    streams_(file_begin, file_end, 1), // 1: parse one file at a time
    parser_(mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    for(mer_iterator mers(parser_, canonical_) ; mers; ++mers)
      ary_.add(*mers, 1);
    ary_.done();
  }
};


class JellyReadIndex : public ReadIndex {

  mer_hash *ary;

  public:

    /*
     * Constructut a ReadIndex from a GkArrays instance
     * @param gk the gkArrays object to use
     */
  JellyReadIndex(uint k, size_t initial_hash_size, bool canonical, uint nb_threads, char *reads, char *reads2 = NULL) {
      char *files[2] = { reads, reads2 };
      uint nb_files  = reads2 ? 2 : 1;

      mer_dna::k(k);

      // 7 is the length in bits of the counter field
      // 126 is the maximum reprobe value
      ary = new mer_hash(initial_hash_size, mer_dna::k() * 2, 7, nb_threads, 126);
      // Count the k-mers in the reads and store then in the hash array 'ary'
      mer_counter counter(nb_threads, *ary, canonical, files, files + nb_files);
      counter.exec_join(nb_threads);    
    }

    /*
     * Destructor (does nothing), the gkArrays object is not destructed!
     */
    ~JellyReadIndex() {
      delete ary;
    }

    /**
     * @return the factor length used to index the read collection
     */
    virtual uint getFactorLength() const;

    /**
     * @return the number of tags (or reads) indexed in the Gk Arrays
     */
    virtual uint getNbTags() const;

    /**
     * @param r the read object we want to get the support
     * @return an array whose length is it->getLength()-getFactorLength()
     *         and where the value at position i is the number of occurrences
     *         of the k-factor starting at position i in the reads
     *         among all the Pk-factors.
     */
    virtual uint *getSupport(Read *r) const;


};

#endif //JELLYREADINDEX_H
#endif //HAVE_LIBJELLYFISH

