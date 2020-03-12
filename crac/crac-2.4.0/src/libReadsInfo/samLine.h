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

#ifndef SAMLINE_H
#define SAMLINE_H
#include <string>
#include <map>
#include <ostream>
#include <vector>

#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "Cigar.h"

using namespace std;

class OptionalField {
private:
  char tag[2];
  char type;
public:

  OptionalField(const char tag_name[2], char type);

  OptionalField(const OptionalField& o);

  virtual OptionalField* clone() = 0; 

  virtual ~OptionalField() {}

  const char *getTag() const;

  char getType() const;

  void saveInBamRecord(kstring_t *str);

  virtual void saveValueInBamRecord(kstring_t *str) = 0;

  virtual std::ostream& printValue(std::ostream& os) const = 0;
};

std::ostream& operator<<(std::ostream& os, const OptionalField& of);

class CharOptionalField : public OptionalField {
private:
  char value;
public:
  CharOptionalField(const char tag[2], char value);
  CharOptionalField(const CharOptionalField& o);
  ~CharOptionalField();
  virtual OptionalField* clone(); 
  char getValue() const;
  virtual void saveValueInBamRecord(kstring_t *str);
  virtual std::ostream& printValue(std::ostream& os) const;
};

class IntOptionalField : public OptionalField {
private:
  int value;
public:
  IntOptionalField(const char tag[2], int value);
  IntOptionalField(const char tag[2], uint value);
  IntOptionalField(const IntOptionalField& o);
  ~IntOptionalField();
  virtual OptionalField* clone(); 
  int getValue() const;
  virtual void saveValueInBamRecord(kstring_t *str);
  virtual std::ostream& printValue(std::ostream& os) const;
};

class StringOptionalField : public OptionalField {
private:
  string value;
public:
  StringOptionalField(const char tag[2], string value);
  StringOptionalField(const char tag[2], char* value);
  StringOptionalField(const StringOptionalField& o);
  ~StringOptionalField();
  virtual OptionalField* clone(); 
  const string& getValue() const;
  virtual void saveValueInBamRecord(kstring_t *str);
  virtual std::ostream& printValue(std::ostream& os) const;
};

// TODO Implement an IntArray optional field sub-class
// This could be a good way to store p_support and p_loc fields
//class IntArrayOptionalField : public OptionalField {
//  private:
//    int *value;
//public:
//
//  IntArrayOptionalField(const char *tag, int* value) 
//  : OptionalField(tag,'B')
//  {
//    value = value; 
//  };
//
//  virtual ~IntArrayOptionalField();
//
//};

typedef map<string, OptionalField* > OF_map; 

class SamLine {
private:
  //GenomeInfo *genome;
  string rname;
  string qname;
  string rnext;
  string seq;
  string qual;
  int rid;
  int ridnext;
  uint flag;
  uint pos;
  uint mapQ;
  uint pnext;
  int tlen;
  Cigar cigar;
  //map<string, string> optionalFields;
  OF_map optionalFields;

  bool isFlagBitSet(uint bit) const;
  void setFlagBit(uint bit);
  void unsetFlagBit(uint bit);

public:
  /*
   * Default constructor that initialize all fields to default value
   */
  SamLine();

  /*
   * Recopy constructor
   */
  SamLine(const SamLine& s);

  // Desctructor
  ~SamLine();

  // Query template NAME
  void setQname(const string& qname);
  const string& getQname() const;

  // bitwize FLAG
  void setFlag(uint flag);
  uint getFlag() const;

  // Flag operations
  // Bit Description
  //   0x1 template having multiple segments in sequencing
  //   0x2 each segment properly aligned according to the aligner
  //   0x4 segment unmapped
  //   0x8 next segment in the template unmapped
  //   0x10 SEQ being reverse complemented
  //   0x20 SEQ of the next segment in the template being reversed
  //   0x40 the rst segment in the template
  //   0x80 the last segment in the template
  //   0x100 secondary alignment
  //   0x200 not passing quality controls
  //   0x400 PCR or optical duplicate
  //   0x800 supplementary alignment
  // bit 1
  bool isTemplateHavingMultipleSegments() const;
  void setTemplateHavingMultipleSegments();
  void unsetTemplateHavingMultipleSegments();

  // bit 2
  bool isEachSegmentsMapped() const;
  void setEachSegmentsMapped();
  void unsetEachSegmentsMapped();

  // bit 4
  bool isSegmentUnmapped() const;
  void setSegmentUnmapped();
  void unsetSegmentUnmapped();

  // bit 8
  bool isNextSegmentUnmapped() const;
  void setNextSegmentUnmapped();
  void unsetNextSegmentUnmapped();

  // bit 16
  bool isSeqReverseComplemented() const;
  void setSeqReverseComplemented();
  void unsetSeqReverseComplemented();

  // bit 32
  bool isNextSeqReverseComplemented() const;
  void setNextSeqReverseComplemented();
  void unsetNextSeqReverseComplemented();

  // bit 64
  bool isFirstSegmentInTheTemplate() const;
  void setFirstSegmentInTheTemplate();
  void unsetFirstSegmentInTheTemplate();

  // bit 128
  bool isLastSegmentInTheTemplate() const;
  void setLastSegmentInTheTemplate();
  void unsetLastSegmentInTheTemplate();

  // bit 256
  bool isSecondaryAlignement() const;
  void setSecondaryAlignement();
  void unsetSecondaryAlignement();

  // bit 512
  bool isFailingQualityControl() const;
  void setFailingQualityControl();
  void unsetFailingQualityControl();

  // bit 1024
  bool isPCRDuplicated() const;
  void setPCRDuplicated();
  void unsetPCRDuplicated();

  // bit 2048
  bool isChimericAlignement() const;
  void setChimericAlignement();
  void unsetChimericAlignement();

  // Reference sequence NAME
  // This id correspond to the ID given for this reference sequence in
  // the SAM headers
  void setRname(const string& rname, const int rid);
  void setUnknownRname();
  const string& getRname() const;
  int getRid() const;

  // 1-based leftmost mapping POSition
  void setPos(uint pos);
  void setUnknownPos();
  uint getPos() const;

  // MAPping Quality
  void setMapQ(uint mapQ); 
  uint getMapQ() const;

  // CIGAR string
  void setCigar(const Cigar& cigar);
  const Cigar& getCigar() const;

  // Ref. name of the mate/next segment
  void setRnext(const string& rnext, const int ridnext);
  const string& getRnext() const;
  int getRidnext() const;

  // Position of the mate/next segment
  void setPnext(uint pnext);
  uint getPnext() const;

  // observed Template LENgth
  void setTlen(int tlen);
  int getTlen() const;

  // segment SEQuence
  void setSeq(const string& seq);
  const string& getSeq() const;
  void reverseComplementeSeq();

  // ASCII of Phred-scaled base QUALity+33
  void setQual(const string& qual);
  const string& getQual() const;
  void reverseQual();

  // Optional fields
  void addOptionalField(const char tag[], char val);
  void addOptionalField(const char tag[], uint val);
  void addOptionalField(const char tag[], int val);
  void addOptionalField(const char *tag, const string& val);
  void addOptionalField(const char *tag, const char* val);

  // How could we re-implement this methode with the new
  // abstract class OptionalField?
  //const string& getOptionalField(const char *tag);
  const OF_map& getOptionalFields() const;
  bool isOptionalFieldDefined(const char *tag);

  // Remove all fields that are user reserved.
  // That means the tag used matches : X?, Y? or Z?
  void removeAllUserOptionalFields();

  // write the SAM line in the file
  ostream &writeLine(ostream &os);

  int writeBamRecord(samFile *out, const bam_hdr_t *h);

private:
};

#endif
