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

#include "samLine.h"
#include "bamUtils.h"

#include "../libSSA/utils.h"
#include <sstream>
#include <cstring>

OptionalField::OptionalField(const char tag_name[2], char t) {
  strncpy(tag,tag_name,2);
  type = t;
}

OptionalField::OptionalField(const OptionalField& o)
: type(o.getType()) {
  strncpy(tag,o.getTag(),2);
}

const char *OptionalField::getTag() const {
  return tag;
}

char OptionalField::getType() const {
  return type;
}

void OptionalField::saveInBamRecord(kstring_t *str) {
  kputsn_(getTag(),2,str);
  kputc_(type, str);
  saveValueInBamRecord(str);
}

std::ostream& operator<<(std::ostream& os, const OptionalField& of)
{
  os << of.getTag()[0] << of.getTag()[1] << ":" << of.getType() << ":";
  of.printValue(os);
  return os;
}

CharOptionalField::CharOptionalField(const char tag[2], char v) 
: OptionalField(tag,'A'), value(v) {}

CharOptionalField::CharOptionalField(const CharOptionalField& o)
: OptionalField(o),value(o.getValue()) {}

CharOptionalField::~CharOptionalField() {}

OptionalField* CharOptionalField::clone() {
  return new CharOptionalField(*this);
}

char CharOptionalField::getValue() const {
  return value;
}

void CharOptionalField::saveValueInBamRecord(kstring_t *str) {
  uint8_t y = value;
  kputsn_(&y, 1, str);
}

std::ostream& CharOptionalField::printValue(std::ostream& os) const {
  os << value;
  return os;
}

IntOptionalField::IntOptionalField(const char tag[2], int v) 
: OptionalField(tag,'i'), value(v) {}


IntOptionalField::IntOptionalField(const char tag[2], uint v) 
: OptionalField(tag,'i'), value(v) {}

IntOptionalField::IntOptionalField(const IntOptionalField& o)
: OptionalField(o),value(o.getValue()) {}

IntOptionalField::~IntOptionalField() {}

OptionalField* IntOptionalField::clone() {
  return new IntOptionalField(*this);
}

int IntOptionalField::getValue() const {
  return value;
}

void IntOptionalField::saveValueInBamRecord(kstring_t *str) {
  uint32_t y = value;
  kputsn_(&y, 4, str);
}

std::ostream& IntOptionalField::printValue(std::ostream& os) const {
  os << value;
  return os;
}

StringOptionalField::StringOptionalField(const char tag[2], string v) 
: OptionalField(tag,'Z'), value(v) {}

StringOptionalField::StringOptionalField(const char tag[2], char* v) 
: OptionalField(tag,'Z'), value(v) {}

StringOptionalField::StringOptionalField(const StringOptionalField& o) 
: OptionalField(o), value(o.getValue()) {}

StringOptionalField::~StringOptionalField() {}

OptionalField* StringOptionalField::clone() {
  return new StringOptionalField(*this);
}

const string& StringOptionalField::getValue() const {
  return value;
}

void StringOptionalField::saveValueInBamRecord(kstring_t *str) {
  kputsn_(value.c_str(), value.length()+1, str);
}

std::ostream& StringOptionalField::printValue(std::ostream& os) const {
  os << value;
  return os;
}


SamLine::SamLine() {
  // init SAM fields with default values
  rname = "*";
  qname = "*";
  rnext = "*";
  seq = "*";
  qual = "*";
  rid = -1;
  ridnext = -1;
  flag = 0;
  pos = 0;
  mapQ = 255;
  //mapQ = 0;
  pnext = 0;
  tlen = 0;
  //is_mapped = true;
}

SamLine::~SamLine() {
  for (OF_map::iterator it=optionalFields.begin(); it!=optionalFields.end(); ++it) {
    delete it->second;
  }
}

SamLine::SamLine(const SamLine& s) 
: rname(s.getRname()),qname(s.getQname()),rnext(s.getRnext()),
  seq(s.getSeq()),qual(s.getQual()),rid(s.getRid()),ridnext(s.getRidnext()),
  flag(s.getFlag()),pos(s.getPos()),mapQ(s.getMapQ()),pnext(s.getPnext()),
  tlen(s.getTlen()),cigar(s.getCigar())
{
  OF_map tmp = s.getOptionalFields();
  for (OF_map::iterator it=tmp.begin(); it!=tmp.end(); ++it) {
    optionalFields[it->first] = it->second->clone();
  }
}

// Query template NAME
void SamLine::setQname(const string& qname) {
  this->qname = qname;
}
const string& SamLine::getQname() const {
  return qname;
}

// bitwize FLAG
void SamLine::setFlag(uint flag){
  this->flag = flag;
}
uint SamLine::getFlag() const {
  return flag;
}

// generic flag bit methods
bool SamLine::isFlagBitSet(uint bit) const {
  return flag & bit;
}
void SamLine::setFlagBit(uint bit) {
  flag = flag | bit;
}
void SamLine::unsetFlagBit(uint bit) {
  flag = flag^bit;
}

// Flag operations
// bit 1
bool SamLine::isTemplateHavingMultipleSegments() const {
  return isFlagBitSet(1);
}
void SamLine::setTemplateHavingMultipleSegments() {
  setFlagBit(1);
}
void SamLine::unsetTemplateHavingMultipleSegments() {
  unsetFlagBit(1);
}

// bit 2
bool SamLine::isEachSegmentsMapped() const {
  return isFlagBitSet(2);
}
void SamLine::setEachSegmentsMapped() {
  setFlagBit(2);
}
void SamLine::unsetEachSegmentsMapped() {
  unsetFlagBit(2);
}

// bit 4
bool SamLine::isSegmentUnmapped() const {
  return isFlagBitSet(4);
}
void SamLine::setSegmentUnmapped() {
  setFlagBit(4);
  mapQ = 0;
}
void SamLine::unsetSegmentUnmapped() {
  unsetFlagBit(4);
}

// bit 8
bool SamLine::isNextSegmentUnmapped() const {
  return isFlagBitSet(8);
}
void SamLine::setNextSegmentUnmapped() {
  setFlagBit(8);
}
void SamLine::unsetNextSegmentUnmapped() {
  unsetFlagBit(8);
}

// bit 16
bool SamLine::isSeqReverseComplemented() const {
  return isFlagBitSet(16);
}
void SamLine::setSeqReverseComplemented() {
  setFlagBit(16);
}
void SamLine::unsetSeqReverseComplemented() {
  unsetFlagBit(16);
}

// bit 32
bool SamLine::isNextSeqReverseComplemented() const {
  return isFlagBitSet(32);
}
void SamLine::setNextSeqReverseComplemented() {
  setFlagBit(32);
}
void SamLine::unsetNextSeqReverseComplemented() {
  unsetFlagBit(32);
}

// bit 64
bool SamLine::isFirstSegmentInTheTemplate() const {
  return isFlagBitSet(64);
}
void SamLine::setFirstSegmentInTheTemplate() {
  setFlagBit(64);
}
void SamLine::unsetFirstSegmentInTheTemplate() {
  unsetFlagBit(64);
}

// bit 128
bool SamLine::isLastSegmentInTheTemplate() const {
  return isFlagBitSet(128);
}
void SamLine::setLastSegmentInTheTemplate() {
  setFlagBit(128);
}
void SamLine::unsetLastSegmentInTheTemplate() {
  unsetFlagBit(128);
}

// bit 256
bool SamLine::isSecondaryAlignement() const {
  return isFlagBitSet(256);
}
void SamLine::setSecondaryAlignement() {
  setFlagBit(256);
}
void SamLine::unsetSecondaryAlignement() {
  unsetFlagBit(256);
}

// bit 512
bool SamLine::isFailingQualityControl() const {
  return isFlagBitSet(512);
}
void SamLine::setFailingQualityControl() {
  setFlagBit(512);
}
void SamLine::unsetFailingQualityControl() {
  unsetFlagBit(512);
}

// bit 1024
bool SamLine::isPCRDuplicated() const {
  return isFlagBitSet(1024);
}
void SamLine::setPCRDuplicated() {
  setFlagBit(1024);
}
void SamLine::unsetPCRDuplicated() {
  unsetFlagBit(1024);
}

// bit 2048
bool SamLine::isChimericAlignement() const {
  return isFlagBitSet(2048);
}
void SamLine::setChimericAlignement() {
  setFlagBit(2048);
}
void SamLine::unsetChimericAlignement() {
  unsetFlagBit(2048);
}

// Reference sequence NAME
void SamLine::setRname(const string& rname, const int rid) {
  this->rname = rname;
  this->rid = rid;
}
void SamLine::setUnknownRname() {
  this->rname = "*";
  this->rid = -1;
}
const string& SamLine::getRname() const {
  return rname;
}

int SamLine::getRid() const {
  return rid;
}

// 1-based leftmost mapping POSition
void SamLine::setPos(uint pos) {
  this->pos = pos;
}
void SamLine::setUnknownPos() {
  this->pos = 0;
}
uint SamLine::getPos() const {
  return pos;
}

// MAPping Quality
void SamLine::setMapQ(uint mapQ) {
  this->mapQ = mapQ;
}
uint SamLine::getMapQ() const {
  return mapQ;
}

// CIGAR string
void SamLine::setCigar(const Cigar &cigar) {
  this->cigar = cigar;
}
const Cigar& SamLine::getCigar() const {
  return cigar;
}

// Ref. name of the mate/next segment
void SamLine::setRnext(const string& rnext, const int ridnext) {
  this->rnext = rnext;
  this->ridnext = ridnext;
}
const string& SamLine::getRnext() const{
  return rnext;
}
int SamLine::getRidnext() const{
  return ridnext;
}

// Position of the mate/next segment
void SamLine::setPnext(uint pnext) {
  this->pnext = pnext;
}
uint SamLine::getPnext() const{
  return pnext;
}

// observed Template LENgth
void SamLine::setTlen(int tlen) {
  this->tlen = tlen;
}
int SamLine::getTlen() const{
  return tlen;
}

// segment SEQuence
void SamLine::setSeq(const string& seq) {
  this->seq = seq;
}
const string& SamLine::getSeq() const {
  return seq;
}
void SamLine::reverseComplementeSeq() {
  string new_seq;
  for (int i = seq.length() - 1; i >= 0; i--) {
    new_seq.push_back(complementDNA(seq[i]));
  }
  setSeq(new_seq);
}

// ASCII of Phred-scaled base QUALity+33
void SamLine::setQual(const string& qual) {
  this->qual = qual;
}
const string& SamLine::getQual() const {
  return qual;
}
void SamLine::reverseQual() {
  string new_qual;
  for (int i = qual.length() - 1; i >= 0; i--) {
    new_qual.push_back(qual[i]);
  }
  setQual(new_qual);
}

// Optional fields
void SamLine::addOptionalField(const char tag[2], char val) {
  OptionalField *new_optional_field = new CharOptionalField(tag,val);
  map<string,OptionalField*>::iterator it = optionalFields.find(tag);
  if(it != optionalFields.end()) {
    delete it->second;
    optionalFields.erase(it);
  }
  optionalFields[tag] = new_optional_field;
}
void SamLine::addOptionalField(const char tag[2], int val) {
  OptionalField *new_optional_field = new IntOptionalField(tag,val);
  map<string,OptionalField*>::iterator it = optionalFields.find(tag);
  if(it != optionalFields.end()) {
    delete it->second;
    optionalFields.erase(it);
  }
  optionalFields[tag] = new_optional_field;
}
void SamLine::addOptionalField(const char tag[2], uint val) {
  OptionalField *new_optional_field = new IntOptionalField(tag,val);
  map<string,OptionalField*>::iterator it = optionalFields.find(tag);
  if(it != optionalFields.end()) {
    delete it->second;
    optionalFields.erase(it);
  }
  optionalFields[tag] = new_optional_field;
}
void SamLine::addOptionalField(const char *tag, const string& val) {
  OptionalField *new_optional_field = new StringOptionalField(tag,val);
  map<string,OptionalField*>::iterator it = optionalFields.find(tag);
  if(it != optionalFields.end()) {
    delete it->second;
    optionalFields.erase(it);
  }
  optionalFields[tag] = new_optional_field;
}
void SamLine::addOptionalField(const char *tag, const char* val) {
  OptionalField *new_optional_field = new StringOptionalField(tag,val);
  map<string,OptionalField*>::iterator it = optionalFields.find(tag);
  if(it != optionalFields.end()) {
    delete it->second;
    optionalFields.erase(it);
  }
  optionalFields[tag] = new_optional_field;
}

const OF_map& SamLine::getOptionalFields() const {
  return optionalFields;
}

//const string& SamLine::getOptionalField(const char *tag) {
//  return optionalFields[tag];
//}

bool SamLine::isOptionalFieldDefined(const char *tag) {
  return optionalFields.count(tag) > 0;
}

void SamLine::removeAllUserOptionalFields() {
  OF_map::iterator it=optionalFields.begin();
  while(it != optionalFields.end()) {
    if(it->first[0] == 'X' || it->first[0] == 'Y' || it->first[0] == 'Z') {
      delete it->second;
      optionalFields.erase(it++);
    } else {
      ++it;
    }
  }
}

// write the SAM line in the file
ostream &SamLine::writeLine(ostream &os){
  os << getQname()  << "\t"
     << getFlag()   << "\t"
     << getRname()  << "\t"
     << getPos()    << "\t"
     << getMapQ()   << "\t"
     << getCigar()  << "\t"
     << getRnext()  << "\t"
     << getPnext()  << "\t"
     << getTlen()   << "\t"
     << getSeq()    << "\t"
     << getQual();
  for (OF_map::iterator it=optionalFields.begin(); it!=optionalFields.end(); ++it)
        os << "\t" << *it->second;
  os << endl;
  return os; 
}

// write the BAM record in the file
int SamLine::writeBamRecord(samFile *out, const bam_hdr_t *h) {
#define _get_mem(type_t, _x, _s, _l) ks_resize((_s), (_s)->l + (_l)); *(_x) = (type_t*)((_s)->s + (_s)->l); (_s)->l += (_l)

  int i;
  uint nb_cigar_ops = getCigar().count();

  // First we create a bam record that will hold our alignement informations
  bam1_t *b = bam_init1();

  // We will put everything that is not on the bam core structure here
  kstring_t str;

  str.l = b->l_data = 0;
  str.s = (char*)b->data; str.m = b->m_data;

  // First we fill the bam core structure
  bam1_core_t *c = &b->core;
  memset(c, 0, 32);

  c->tid = getRid();
  c->pos = getPos() - 1;
  c->bin = bam_reg2bin(getPos(),getPos() + getCigar().getReferenceAlignementLength());
  c->qual = getMapQ();
  c->l_qname = getQname().length() + 1; 
  c->flag = getFlag();
  c->n_cigar = (size_t) nb_cigar_ops;
  c->l_qseq = getSeq().length() > 0 && getSeq()[0] != '*'? getSeq().length() : 0;
  c->mtid = getRidnext();
  c->mpos = getPnext() - 1;
  c->isize = getTlen();

  // Now lets fill the rest of the data
  /**
   * QNAME
   *
   * First the read name (Qname)
   */
  kputsn_(getQname().c_str(), getQname().length()+1, &str);

  /**
   * CIGAR
   *
   * Now we set the cigar
   */
  if(nb_cigar_ops > 0) {
    uint32_t *packed_cigar;
    //uint32_t *packed_cigar = new uint32_t[nb_cigar_ops];
    // We set memory spaces on str and point the cigar array to that memory space
    _get_mem(uint32_t, &packed_cigar, &str, nb_cigar_ops<<2);
    for(i = 0; i < (int)nb_cigar_ops; i++) {
      // We set the nb and we shift 4 bits to mask the uint with the cigar op char
      // encoded of these 4 bits
      packed_cigar[i] = getCigar().get(i).nb << BAM_CIGAR_SHIFT;
      //                             0123456789
      // This is the current BAM_TAB MIDNSHP=XB
      // TODO use bam_hdr_t in order to get a cigar_tab that is a more nice way to do this
      uint cigar_op_mask = cigarop2int(getCigar().get(i).type);
      packed_cigar[i] |= cigar_op_mask;
      //c->bin = hts_reg2bin(c->pos, c->pos + bam_cigar2rlen(nb_cigar_ops,packed_cigar), 14, 5);
    }
    //kputsn_(packed_cigar, nb_cigar_ops, &str)
    //delete [] packed_cigar;
  } //else {
  //  c->bin = hts_reg2bin(c->pos, c->pos + 1, 14, 5);
  //}
  
  /**
   * SEQUENCE
   */
  if(c->l_qseq > 0) {
    uint8_t *packed_seq;
    i = (c->l_qseq + 1) >> 1;
    _get_mem(uint8_t, &packed_seq, &str, i);
    memset(packed_seq, 0, i);
    for (i = 0; i < c->l_qseq; ++i)
      packed_seq[i>>1] |= seq_nt16_table[(int)getSeq()[i]] << ((~i&1)<<2);
  }

  /**
   * QUALITY
   */
  uint8_t *packed_qual;
  _get_mem(uint8_t, &packed_qual, &str, c->l_qseq);
  if (strcmp(getQual().c_str(), "*")) {
    for (i = 0; i < c->l_qseq; ++i) packed_qual[i] = getQual()[i] - 33;
  } else memset(packed_qual, 0xff, c->l_qseq);

  /**
   * OPTIONAL FIELDS
   */
  for (OF_map::iterator it=optionalFields.begin(); it!=optionalFields.end(); ++it) {
    it->second->saveInBamRecord(&str);
  }

  // We assing kstring structure values to the bam entry
  b->data = (uint8_t*)str.s; b->l_data = str.l; b->m_data = str.m;

  int r = sam_write1(out,h,b);

  bam_destroy1(b);

  return r;
}

