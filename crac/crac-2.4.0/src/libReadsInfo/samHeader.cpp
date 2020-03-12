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
#include "samHeader.h"
#include <sstream>
#include <sys/time.h> // time()

/*************************
 * SamReference methods
 *************************/

SamReference::SamReference(const char *n, uint l)
  : name(n),length(l) {}

const string& SamReference::getName() const {
  return name;
}

uint SamReference::getLength() const {
  return length;
}


/*************************
 *   SamHeader methods
 *************************/

SamHeader::SamHeader() {
  read_group_id = "1";
}

void SamHeader::addReference(const char *name, uint length) {
  references.push_back(SamReference(name,length)); 
}

void SamHeader::addArgument(const char *arg) {
  arguments.push_back(string(arg));
}

void SamHeader::setVersionNumber(const char *vn) {
  version_number = vn;
}

void SamHeader::setSampleName(const char* sn) {
  sample_name = sn;
}


void SamHeader::setReadGroupId(const char* rg) {
  read_group_id = rg;
}

/*
 * Create a new bam_hdr_t object and return it
 */
bam_hdr_t *SamHeader::getBamHeader() {

  bam_hdr_t *h = bam_hdr_init();

  ostringstream ss;
  // Set Header text
  //ss << "@HD\tVN:" << hts_version() << "\tSO:unsorted" << endl;
  ss << "@HD\tVN:" << "1.5" << "\tSO:unsorted" << endl;
  ss << "@RG\tID:" << read_group_id;
  if(!sample_name.empty())
    ss << "\tSM:" << sample_name;
  struct tm tm;
  char tmp[120]; // large enough buffer
  time_t t;
  time(&t);
  localtime_r(&t,&tm); 
  sprintf(tmp,"%i-%02i-%02iT%02i:%02i:%02i", tm.tm_year + 1900, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  ss << "\tDT:" << tmp << endl;
  ss << "@PG\tID:1\tPN:crac\tVN:" << version_number;
  if(arguments.size() > 0) {
    ss << "\tCL:";
    int k = 0;
    for (std::vector<string>::iterator it = arguments.begin() ; it != arguments.end(); ++it) {
      if(k > 0)
        ss << " ";
      ss << (*it);
      k++;
    }
  }
  ss << endl;

  string text = ss.str();
  h->text = new char[text.length()+1];
  strncpy(h->text,text.c_str(),text.length());
  h->text[text.length()] = '\0';
  h->l_text = text.length();

  // Set Header references
  h->n_targets = references.size();
  h->target_len = (uint32_t*)malloc(sizeof(uint32_t) * h->n_targets);
  h->target_name = (char**)malloc(sizeof(char*) * h->n_targets);
  int i = 0;
  for (std::vector<SamReference>::iterator it = references.begin() ; it != references.end(); ++it) {
    h->target_name[i] = (char*)malloc(sizeof(char) * (it->getName().length()+1));
    strcpy(h->target_name[i],it->getName().c_str());
    h->target_len[i] = it->getLength();
    i++;
  }
  return h;
}

/*
 * Similar as getBamHeader but also write the header
 * to the samFile output before returning it
 */
bam_hdr_t *SamHeader::writeBamHeader(samFile *out) {
 bam_hdr_t *h = getBamHeader();
 sam_hdr_write(out,h);
 return h;
}
