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

#include "GenomeInfo.h"
#include <cassert>
#include <fstream>

  /**
   * Array containing the sequence names
   */
  uchar **chr_names;
  /**
   * Number of sequences
   */
  ulong nb_chr;
  /**
   * Array containing the length (in nt) of each sequence
   */
  ulong *length_chr;
  /**
   * Array containing the cumulative length of each sequence.
   * First cell's value is 0
   */
  ulong *total_length;


  /**
   * Read the file whose name is stored in <basename>.conf.
   * Retrieves the information from the file and store them 
   * in the class.
   * @param basename: name of the file without the .conf extension
   */
GenomeInfo::GenomeInfo(char *basename) {
  char *genomeConf = new char[strlen(basename)+6];
  fstream file;
  sprintf(genomeConf, "%s.conf", basename);
     
  file.open(genomeConf, ios::in);
  if (! file.is_open()) {
    cerr << "Cannot open chromosomes file: " << genomeConf << endl;
    exit(2);
  }

  delete [] genomeConf;

  // Store the number of sequences
  file >> skipws >> nb_chr;

  // Allocate memory for the arrays
  chr_names = new uchar*[nb_chr];
  total_length = new ulong[nb_chr+1];
  length_chr = new ulong[nb_chr];
  total_length[0] = 0;

  int i = 0;
  file.width(MAX_LENGTH_CHR_NAME);
  while (! file.eof()) {
    chr_names[i] = new uchar[MAX_LENGTH_CHR_NAME];
    file >> skipws >> chr_names[i];
    file >> skipws >> length_chr[i];
    total_length[i+1] = total_length[i]+ length_chr[i];
    i++;
  }
  file.close();
}

GenomeInfo::GenomeInfo(uchar **chr_names, ulong nb_chr, ulong *length_chr):
  chr_names(chr_names), nb_chr(nb_chr), length_chr(length_chr) {
  total_length = new ulong[nb_chr+1];
  total_length[0] = 0;
  for (ulong i = 0; i < nb_chr; i++) {
    total_length[i+1] = total_length[i] + length_chr[i];
  }
}

GenomeInfo::~GenomeInfo() {
  for (uint i = 0; i < nb_chr; i++)
    delete [] chr_names[i];
  delete [] chr_names;
  delete [] total_length;
  delete [] length_chr;
}

uchar *GenomeInfo::getChrName(ulong i) {
  assert(i >= 0 && i < getNbChr());
  
  return chr_names[i];
}

ulong GenomeInfo::getChrLength(ulong i) {
  assert(i >= 0 && i < getNbChr());
  
  return length_chr[i];
}

ulong GenomeInfo::getGenomeLength() {
  return total_length[getNbChr()];
}

ulong GenomeInfo::getGenomePosition(uchar *chr, ulong position) {
  return getGenomePosition(getNumFromName(chr), position);
}

ulong GenomeInfo::getGenomePosition(ulong id, ulong position) {
  if (id < getNbChr() && position < length_chr[id]) {
    return total_length[id] + position;
  }
  return ~0;
}

ulong GenomeInfo::getNbChr() {
  return nb_chr;
}

ulong GenomeInfo::getNumFromName(uchar *name) {
  ulong num_chr = 0;
  while (num_chr < getNbChr()
         && strcmp((char *)chr_names[num_chr], (char *)name) != 0){
    num_chr++;
  }
  if (num_chr >= getNbChr()) {
    return getNbChr();
  }
  return num_chr;
}

ChrPosition *GenomeInfo::getChrPositionFromPosition(ulong pos) {
  pair<ulong, ulong> tmp_pair = getNumFromPosition(pos);
  int chr_id = (int) tmp_pair.first;
  uchar *chr_name;
  if (tmp_pair.first >= getNbChr()) {
    chr_name = NULL;
    chr_id = -1;
    cerr << "Erorr in chr position" << endl;
  } else {
    chr_name = getChrName(tmp_pair.first);
  }
  return new ChrPosition(tmp_pair.second, (char *)chr_name, chr_id);
}

pair<const ulong, const ulong> GenomeInfo::getNumFromPosition(ulong pos) {
  ulong num_chr = 0;
  
  while (num_chr < getNbChr() && pos >= total_length[num_chr+1]) {
    num_chr++;
  }
  return pair<const ulong, const ulong>(num_chr, pos - total_length[num_chr]);
}

bool GenomeInfo::isValidPosition(ulong pos) {
  return pos < getGenomeLength();
}


void GenomeInfo::save(char *name, const char *ext) {
  char *confName = new char[strlen(name) + strlen(ext)];
  ofstream conf;
  sprintf(confName, "%s%s", name, ext);
  conf.open(confName, ios_base::out | ios_base::binary);
  if (! conf.is_open()) {
    cerr << "Can't write file " << confName << endl;
    exit(3);
  }
  delete [] confName;
  
  // Writing configuration file
  conf << *this;
  conf.close();

}

ostream &operator<<(ostream &os, GenomeInfo &g) {
  os << g.getNbChr();
  for (uint i=0; i < g.getNbChr(); i++) {
    os << endl << g.getChrName(i) 
       << endl << g.getChrLength(i);
  }
  return os;
}

/* DEPREACTED METHODS */

pair<const uchar *,const ulong> GenomeInfo::getChrNameFromPosition(ulong pos) {
  pair<ulong, ulong> tmp_pair = getNumFromPosition(pos);
  uchar *chr_name;

  if (tmp_pair.first >= getNbChr()) {
    chr_name = NULL;
  } else {
    chr_name = getChrName(tmp_pair.first);
  }
  return (pair<const uchar *, const ulong>(chr_name, tmp_pair.second));
}
