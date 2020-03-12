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

#include "Cigar.h"
#include <cassert>

Cigar::Cigar(bool use_x_in_cigar):nb_M(0), nb_I(0), nb_D(0), nb_S(0), nb_P(0), nb_H(0), nb_EQ(0), nb_X(0), nb_N(0), nb_chimeras(0){
  this->use_x_in_cigar = use_x_in_cigar;
}


void Cigar::append(uint nb, char op) {
  assert(op == 'M' || op == 'I' || op == 'D' || op == 'S' || op == '='
         || op == 'P' || op == 'H' || op == 'X' || op == 'N');
  cigar_type ct = {nb, op};
  append(ct);
}

void Cigar::append(cigar_type cigar) {
  update_counters(cigar);
  this->cigar.push_back(cigar);
}

void Cigar::filter() {
  char old_type = 0;
  uint old_nb = 0;
  vector<cigar_type>::iterator it = cigar.begin();
  while (it < cigar.end()) {
    // We remove X operators if the option is not activated
    if (it->type == 'X' && !use_x_in_cigar) {
      it->type = 'M';
    }
    // Two identical consecutive types
    if (it->type == old_type 
        && (it->type != 'N' || (it->nb > 0 && old_nb > 0))) {
      it->nb += old_nb;
      it = cigar.erase(it-1);
    } else if (it->type != 'N' && it->nb == 0) {
      it = cigar.erase(it)-1;
    }

    old_type = it->type;
    old_nb = it->nb;
    it++;
  }
}

void Cigar::prepend(uint nb, char op) {
  assert(op == 'M' || op == 'I' || op == 'D' || op == 'S' || op == '='
         || op == 'P' || op == 'H' || op == 'X' || op == 'N');
  cigar_type ct = {nb, op};
  prepend(ct);
}

void Cigar::prepend(cigar_type cigar) {
  update_counters(cigar);
  this->cigar.insert(this->cigar.begin(), cigar);
}

  /* Queries */

uint Cigar::count() const {
  return cigar.size();
}

const cigar_type& Cigar::get(uint i) const{
  return cigar[i];
}

uint Cigar::getQueryAlignementLength() const {
  return nb_S + nb_M + nb_I + nb_EQ + nb_X;
}

uint Cigar::getReferenceAlignementLength() const {
  return nb_M + nb_D + nb_N + nb_EQ + nb_X;
}

bool Cigar::isReferenceBasedElement(uint i) const {
  return cigar[i].type == 'M' || cigar[i].type == '=' || cigar[i].type == 'X' || cigar[i].type == 'D' || cigar[i].type == 'N';
}

/**
 * @return true is the operator assigned to the ith element correspond to 
 * the read sequence (ie. M, =, X, I, S, H)
 */
bool Cigar::isQueryBasedElement(uint i) const {
  return cigar[i].type == 'M' || cigar[i].type == '=' || cigar[i].type == 'X' || cigar[i].type == 'I' || cigar[i].type == 'S' || cigar[i].type == 'H';
}

uint Cigar::getNbDeletions() const {
  return nb_D + nb_N;
}

uint Cigar::getNbPositionsRead() const {
  return nb_M + nb_EQ + nb_X + nb_I + nb_S;
}

uint Cigar::getEditDistance() const {
  return nb_X + nb_I + nb_D + nb_N;
}

uint Cigar::getNbChimeras() const {
  return nb_chimeras;
}

void Cigar::reverse() {
  vector<cigar_type> reverse_cigar;
  for(int i = cigar.size(); i > 0; i--) {
    reverse_cigar.push_back(cigar[i-1]);
  }
  cigar = reverse_cigar;
}

cigar_type &Cigar::operator[](uint i) {
  return cigar[i];
}


ostream &operator<<(ostream &os, const Cigar &c) {
  if (c.count() > 0) {
    for (uint i=0; i < c.count(); i++) {
      os << c.get(i);
    }
  } else {
    os << "*";
  }
  return os;
}

void Cigar::update_counters(const cigar_type &cigar) {
  if(cigar.type == 'M') {
    nb_M += cigar.nb;
  } else if(cigar.type == '=') {
    nb_EQ += cigar.nb;
  } else if(cigar.type == 'X') {
    nb_X += cigar.nb;
  } else if(cigar.type == 'I') {
    nb_I += cigar.nb;
  } else if(cigar.type == 'D') {
    nb_D += cigar.nb;
  } else if(cigar.type == 'N') {
    if(cigar.nb == 0) {
      nb_chimeras += 1;
    } else {
      nb_N += cigar.nb;
    }
  } else if(cigar.type == 'S') {
    nb_S += cigar.nb;
  } else if(cigar.type == 'H') {
    nb_H += cigar.nb;
  }
}
