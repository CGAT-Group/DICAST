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

#include "CandidatBreak.h"
#include "SupportBreak.h"
#include <math.h>

using namespace std;

CandidatBreak::CandidatBreak(uint pos_start, uint loc_start, int strand_start, 
			     uint pos_end, uint loc_end, int strand_end,
			     SupportBreak *sb, Parameters *p 
			     ):pos_start_break(pos_start), loc_start_break(loc_start), strand_start_break(strand_start), pos_end_break(pos_end), loc_end_break(loc_end), strand_end_break(strand_end), supportBreak(sb), parameters(p) {
  is_single_consistent = false;
  single_distance = ~0;
  is_duplicated = false;
} 

CandidatBreak::CandidatBreak(const CandidatBreak &c):pos_start_break(c.pos_start_break), loc_start_break(c.loc_start_break), strand_start_break(c.strand_start_break), pos_end_break(c.pos_end_break), loc_end_break(c.loc_end_break), strand_end_break(c.strand_end_break), supportBreak(c.supportBreak), parameters(c.parameters) {
  is_single_consistent = c.is_single_consistent;
  is_duplicated = c.is_duplicated;
  single_distance = c.single_distance;
}

CandidatBreak::~CandidatBreak(){}


void CandidatBreak::checkSingleCorrespondence(uint pos_single, uint loc_single, int strand_single){
  // case of loc_single does not exist
  if (strand_single == 0){
    is_single_consistent = false;
    single_distance = ~0;
  }else{
    ulong chr_single = getGenome()->getIdChromosome(loc_single);
    // if pos_single is before break, checking of strand and chr
    if (strand_single == strand_start_break 
	&& getPosEndBreak() > pos_single
	&&  chr_single == getChrIdStartBreak()){
      // compute the distance between loc_single and loc_start_break
      if ((strand_single == 1 && pos_single <= getPosStartBreak()) 
	  || (strand_single == -1 && pos_single >= getPosStartBreak()))
	single_distance = (long long int) getLocStartBreak() - loc_single;
      else
	single_distance = (long long int) loc_single - getLocStartBreak();
      // is_single consistent iff single_distance >= 0 && single_distance <= parameters->max_splice_length
      if (single_distance >= 0 && single_distance <= parameters->max_splice_length) 
	is_single_consistent = true;
      else{
	is_single_consistent = false;
	single_distance = ~0;
      }
    }
    // if pos_single is after break, checking of strand and chr
    else if (strand_single == strand_end_break 
	     && getPosStartBreak() < pos_single
	     && chr_single == getChrIdEndBreak()){
      // compute the distance between loc_single and loc_start_break
      if ((strand_single == 1 && pos_single >= getPosEndBreak()) 
	  || (strand_single == -1 && pos_single <= getPosEndBreak()))
	single_distance = (long long int) loc_single - getLocEndBreak();
      else
	single_distance = (long long int) getLocEndBreak() - loc_single;
      // is_single consistent iff single_distance >= 0
      if (single_distance >= 0) 
	is_single_consistent = true;
      else{
	is_single_consistent = false;
	single_distance = ~0;
      }
    }
    // other cases: not single consistency
    else{
      is_single_consistent = false;
      single_distance = ~0;
    }
  }
}


const uchar *CandidatBreak::getChrEndBreak(){
  return getGenome()->getChromosome(loc_end_break);
}


ChrPosition *CandidatBreak::getChrPositionEndBreak(){
  return getGenome()->getChrPos(loc_end_break,strand_end_break);
}


ChrPosition *CandidatBreak::getChrPositionStartBreak(){
  return getGenome()->getChrPos(loc_start_break,strand_start_break);
}

ulong CandidatBreak::getChrIdEndBreak(){
  return getGenome()->getIdChromosome(loc_end_break);
}

ulong CandidatBreak::getChrIdStartBreak(){
  return getGenome()->getIdChromosome(loc_start_break);
}

const uchar *CandidatBreak::getChrStartBreak(){
    return getGenome()->getChromosome(loc_start_break);
}

LocateOnGenome *CandidatBreak::getGenome(){
  return supportBreak->getGenome();
}

gap_size_t CandidatBreak::getGenomeGapLength(){
  if (strand_start_break != strand_end_break
      || (getChrIdStartBreak() != getChrIdEndBreak()))
    return min(labs(loc_end_break - loc_start_break),
	       labs(loc_start_break - loc_end_break)) - 1;
  if (strand_end_break == -1)
    return (gap_size_t)loc_start_break - (gap_size_t)loc_end_break-1;
  return (gap_size_t)loc_end_break - (gap_size_t)loc_start_break-1;
}

uint CandidatBreak::getLength(){
  return supportBreak->getLength();
}

uint CandidatBreak::getLocEndBreak(){
  return loc_end_break;
}

uint CandidatBreak::getLocStartBreak(){
  return loc_start_break;
}

uint CandidatBreak::getPosEndBreak(){
  return pos_end_break;
}

uint CandidatBreak::getPosStartBreak(){
  return pos_start_break;
}

uint CandidatBreak::getReadBreakLength(){
  return (pos_end_break - pos_start_break + 1);
}

long long int CandidatBreak::getSingleDistance(){
  return single_distance;
}

int CandidatBreak::getStrandEndBreak(){
  return strand_end_break;
}

int CandidatBreak::getStrandStartBreak(){
  return strand_start_break;
}

bool CandidatBreak::hasNoEndBreak(){
  return (pos_end_break == (getLength() - 1)); 
}

bool CandidatBreak::hasNoStartBreak(){
  return (pos_start_break == 0);
}

bool CandidatBreak::isDuplicated(){
  return is_duplicated;
}

bool CandidatBreak::isNiceCandidat(){
  return (isSameChr() && isSameStrand()
	  && (
	      (getStrandStartBreak() == 1 && (getLocStartBreak() + supportBreak->getThreshold() - 1) < getLocEndBreak())
	      || (getStrandStartBreak() == -1 && (getLocEndBreak() + supportBreak->getThreshold() - 1) < getLocStartBreak())
	      )
	  //&& (getGenomeGapLength() <= parameters->max_splice_length)
	  );
}

bool CandidatBreak::isSameChr(){
  return (getChrIdStartBreak() == getChrIdEndBreak());
}

bool CandidatBreak::isSameStrand(){
  return (strand_start_break == strand_end_break);
}

bool CandidatBreak::isSingleConsistent(){
  return is_single_consistent;
}
 
void CandidatBreak::setDuplicated(bool flag){
  is_duplicated = flag;
}

void CandidatBreak::setLocEndBreak(uint loc_end){
  loc_end_break = loc_end;  
}

void CandidatBreak::setLocStartBreak(uint loc_start){
  loc_start_break = loc_start;
}

void CandidatBreak::setPosEndBreak(uint pos_end){
  pos_end_break = pos_end;
}

void CandidatBreak::setPosStartBreak(uint pos_start){
  pos_start_break = pos_start;
}

void CandidatBreak::setSingleConsistent(bool flag){
  is_single_consistent = flag;
}

void CandidatBreak::setStrandEndBreak(int flag){
  strand_end_break = flag;
}

void CandidatBreak::setStrandStartBreak(int flag){
  strand_start_break = flag;
}
