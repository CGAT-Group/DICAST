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

#include "utils.h"

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

namespace gkarrays {

  void dnaFiltration(char *dna, uintSA length) {
    for (uintSA i = 0; i < length; i++) {
      dna[i] = intToNuc(convNuc(dna[i]));
    }
  }

  uintSA DNAtoInt(char *dna, uintSA dna_length){
    uintSA dna_int = 0;
    for (uintSA i=0 ; i< dna_length ; i++){
      dna_int <<= 2;      
      dna_int |= convNuc(dna[i]);
    }
    return dna_int;
  }

  int comparUint(const void *a1, const void* a2){
    return (int) (* (uint *)a1 - * (uint *) a2);
  }

  int comparUintSA(const void *a1, const void* a2){
    return (int) (* (uintSA *)a1 - * (uintSA *) a2);
  }

  uint convNuc(char nuc){
    switch (nuc){
    case 'a' : case 'A' : return 0 ;
    case 'c' : case 'C' : return 1 ;
    case 'g' : case 'G' : return 2 ;
    case 't' : case 'T' : return 3 ;
    default : return 0;
    }  
    return 0;
  }

  char intToNuc(uint c) {
    switch(c) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    default: cerr << "Invalid number (0-3): " << c <<endl; exit(5);
    }
  }

  void intToDNA(uint64 code, uint dna_length, char *dna) {
    uint64 mask = 3;
    for (uint i=0; i < dna_length; i++) {
      dna[dna_length-i-1] = intToNuc(code & mask);
      code >>=2;
    }
  }

  uint64 factorsToInt(uintSA pos, uchar *dna, uint length) {
    assert(length >= 4 && length <= 32);

    uintSA posFirstFactor = 4 - (pos % 4);
    uintSA numFirstFactor = pos / 4;
    uintSA posLastFactor = (length - posFirstFactor) % 4;
    uint64 factor = dna[numFirstFactor] & ((1 << (posFirstFactor*2)) - 1);
    for(uintSA i = 1 ; i <= (length-posFirstFactor)/4 ; i++){
      factor <<= 8;
      factor |= dna[numFirstFactor+i];
    }
    factor <<= posLastFactor*2;

    factor |= dna[numFirstFactor+(length-posFirstFactor)/4+1] >> (8-2*posLastFactor) & ((1 << (posLastFactor*2)) - 1);
    
    return factor;
  }

  //identique a factorsToInt pour le cas stranded
  uint64 factorsToIntNoStranded(uintSA pos, uchar *dna, uint length) {
    uint64 f1=factorsToInt(pos,dna,length);
    return min(f1, intRevcomp(f1, length));
  }

  uint64 intRevcomp(uint64 factor, uint length) {
    uint64 mask;
    if (length == 32)
      mask = ~0;
    else
      mask =  ((uint64) 1 << (2*length)) - 1;

    factor ^= mask;

    uint64 mask_lsb;
    // Corresponds to the rightmost nucleotide
    mask_lsb = 3;
    uint64 shift = 0;
    uint64 result=0;
    for(uint j(0);j<length;j++){
      result <<= 2;
      // get the leftmost nucleotide and put it at the end
      result |= (factor & mask_lsb) >> shift;
      mask_lsb <<= 2;
      shift += 2;
    }

    return result;
  }
}
