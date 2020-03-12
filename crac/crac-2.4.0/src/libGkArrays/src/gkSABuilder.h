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

#ifndef GKSABUILDER_H
#define GKSABUILDER_H

#include "gkArraysTypes.h"
#include "gkArrays.h"
#include "solArray.h"
#include <semaphore.h>

using namespace std;

namespace gkarrays {

struct qs
{
  SolArray *gkSA;
  SolArray *values;
  int64_t start;
  int64_t end;
  int threshold_sort;
  sem_t *s_threads;
};

/**
 * Initializes a pointer to a qs data structure.
 */
 void init_qs_ptr(struct qs *data, SolArray *sa, SolArray *val, int64_t start,
                  int64_t end, int threshold_sort, sem_t *s_threads);

struct subGkSA
{
  //general features
  unsigned char *cr;
  SolArray *gkSA;
  SolArray *values;
  int m;
  int k;
  gkArrays *gk;
  //subGkSA features
  uintSA start;
  uintSA end;
  uintSA offset;
};

class GkSABuilder {
  
 private:
  SolArray *gkSA;
  
 public:
 
  GkSABuilder();
  ~GkSABuilder();

  /**
   * Build a GkSA
   * @param values: array storing the integer values of the k-mers
   * @param cr: the concatenation of reads
   * @param length: length of the read concatenation
   * @param k: length of the k-mers
   * @param gk: gkArrays to be used for computing gkIFA
   * @param threshold_sort: threshold below which an insertion sort is performed
   *        rather than a quicksort
   * @param type: the type of SolArray to build
   * @param nb_threads: the number of threads that can be used
   */
  void build(SolArray *values, unsigned char* cr, uintSA length, int k, 
             gkArrays *gk, int threshold_sort, 
             array_type type, uint nb_threads);  

  /**
   * @return the last GkSA that has been built.
   */
  SolArray *getGkSA();
  
};

/**
 * Sort the elements in sa according to the values stored in values. Both arrays are sorted at the same time
 * (actually the operations performed on values are also performed on sa).
 * @param sa: the array to be sorted, its values are not used for sorting
 * @param values: the array whose values are used for sorting
 * @param start: where to start the sort in the input arrays.
 * @param end: where the sort should end in the input arrays.
 * @param threshold_sort: Threshold below which an insertion sort is used rather than a quicksort.
 *                        Insertion sort is quicker for small arrays. TODO: replace insertion sort by heapsort.
 * @param s_threads: the semaphore for the threads
 */
 void quickSort(SolArray *sa, SolArray values[], int64_t start, int64_t end,int threshold_sort,sem_t *s_threads, bool is_threaded=false);

/**
 * The function called by each thread when they're created.
 * This function then direcly calls the function quickSort with the values
 * stored in the data parameter.
 * @return NULL
 */
void* quickSortT(void* data);

/**
 * Partition step in quicksort algorithm
 * @param sa: the arrays to be sorted, but its values are not considered for sorting
 * @param values: the arrays to be sorted, the exact same operations are performed on sa
 * @param start: start pos.
 * @param end: end pos.
 * @param p: Pivot value.
 * @return the pivot's position in the array values.
 */
int64_t partition(SolArray *sa,SolArray *values,uintSA start, uintSA end,uintSA p);

/**
 * Performs an insertion sort on the given arrays, at the given start and end
 * positions
 */
void insertion_sort(SolArray *sa , SolArray *values, int64_t &start, int64_t &end);

//return the value at pos p on the strand forward
//the last parameter is just for the compatibility with valns
uintSA vals(uintSA p,unsigned char dna[],int k,int64_t *last,uintSA *valfwd, uintSA *val=NULL);

  //return the minumum value of the k-mer at pos p between strand rev and stran fwd
uintSA valns(uintSA p, unsigned char dna[],int k,int64_t *last,uintSA *valfwd,uintSA *valrev);

// It is a function that is common to calcs and calcns
// It performs all the computations for putting the values in GkSA, and dedicates
// the GkSA value calculation to the function given in parameter
 void *calc_gkSASubPart(void *data);
 
 inline void swap(SolArray *sa,SolArray *values, uintSA &a, uintSA &b){
   if(a!=b){
     uintSA temp = values->get(a);
     uintSA temp2= sa->get(a);
     values->set(a, values->get(b));
     values->set(b, temp);
     sa->set(a,sa->get(b));
     sa->set(b,temp2);
   }
 }

}
#endif
