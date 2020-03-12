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

#include <gkSABuilder.h>
#include "basic.h"
#include "utils.h"
#include "solArray.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <pthread.h>
#include <errno.h>
#include <gkCFPSBuilder.h>
#include <gkIFABuilder.h>
#include <const.h>
#ifdef HAVE_LIBPROGRESSBAR
#include <libProgressBar/progressBar.h>
uintSA PB_thread_master_start_pos = 0;
#include <unistd.h>
#endif

using namespace std;

namespace gkarrays {

  GkSABuilder::GkSABuilder():gkSA(NULL){}

  //When nb_kmers < threshold_sorting by insertion is more efficient and
  //quicksort is more efficient otherwise
  void GkSABuilder::build(SolArray *values, unsigned char *cr, uintSA cr_length, 
                          int k, gkArrays *gk, int threshold_sort, array_type type, uint nb_threads){

    uintSA nb_reads=gk->getNbTags();    
    uintSA nb_kmers=gk->getNbPposition(nb_reads); 

    // Creates the structure
    gkSA = new SolArray(nb_kmers, cr_length, type);
    
    // Building GkSA
    // gkSA is divided by nb_threads    
    pthread_t T[nb_threads];
    uintSA start=0;
    uintSA step=nb_reads/nb_threads;
    subGkSA sgkSA[nb_threads];
    void *fus;    

    // We just initialize SA with an increasing sequence of numbers (ie. after
    // that gkSA[i] == i) and values with the corresponding values. After
    // that, values[i] == the integer value of the k-mer starting at
    // P-position i
    for(uintSA i=0;i<nb_threads-1;i++){
      sgkSA[i].cr = cr;
      sgkSA[i].gkSA = gkSA;
      sgkSA[i].values = values;
      sgkSA[i].k = k;
      sgkSA[i].gk = gk;
      sgkSA[i].start = gk->getStartPosOfTagNum(start);
      sgkSA[i].end = gk->getEndPosOfTagNum(start+step);
      sgkSA[i].offset = gk->getStartQPosOfTagNum(i*step);
      pthread_create(&T[i], NULL, calc_gkSASubPart, (void *)&sgkSA[i]);
      start=start+step;
    }

    // Plus the final one.
    subGkSA s;
    s.cr = cr;
    s.gkSA = gkSA;
    s.values = values;
    s.k = k;
    s.gk = gk;
    s.start = gk->getStartPosOfTagNum(start);
    s.end = cr_length-1;
    s.offset = gk->getStartQPosOfTagNum((nb_threads-1)*step);
    pthread_create(&T[nb_threads-1], NULL, calc_gkSASubPart, (void*)&s);
    
    // Wait for them to finish their job
    for(uint i=0;i<nb_threads;i++){
      pthread_join(T[i], &fus);
    }

    sem_t thread_semaphore;     // Semaphore used for counting the number of
                                // threads launched and for avoiding launching
                                // more threads than needed.
    if (sem_init(&thread_semaphore, 0, nb_threads-1) == -1) {
      perror("Semaphore initialization for multi-threading\n");
      exit(2);
    }

#ifdef HAVE_LIBPROGRESSBAR
    clock_t start_sort = clock();
    int width = (gkArrays::GetNbColumns() ? gkArrays::GetNbColumns() : 100);
    if (gkArrays::ShowProgressBar()) {
      cerr << "Sorting SA       [This may take a while"
	   << setw(width - 54)
	   << setfill('.') << "]"
	   << setw(width - 1) << setfill('\b') << '\b' << flush;
    }
#endif
    quickSort(gkSA,values,0,nb_kmers-1,threshold_sort,&thread_semaphore, true);
    // The current thread is going to wait, so we say that we have another
    // thread free for computations
    sem_post(&thread_semaphore);

    int val = 0;
    useconds_t delay = 1000;
    // Wait for all the threads to finish
    while (val < (int)nb_threads) {
      sem_getvalue(&thread_semaphore, &val);
      usleep(delay);
      if (delay < 100000)
        delay *= 2;
    }
#ifdef HAVE_LIBPROGRESSBAR
    if (gkArrays::ShowProgressBar()) {
      double temps = (clock() - start_sort)/double(CLOCKS_PER_SEC);
      stringstream tps;
      tps << " 100 % (";
      if (temps < 60.0) {
        tps.precision((temps < 10.0)?2:1);
        tps.setf(stringstream::fixed,stringstream::floatfield);
        tps << temps << "s";
      } else {
        tps.precision(0);
        if (temps < 3600.0) {
	  tps << int(temps)/60 << "m"
	      << setfill('0') << setw(2) << int(temps)%60;
        } else {
          if (temps < 86400.0) {
	    tps << int(temps)/3600 << "h"
	        << setfill('0') << setw(2) << int(temps)%3600/60;
          } else {
	    tps << int(temps)/86400 << "d"
	        << setfill('0') << setw(2) << int(temps)%86400/3600;
          }
        }
      }
      tps << ")";
      cerr << "Sorting SA       ["
	   << setw(width - 19 - tps.str().length())
	   << setfill('=') << "]" << tps.str() << endl;
    }
#endif

  }

  GkSABuilder::~GkSABuilder() {}

  SolArray *GkSABuilder::getGkSA(){
    return gkSA;
  }

  //PRIVATE

  //tri rapide f est le seuil pour lancer le tri par insertion, r est le nombre de thread deja lancés et t est le nombre de thread max
  void quickSort(SolArray *sa, SolArray *values, int64_t start, int64_t end,
                 int threshold_sort, sem_t *s_threads, bool is_threaded){
    if(start>=end){
      return;
    }
    if((end-start)<threshold_sort){
      insertion_sort(sa,values,start,end);
      return;
    }
    // Choose the pivot
    uintSA pivot=start + (rand() % (end - start + 1));
    // Partition the array.
    int64_t m=partition(sa,values,start,end,pivot);
    
    pthread_t thread, thread2;
    qs *valeur = (qs*)malloc(sizeof(qs)), *valeur2 = (qs*)malloc(sizeof(qs)),
      *biggest, *smallest, *remaining;
    init_qs_ptr(valeur, sa, values, m+1, end, threshold_sort, s_threads);
    init_qs_ptr(valeur2, sa, values, start, m-1, threshold_sort, s_threads);
    
    // In which part do we have more work?
    if (end - m > m - start) {
      biggest = valeur;
      smallest = valeur2;
    } else {
      biggest = valeur2;
      smallest = valeur;
    }
      
    if (smallest->end - smallest->start < MIN_SIZE_FOR_THREAD_CALL ||
        sem_trywait(s_threads) == -1) {
      quickSort(sa, values, smallest->start, smallest->end, threshold_sort,
                s_threads);
      remaining = biggest;
      free(smallest);
    } else {
      pthread_create(&thread, NULL, quickSortT, (void*)biggest);
      // Let the thread go, we won't join it, otherwise that poses
      // some difficult challenges to wait for the threads while not
      // counting a resource (a thread) for that task (waiting).
      pthread_detach(thread);
      remaining = smallest;
    }

    if (remaining->end - remaining->start < MIN_SIZE_FOR_THREAD_CALL ||
        sem_trywait(s_threads) == -1) {
      quickSort(sa, values, remaining->start, remaining->end, threshold_sort,
                s_threads);
      free(remaining);
    } else {
      pthread_create(&thread2, NULL, quickSortT, (void*)remaining);
      pthread_detach(thread2);
    }
  
  }

  void* quickSortT(void* data){
    qs *dat=(qs *)data;
    quickSort(dat->gkSA,dat->values, dat->start, dat->end,dat->threshold_sort,dat->s_threads, true);
    // If we launched that quicksort in a new thread, we can tell that
    // it is done. Indeed, we have nothing else to do!
    sem_post(dat->s_threads);
    free(dat);
    return 0;
  }

  

  //partitionne le tableau et renvoie la position du pivot
  int64_t partition(SolArray *sa,SolArray *values,uintSA start, uintSA end,uintSA p){
    uintSA i = start;
    uintSA j = end+1;
    uintSA pivotv=values->get(p);
    uintSA pivotp=sa->get(p);
    uintSA tmp;
    swap(sa,values,p,start);
    while(true)
      {do{
	  j--;
	}
	while(((j >= start) && (pivotv < values->get(j))) || ((pivotv==values->get(j)) && (pivotp<sa->get(j))));
	do{
	  i++;
	}
	while(i <= end && values->get(i) < pivotv||((pivotv==values->get(i)) &&(pivotp>sa->get(i))));
	if(i<j)
	  {
	    tmp = values->get(i);
	    values->set(i, values->get(j));
	    values->set(j,tmp);
	    tmp = sa->get(i);
	    sa->set(i,sa->get(j));
	    sa->set(j,tmp);
	  }
	else{
	  swap(sa,values,j,start);
	  return j;
	}
      }
    return 0;
  }


  void insertion_sort(SolArray *sa , SolArray *values, int64_t &start, int64_t &end){
    int64_t i, j;
    for (i = start; i <= end; ++i) {
      uintSA elem2=sa->get(i);
      uintSA elem = values->get(i);
      for (j = i; j > 0 && (values->get(j-1) > elem||(values->get(j-1)==elem&&sa->get(j-1)>elem2)); j--){
	values->set(j,values->get(j-1));
	sa->set(j,sa->get(j-1));
      }
      values->set(j,elem);
      sa->set(j,elem2);
    }
  }

  //return the value at pos p on the strand forward
  uintSA vals(uintSA p,unsigned char dna[],int k,int64_t *last,uintSA *valfwd, uintSA *valrev){
    int e=p-*last;
    if(e!=1){
      *last=p;
      *valfwd = factorsToInt(p, dna, k);
    }else{
      // Compute the new value from the previous one.
      uintSA np=(p+k-1)/4;
      int rp=(p+k-1)%4;
      uintSA m=1;
      *valfwd%=m<<(2*k-2);
      *valfwd<<=2;
      int i=dna[np]%(1<<(8-2*rp));
      i/=(1<<(6-2*rp));
      *valfwd+=i;
      *last=p;
    }
    return *valfwd;
  }

  //return the minumum value of the k-mer at pos p between strand rev and stran fwd
  uintSA valns(uintSA p, unsigned char dna[],int k,int64_t *last,uintSA *valfwd,uintSA *valrev){
    int e=p-*last;
    if(e!=1){
      *last=p;
      *valfwd=factorsToInt(p, dna, k);
      *valrev=intRevcomp(*valfwd, k);
    }else{
      // Compute the new value from the previous one.
      uintSA np=(p+k-1)/4;
      int rp=(p+k-1)%4;
      uintSA m=1;
      *valfwd%=m<<(2*k-2);
      *valfwd<<=2;
      int i=dna[np]%(1<<(8-2*rp));
      i/=(1<<(6-2*rp));
      *valfwd+=i;
      *last=p;
      *valrev/=1<<(2);
      *valrev+=complementaryNuc(i)<<(2*k-2);
    }
    return min(*valfwd,*valrev);
  }


  //compute gkSA values on a subpart depending on a function given as parameter
  void *calc_gkSASubPart(void *data) { 
    subGkSA sgkSA=*(subGkSA *)data;
    uintSA gkSA_pos = sgkSA.offset;
    uintSA valrev,valfwd;
    int64_t last(-3);
    bool isPairedEnd = sgkSA.gk->getReads()->isPairedEnd();
    bool isStranded = sgkSA.gk->isStranded();

#ifdef HAVE_LIBPROGRESSBAR
    bool show_pb = (gkSA_pos == PB_thread_master_start_pos) && gkArrays::ShowProgressBar();
    DoccY::ProgressBar *PB = NULL;
    unsigned int maxval = 0;
    unsigned int reset_cpt = 0;
    unsigned int cpt = 0;
    if (show_pb) {
      maxval = gkArrays::GetNbColumns() ? gkArrays::GetNbColumns() : 100;
      reset_cpt = (sgkSA.end - sgkSA.start) / maxval;
      PB = new DoccY::ProgressBar ("Building SA     ", maxval, gkArrays::GetNbColumns(), cerr, false);
      PB->ShowPercent();
      PB->ShowTime();
      PB->update();
    }
#endif

    for(uintSA i=sgkSA.start;i<sgkSA.end;i++){
      if(sgkSA.gk->isPposition(i)){
        if(isStranded) {
          if(isPairedEnd && sgkSA.gk->isTheFirstMemberOfPair(sgkSA.gk->getTagNum(i))) {
            valns(i,sgkSA.cr, sgkSA.k, &last,&valfwd,&valrev);
            sgkSA.values->set(gkSA_pos,valrev);
          } else {
	    sgkSA.values->set(gkSA_pos,vals(i,sgkSA.cr, sgkSA.k, &last,&valfwd,&valrev));
          }
        } else {
	  sgkSA.values->set(gkSA_pos,valns(i,sgkSA.cr, sgkSA.k, &last,&valfwd,&valrev));
        } 
        sgkSA.gkSA->set(gkSA_pos,i);
        gkSA_pos++;
      }
#ifdef HAVE_LIBPROGRESSBAR
      if ((cpt++ == reset_cpt) && show_pb) {
	PB->Step();
	cpt = 0;
      }
#endif
    }

#ifdef HAVE_LIBPROGRESSBAR
    if (show_pb) {
      PB->SetVal(maxval);
      PB->update(false);
      delete PB;
    }
#endif
    return 0;
   } 

  void init_qs_ptr(struct qs *data, SolArray *sa, SolArray *val, int64_t start,
                   int64_t end, int threshold_sort, sem_t *s_threads) {
    data->gkSA = sa;
    data->values = val;
    data->start = start;
    data->end = end;
    data->threshold_sort = threshold_sort;
    data->s_threads = s_threads;
  }
}
