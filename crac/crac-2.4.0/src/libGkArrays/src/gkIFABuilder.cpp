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

#include <cstdlib>
#include <pthread.h>
#include <gkIFABuilder.h>
#include "basic.h"
#include "utils.h"
#include "solArray.h"
#ifdef HAVE_LIBPROGRESSBAR
#include <libProgressBar/progressBar.h>
#endif

using namespace std;

namespace gkarrays {

  GkIFABuilder::GkIFABuilder():gkIFA(NULL) {}
  GkIFABuilder::~GkIFABuilder() {}

  void GkIFABuilder::build(gkArrays *gk, SolArray *gkSA, SolArray *gkCFPS, 
                           array_type t, uint nb_threads){
    uintSA nb_kmers = gkSA->length();
    struct ifa_builder_thread ibt[nb_threads];
    pthread_t threads[nb_threads];
    // Building GkIFA
    gkIFA = new SolArray(nb_kmers, t);

    // Monothread and multithread are equivalent
    for (uint i = 0; i < nb_threads; i++) {
      ibt[i].gk = gk;
      ibt[i].gkSA = gkSA;
      ibt[i].gkCFPS = gkCFPS;
      ibt[i].gkIFA = gkIFA;
      ibt[i].thread_id = i;
      ibt[i].nb_threads = nb_threads;

      pthread_create(&threads[i], NULL, thread_work_ifa, (void *)&ibt[i]);
    }

    for (uint i = 0; i < nb_threads; i++) {
      pthread_join(threads[i], NULL);
    }
  }
  
  SolArray *GkIFABuilder::getGkIFA(){
    return gkIFA;
  }

  void *thread_work_ifa(void *data) {
    struct ifa_builder_thread *ibt = (struct ifa_builder_thread *)data;
    gkArrays *gk = ibt->gk;
    SolArray *gkSA = ibt->gkSA;
    SolArray *gkCFPS = ibt->gkCFPS;
    SolArray *gkIFA = ibt->gkIFA;
    uintSA len_gkCFPS = gkCFPS->length();
    uintSA pos_gkCFPS = 1 + ibt->thread_id * (len_gkCFPS / ibt->nb_threads);
    uintSA end = pos_gkCFPS + len_gkCFPS / ibt->nb_threads;
    uintSA i = gkCFPS->get(pos_gkCFPS - 1);

    if (ibt->thread_id == ibt->nb_threads - 1)
      end = len_gkCFPS;

#ifdef HAVE_LIBPROGRESSBAR
    DoccY::ProgressBar *PB = NULL;
    unsigned int maxval = 0;
    unsigned int reset_cpt = 0;
    unsigned int cpt = 0;
    bool show_pb = !ibt->thread_id && gkArrays::ShowProgressBar();
    if (show_pb) {
      maxval = gkArrays::GetNbColumns() ? gkArrays::GetNbColumns() : 100;
      reset_cpt = (end - pos_gkCFPS) / maxval;
      PB = new DoccY::ProgressBar ("Building IFA    ", maxval, gkArrays::GetNbColumns(), cerr, false);
      PB->ShowPercent();
      PB->ShowTime();
      PB->update();
    }
#endif
    while (pos_gkCFPS < end) {
      uintSA nb_in_bucket = gkCFPS->get(pos_gkCFPS) - gkCFPS->get(pos_gkCFPS-1);
      for (uintSA j = 0; j < nb_in_bucket; j++) {
        gkIFA->set(gk->convertPposToQpos(gkSA->get(i+j)), pos_gkCFPS - 1);
      }
      pos_gkCFPS++;
      i += nb_in_bucket;
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
    return NULL;
  }
}
