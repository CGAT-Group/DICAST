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
#include <gkCFPSBuilder.h>
#include "basic.h"
#include "utils.h"
#include "solArray.h"
#include <iostream>
#include "gkArrays.h"
#ifdef HAVE_LIBPROGRESSBAR
#include <libProgressBar/progressBar.h>
#endif
using namespace std;

namespace gkarrays {

  GkCFPSBuilder::GkCFPSBuilder():gkCFPS(NULL) {}
  GkCFPSBuilder::~GkCFPSBuilder() {}

  void GkCFPSBuilder::build(SolArray *values, uintSA length, array_type t, 
                            uint nb_threads){

    // Building GkCFPS
    gkCFPS = values;    

    // It is so quick that we don't need to bother with multithread construction.
    monothreadConstruction(values, length);
  }
  
  SolArray *GkCFPSBuilder::getGkCFPS(){
    return gkCFPS;
  }

  void GkCFPSBuilder::monothreadConstruction(SolArray *values, uintSA length) {
    uintSA bucket_id = ~0;
    uintSA nb_buckets = 0;
    uintSA pos_gkCFPS = 0;

#ifdef HAVE_LIBPROGRESSBAR
    DoccY::ProgressBar *PB = NULL;
    unsigned int maxval = 0;
    unsigned int reset_cpt = 0;
    unsigned int cpt = 0;
    if (gkArrays::ShowProgressBar()) {
      maxval = gkArrays::GetNbColumns() ? gkArrays::GetNbColumns() : 100;
      reset_cpt = length / maxval;
      PB = new DoccY::ProgressBar ("Building gkCFPS ", maxval, gkArrays::GetNbColumns(), cerr, false);
      PB->ShowPercent();
      PB->ShowTime();
      PB->update();
    }
#endif

    for (uintSA i = 0; i < length; i++) {
      if (values->get(i) != bucket_id) {
        bucket_id = values->get(i);
        gkCFPS->set(pos_gkCFPS, nb_buckets);
        pos_gkCFPS++;
      } 
      nb_buckets++;
#ifdef HAVE_LIBPROGRESSBAR
      if ((cpt++ == reset_cpt) && gkArrays::ShowProgressBar()) {
	PB->Step();
        cpt = 0;
      }
#endif
    }
    gkCFPS->set(pos_gkCFPS, nb_buckets);
    gkCFPS->realloc(pos_gkCFPS+1);
#ifdef HAVE_LIBPROGRESSBAR
    if (gkArrays::ShowProgressBar()) {
      PB->SetVal(maxval);
      PB->update(false);
      cerr << endl;
      delete PB;
    }
#endif
  }
}
