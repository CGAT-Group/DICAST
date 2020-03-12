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

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <gkArrays.h>
#include <sys/time.h>

#ifndef ullong
#define ullong unsigned long long
#endif

using namespace std;
using namespace gkarrays;

ullong getChrono() {
	struct timeval time;
	time_t sec;
	suseconds_t usec;

	if (gettimeofday(&time,NULL) != -1) {
		sec = time.tv_sec;
		usec = time.tv_usec;
		return (ullong)sec*1000000+usec;
	}
	return 0;
}

int main(int argc, char **argv) {
  uint nb_queries = 100;
  time_t seed = 0;
  uint nb_threads = 1;
  bool pairedEnd = false;
  int param;
  gkArrays *tags;
  uint read_length;
  uint k;
  bool use_bitvector;
  bool is_stranded;



  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " <file> [-p <paired_end_file>] <read_length> <k> <use_bitvector> <is_stranded>  [nb_threads] [nb_queries] [seed]" << endl
         << "       <file> :          text file containing the reads in FASTA/FASTQ format. " << endl
         << "       -p <file> :       second text file containing the reads in FASTA/FASTQ format for paired-end reads. " << endl
         << "       <read_length> :   integer giving the length of the reads." << endl
         << "       <k> :             integer setting the length of indexed k-mers." << endl
         << "       <use_bitvector> : specify if we should store the data structures using" << endl
         << "                         bit vectors rather than conventionnal arrays" << endl
         << "                         (more space efficient, and more time consuming)." << endl
     	   << "       <is_stranded> :   specify if reads are oriented or not." << endl;
    exit(1);
  }

  // If -p option we have paired-end reads
  if(strcmp(argv[2],"-p") == 0) {
    pairedEnd = true;
    param = 4;
  } else {
    param = 2;
  }
  read_length = atoi(argv[param++]);
  k = atol(argv[param++]);
  use_bitvector = atol(argv[param++]);
  is_stranded = atol(argv[param++]);
  if (argc >= param+1)
    nb_threads = atol(argv[param++]);
  if (argc >= param+1) 
    nb_queries = atol(argv[param++]);
  if (argc >= param+1) 
    seed = atol(argv[param++]);

  ullong start = getChrono();
  if(pairedEnd) {
    tags = new gkArrays(argv[1],argv[3], k, use_bitvector, read_length, 
                                  is_stranded, nb_threads);
  } else {
    tags = new gkArrays(argv[1], k, use_bitvector, read_length, 
                                  is_stranded, nb_threads);
  }
  cout << "Creation: " << (getChrono()-start)/1000000. << " s" << endl;

  if (seed == 0)
    seed = time(NULL);

  srand(seed);
//   srand(seed = time(NULL));

  cout << "Pseudo-randomly choosing "<< nb_queries << " k-mers for the queries" <<endl;
  cout << "I.e., pseudo-randomly choosing "<< nb_queries << " tag numbers and tag_positions" <<endl;
  uint *tag_numbers = new uint[nb_queries];
  uint *tag_positions = new uint[nb_queries];
  uint nb_tags = tags->getNbTags();
  for (uint i=0; i < nb_queries; i++) {
    tag_numbers[i] = rand() % nb_tags;
    tag_positions[i] = rand() % (tags->getTagLength(tag_numbers[i]) - k + 1);
  }

  cout << "We have " << tags->getGkCFALength() << " distinct k-mers in the dataset." << endl
       << endl;

  cout << "Running "<< nb_queries << " queries of type Q1, Q2, Q3, Q4 with a k-mer f." <<endl;
  cout << "Query 1: what are the reads where f occurs?" << endl;

  start = getChrono();
  for (uint i=0; i < nb_queries; i++) {
    uint *matched_tags = tags->getTagNumWithFactor(tag_numbers[i],
						   tag_positions[i]);
    free(matched_tags);
  }
  cout << "Query 1: " << (getChrono()-start)*1./nb_queries << " us" << endl;

  cout << "Query 2: In how many reads does f occur?" << endl;

  start = getChrono();
  ullong total = 0;
  for (uint i=0; i < nb_queries; i++) {
    total += tags->getNbTagsWithFactor(tag_numbers[i],
				       tag_positions[i]);
  }
  cout << "Query 2: " << (getChrono()-start)*1./nb_queries << " us" << endl;
  cout << "        (total: " << total << " tags)" << endl;

  cout << "Query 3: What are the occurrence positions of f in the reads?" << endl;

  start = getChrono();
  for (uint i=0; i < nb_queries; i++) {
    pair<uint,uint> *matched_tags = tags->getTagsWithFactor(tag_numbers[i],
							    tag_positions[i]);
    delete [] matched_tags;
  }
  cout << "Query 3: " << (getChrono()-start)*1./nb_queries << " us" << endl;

  cout << "Query 4: What is the number of occurrences of f in the reads?" << endl;

  start = getChrono();
  total = 0;
  for (uint i=0; i < nb_queries; i++) {
    total += tags->getNbTagsWithFactor(tag_numbers[i],
				       tag_positions[i],
				       true);
  }
  cout << "Query 4: " << (getChrono()-start)*1./nb_queries << " us" << endl;
  cout << "        (total: " << total << " tags)" << endl;

  char *tagFactor = tags->getTagFactor(tag_numbers[0],tag_positions[0],k);
  cout << tags->getNbTagsWithFactor(tag_numbers[0], tag_positions[0]) 
       << " tags share the factor " 
       << tagFactor << " (read " << tag_numbers[0] << ", position " 
       << " " << tag_positions[0] << ")" 
       << " and occurs " << tags->getNbTagsWithFactor(tag_numbers[0],
                                                      tag_positions[0],
                                                      true)
       << " times in total" << endl;

  delete [] tagFactor;

  cout << "Seed used from pseudo-random generation: " << seed << endl;

  delete [] tag_positions;
  delete [] tag_numbers;
  delete tags;
}
