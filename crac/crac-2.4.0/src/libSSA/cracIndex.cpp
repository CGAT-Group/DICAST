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

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <list>
#include <cstring>
#include <cstdlib>
#include <ctime>
using namespace std;


extern "C" {
#include <sys/types.h>
#include <unistd.h>
}

#include <blockwise_sa.h>
#include <SSA.h>
#include <config.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <cracIndex.h>
#include <GenomeInfo.h>
#include <utils.h>

typedef String<Dna, Packed<> > TStr;

using namespace seqan;

int verbose = 0;

void usage(char *prog) {
  cerr << endl << prog << " version " << PACKAGE_VERSION;
  cerr << "\tCompiled on " <<__DATE__ << endl << endl;
  cerr << "Usage : "<<prog << " [options] <command> <output file> <input file>+"
       << endl << endl 
       << "  command must be one of:" << endl
       << "    index: create an index on the specified input file(s)." 
       << endl << endl
       << "    get:   get a (multi)FASTA file containing the original sequences."
       << endl << endl
       << "  options can be (for the index command only):" << endl << endl
       << "  -b <bucket_size>\t the size of the bucket for the index construction"
       << endl 
       << "                  \t (default " << THIS_BUCKET_SIZE << ")" << endl
       << "  -d <diff-cover> \t parameter for the index construction (default "
       << COVER << ")"
       << endl
       << "  -s <sample dist>\t sample distance (default "
       << DEFAULT_SAMPLE_DIST << ")" << endl
       << "  -v              \t verbose mode"
       << endl << endl 
       << "  Examples: " << endl
       << "\tIndexing:" << endl
       << "\t\t" << prog << " index myIndex sequence1.fa sequence2.fa sequence3.fa"
       << endl
       << "\t\t\tYou can specify FASTA or MultiFASTA file(s)." << endl
       << "\t\t\tIn this example, two files will be created:" << endl
       << "\t\t\t- myIndex.ssa (index storing the compressed sequences)" << endl
       << "\t\t\t- myIndex.conf (information on sequence names and length)" << endl
       << "\tExtracting:" << endl
       << "\t\t" << prog << " get sequences.fa myIndex" << endl
       << "\t\t\tSequences indexed in myIndex will be extracted" << endl
       << "\t\t\tto the sequences.fa file" << endl;
  exit(1);
}


void index(char **input_files, int nb_input_files, char *output_file, 
           uint cover, uint bucket_size, uint sample_distance) {
  unsigned long long n = 0;
  void *index;
  ifstream file;
  ofstream conf;
  TStr *text = new TStr();
  uchar *bwt;
  int i = 0;
  list<pair<uint, String<char> > > sequences;

  for (i = 0; i < nb_input_files; i++) {
    if (verbose) 
      cerr << "Reading file " << input_files[i] << endl;

    file.open(input_files[i], ios_base::in | ios_base::binary);
    if (file.is_open()) {
      String<Dna5> current_seq;
      String<char> fasta_tag;
      uint current_length;
      do {
        readMeta(file, fasta_tag, Fasta());
        read(file, current_seq, Fasta());
        current_length = length(current_seq);
        if (verbose && current_length)
          cerr << "\tSequence of length " << current_length << endl;
        if (current_length) {
          // Search sequence name in fasta_tag
          Iterator<String<char> >::Type it = begin(fasta_tag);
          while (it != end(fasta_tag) && value(it) != ' ')
            it++;
          fasta_tag = prefix(fasta_tag, it);

          if (length(fasta_tag) > MAX_LENGTH_CHR_NAME) {
            cerr << "Warning the name of the sequence is too long and will be cut:" << endl
                 << "\t" << fasta_tag << endl;
            fasta_tag = prefix(fasta_tag, MAX_LENGTH_CHR_NAME-1);
          }
          
          sequences.push_back(pair<uint, String<char> >(current_length,
                                                        fasta_tag));
          
          uint this_seq_length = length(current_seq);
          for (uint i = 0; i < this_seq_length; i++) {
            // Replace N with a random nucleotide... that's not great. I know.
            // But the text is encoded on two bits.
            if (current_seq[i] == 'N') {
              current_seq[i] = randomDNA();
            }
          }

          n += current_length;
          *text += current_seq;
        }
      } while (current_length > 0);
      file.close();
    } else {
      cerr<<"Unable to open the input file " << input_files[i] <<  endl;
      exit(2);
    }
  }


  // Creating configuration file
  char *confName = new char[strlen(output_file)+6];
  sprintf(confName, "%s.conf", output_file);
  conf.open(confName, ios_base::out | ios_base::binary);
  if (! conf.is_open()) {
    cerr << "Can't write configuration file " << confName << endl;
    exit(3);
  }
  delete [] confName;
  
  // Writing configuration file
  list<pair<uint, String<char> > >::iterator it = sequences.begin();
  conf << sequences.size();
  while (it != sequences.end()) {
    conf << endl << (*it).second 
         << endl << (*it).first;
    it++;
  }
  conf.close();
  
  if (n > ((unsigned long long) 1 << 32) - 1) {
    cerr << "Your genome sequence is too large. I cannot index it. Sorry" << endl;
    exit(4);
  }

  if (verbose) 
    cerr << "Initialising KBSA" <<  endl;

  BlockwiseSA<TStr> *kSA = new KarkkainenBlockwiseSA<TStr>(*text, bucket_size, cover);
  
  if (verbose)
    cerr << "Building BWT. Do you know what the BWT is?" << endl
         << "You can have a look at http://enwp.org/Burrows-Wheeler%20transform" << endl;

  // Build the BWT of chars.
  bwt = new unsigned char[n+1];
  for (uint i = 0; i <= n; i++) {
    uint current = kSA->nextSuffix();
    if (current == 0)
      bwt[i] = 255;
    else {
      bwt[i] = convert<char>((*text)[current - 1]);
    }
  }
  delete kSA;
  delete text;

  if (verbose)
    cerr << "Building the BWT-based index." << endl
         << "Do you know what an FM-index is?" << endl
         << "You can have a look at http://enwp.org/FM-index" << endl;

  // build_index has been modified so that it takes a BWT as parameter
  char *params = new char[30];
  sprintf(params, "terminator %c;samplerate %u",255, sample_distance);
  build_index((unsigned char *)bwt, n+1, params, &index);
  if (verbose)
    cerr << "Saving the index in the output file" << endl;
  save_index(index, output_file);
  delete [] params;
  free_index(index);
}


void get(char *input, char *output) {
  void *index;
  unsigned int n;
  uchar *chars;
  uint total_pos;
  uint current_pos;
  ofstream fasta;
  GenomeInfo *gInfo;

  if (verbose)  
    cerr << "Reading informations on sequences" << endl;
  gInfo = new GenomeInfo(input);

  fasta.open(output, ios_base::out | ios_base::binary);
  if (! fasta.is_open()) {
    cerr << "Can't write FASTA file " << output << endl;
    exit(3);
  }

  if (verbose)
    cerr << "Loading index" << endl;
  load_index(input,&index);
  get_length(index, &n);
  if (verbose)
    cerr << "Extracting the sequences from the index" << endl;
  extract(index, 0, n-1, &chars, &n);

  uint i = 0;
  uchar substr[NB_CHAR_FASTA_LINE+1];
  substr[NB_CHAR_FASTA_LINE]=0;
  total_pos = 0;
  while (i < gInfo->getNbChr()) {
    current_pos = 0;
    // Output FASTA header for the sequence
    if (verbose)
      cerr << "Retrieving " << gInfo->getChrName(i) << endl;
    fasta << ">" << gInfo->getChrName(i) << endl;
    // Reading each sequence
    while (current_pos + NB_CHAR_FASTA_LINE < gInfo->getChrLength(i)) {
      // Writing each line separately
      strncpy((char *)substr, (char *)&chars[total_pos], NB_CHAR_FASTA_LINE);
      fasta << substr << endl;
      total_pos += NB_CHAR_FASTA_LINE;
      current_pos += NB_CHAR_FASTA_LINE;
    }
    // Read the remaining (if any)
    if (current_pos < gInfo->getChrLength(i)) {
      strncpy((char *) substr, (char *) & chars[total_pos], gInfo->getChrLength(i) - current_pos);
      substr[gInfo->getChrLength(i) - current_pos] = 0;
      fasta << substr << endl;
      total_pos += gInfo->getChrLength(i) - current_pos;
    }
    i++;
  }
  fasta.close();
  free(chars);

  free_index(index);

  delete gInfo;
}


int main(int argc, char **argv) {
  uint bucket_size = THIS_BUCKET_SIZE;
  uint cover = COVER;
  uint sample_distance = DEFAULT_SAMPLE_DIST;
  char *output_file = NULL;
  int command=-1;                  // O: index, 1: get
  int option_index;
  int opt;

  static const char *short_options = "d:b:s:hv";
  static struct option long_options[] = 
    {{"help", no_argument, NULL, 'h'},
     {"diff-cover", required_argument, NULL, 'd'},
     {"bucket-size", required_argument, NULL, 'b'},
     {"sample-dist", required_argument, NULL, 's'},
     { 0, 0, 0, 0 } //terminator
    };

  while ((opt = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1){
    switch(opt) {
    case 'd' : 
      cover = atoi(optarg);
      break;
    case 'b' : 
      bucket_size = atoi(optarg);
      break;
    case 's' :
      sample_distance = atoi(optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h' : 
    default: 
      usage(argv[0]);
    }
  }

  if (argc - optind < 3) {
    usage(argv[0]);
  }

  // Get the command argument
  if (! strcmp(argv[optind], "index")) {
    command = 0;
  } else if (! strcmp(argv[optind], "get")) {
    command = 1;
  } else {
    cerr << "Unknown command " << argv[optind] << endl;
    usage(argv[0]);
  }
  optind++;

  // Get the output file
  output_file = argv[optind++];

  // Input files are what is remaining...

  if (command == 0) {
    // index
    index(&argv[optind], argc - optind, output_file, cover, bucket_size,
          sample_distance);
  } else {
    // get
    get(argv[optind], output_file);
  }

    exit(0);
}
