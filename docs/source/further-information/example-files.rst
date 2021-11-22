Examples
========

Here we provide the detailed description of possible workflows. We recommend to run analysis using a terminal multiplexer, e.g. tmux or screen.

Running multiple mapping tools (E.g., STAR, HISAT2 and bbmap)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Make sure you followed the steps described in the :doc:`setup <../get-started/setup>` section carefully.

2. Before getting started make sure to activate the snakemake conda environment:

.. prompt:: bash

  conda activate dicast-snakemake

3. Create the `input` folder:

.. prompt:: bash

   cd /path/to/DICAST/
   mkdir input

4. Create the directory structure as in the `sample_output`:

.. prompt:: bash

   cd input
   mkdir controldir
   cd controldir
   mkdir fastqdir

5. Download or copy the genome fasta file into the `input` folder. Dont't forget to uncompress it. E.g.:

.. prompt:: bash

   cd /path/to/DICAST/input
   wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

6. Download or copy the genome gtf annotation into the `input` folder. Dont't forget to uncompress it. E.g.:

.. prompt:: bash

   wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
   gunzip Homo_sapiens.GRCh38.104.gtf.gz

7. Download or copy the fastq files you want to align into the `/path/to/DICAST/input/controldir/fastqdir`. Note: we support only paired-end RNA-Seq - fastq files have to be in pairs.

8. Go to `/path/to/DICAST/scripts` and edit `config.sh` according to your run (see :doc:`How to change your config.sh file <../get-started//config/general_parameters>`):

.. prompt:: bash

   cd /path/to/DICAST/scripts
   nano config.sh

In the `config.sh` file edit the following lines:

.. prompt::

   read_length=76
   fastaname=Homo_sapiens.GRCh38.dna.primary_assembly.fa
   gtfname=Homo_sapiens.GRCh38.104.gtf

9. List the mapping tools you want to run:

.. prompt:: bash

   cd /path/to/DICAST/scripts/snakemake/
   nano snakemake_config.yaml

In the `snakemake_config.yaml` file edit the following lines:

.. prompt::

   Mapping_tools:
       What_tools_to_run: 'star, hisat, bbmap'

10. In the `/path/to/DICAST/scripts/snakemake/` folder run:

.. prompt:: bash

   snakemake -j 1 -d /path/to/DICAST/input -s Snakefile -c snakemake_config.yaml

This command will start the mapping tools indicated in the `snakemake_config.yaml` (E.g. STAR, HISAT2 and bbmap).

First, the pipeline will build all necessary dockers. Second, in will create a `/path/to/DICAST/index` folder and put the results of indexing. Finally, the pipeline will create a `/path/to/DICAST/output` folder with the alignment results inside the dedicated folders (e.g., star-output, hisat-output, bbmap-output).

Running multiple alternative splicing event detection tools (E.g., MAJIQ and Whippet)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Make sure you followed the steps described in the :doc:`setup <../get-started/setup>` section carefully.

2. Before getting started make sure to activate the snakemake conda environment:

.. prompt:: bash

  conda activate dicast-snakemake

3. Create the `input` folder:

.. prompt:: bash

   cd /path/to/DICAST/
   mkdir input

4. Create the directory structure as in the `sample_output`:

.. prompt:: bash

   cd input
   mkdir controldir
   cd controldir
   mkdir fastqdir
   mkdir bamdir

5. Download or copy the genome fasta file into the `input` folder. Dont't forget to uncompress it. E.g.:

.. prompt:: bash

   cd /path/to/DICAST/input
   wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

6. Download or copy the genome annotation file into the `input` folder. Dont't forget to uncompress it. E.g.:

.. prompt:: bash

   wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
   gunzip Homo_sapiens.GRCh38.104.gtf.gz

7. Download or copy the genome gff3 annotation into the `input` folder (for MAJIQ). Dont't forget to uncompress it. E.g.:

.. prompt:: bash

   wget http://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz
   gunzip Homo_sapiens.GRCh38.104.gff3.gz

8. Download or copy the fastq files you want to use into the `/path/to/DICAST/input/controldir/fastqdir`. Note: we support only paired-end RNA-Seq - fastq files have to be in pairs.

9. Download or copy the bam files you want to use into the `/path/to/DICAST/input/controldir/bamdir`.

10. Go to `/path/to/DICAST/scripts` and edit `config.sh` according to your run (see :doc:`How to change your config.sh file <../get-started/config/general_parameters>`):

.. prompt:: bash

   cd /path/to/DICAST/scripts
   nano config.sh

In the `config.sh` file edit the following lines:

.. prompt:: bash

   read_length=76
   fastaname=Homo_sapiens.GRCh38.dna.primary_assembly.fa
   gtfname=Homo_sapiens.GRCh38.104.gtf
   gffname=Homo_sapiens.GRCh38.104.gff3

11. List the mapping tools you want to run:

.. prompt:: bash

   cd /path/to/DICAST/scripts/snakemake/
   nano snakemake_config.yaml

In the `snakemake_config.yaml` file edit the following lines:

.. prompt::

   Alternative_splicing_detection_tools:
      What_tools_to_run: 'majiq, whippet'

12. In the `/path/to/DICAST/scripts/snakemake/` folder run:

.. prompt:: bash

   snakemake -j 1 -d /path/to/DICAST/input -s Snakefile -c snakemake_config.yaml

This command will start the mapping tools indicated in the `snakemake_config.yaml` (E.g. MAJIQ, Whippet).

First, the pipeline will build all necessary dockers. Second, the pipeline will create a `/path/to/DICAST/output` folder with the event detecton results inside the dedicated folders (e.g., majiq-output, hisat-output, whippet-output).

.. toctree::
   :maxdepth: 2
