Config files
============

.. warning::

	Something broke while changing the config file? Make sure there is no space between the variable, the equal sign and the value.
	For example:
	
	- Wrong: workdir = "../dockers/"
	- Right: workdir="../dockers/"

To manage the input parameters for multiple tools at once you will have to modify the config files. There is one gerneral config file (config.sh) and one config file for mapping and splicing tools respectively (mapping_config.sh and asevent_config.sh) to handle application specific parameters.

.. note::
	
	Here we will only cover the most important settings. If you want do go further into detail check out our advanced section.


How to change your config.sh file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

	The following explanations assume that you use our directory structure as described in directory structure. If you want to use your own directory structure please look at the information in our advanced section.
	
All config files are in the scripts folder inside of your cloned repository.
If you are using our directory structure for your input you will have to change the following parameters:

Basic Parameters
^^^^^^^^^^^^^^^^

ncores
	| if you want to use more cores
read_length
	| check you fastq files and change to the corresponding parameter
differential
	| set to '0' for alternative splicing event detection; set to '1' for differential splicing (beta mode)

Input Parameters
^^^^^^^^^^^^^^^^
asimulator_gtf
	| the genome gtf annotation that you use to simulate the data. Default: 'Homo_sapiens.GRCh38.104.gtf'
fastaname
	| the genome reference file. Default: 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
gtfname
	| the genome gtf annotation that you use for mapping and alternative splicing analysis
gffname
	| the genome gff3 annotation that you use for mapping and alternative splicing analysis

The reference genome, annotation file and gff3 files could be downloaded from `Ensembl <http://ftp.ensembl.org/pub/release-104/>`_
