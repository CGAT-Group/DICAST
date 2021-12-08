.. Links
.. _manual: https://github.com/williamritchie/IRFinder/wiki
.. |tool| replace:: IRFinder

IRFinder
========

.. sidebar:: |tool| Factsheet
 
 ============  ===============================================================================
 **Toolname**  *irfinder*                                                                     
 **Version**   *1.3.1*                                                                        
 **License**   `MIT Licence <https://github.com/williamritchie/IRFinder/blob/master/LICENSE>`_
 ============  ===============================================================================
 
 **Required Files**
  * :ref:`fastq<fastqSplicing>` , :ref:`fasta<fastaSplicing>` , :ref:`gtf<gtfSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `IRFinder: assessing the impact of intron retention on mammalian gene expression <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4>`_
 

.. note::
 
 |tool| can use both fastq and bam files. To use bamfiles please set the parameter $use_bam_input_files=1, and =0 to use fastq files in the as_config.sh script.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/irfinder/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -r
  Base folder of the index files.
  
  .. code-block::
  
   -r $indexdir
 -d
  Output directory. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -d $outdir
 reads
  After all other options call space separated list of file paths to reads in fastq format. One pair of fastq files for paired-end reads.
  
  .. code-block::
  
   *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
