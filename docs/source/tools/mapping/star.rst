.. Links
.. _manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
.. |tool| replace:: STAR

STAR
====

.. sidebar:: |tool| Factsheet
 
 ============  ======================================================================
 **Toolname**  *star*                                                                
 **Version**   *2.7.5*                                                               
 **License**   `MIT Licence <https://github.com/alexdobin/STAR/blob/master/LICENSE>`_
 ============  ======================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`fasta<fastaMapping>`, :ref:`gtf<gtfMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/aspli`
  * :doc:`../splicing/eventpointer`
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/majiq`
  * :doc:`../splicing/sgseq`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `STAR: ultrafast universal RNA-seq aligner <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$star_index/$indexname/genomeParameters.txt` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/star-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/star/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --sjdbGTFfile
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --sjdbGTFfile $gtf
  

 --readFilesIn
  Space separated list of file paths to reads in fastq format. One pair of fastq files for paired-end mapping
  
  .. code-block::
  
   --readFilesIn *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
  

 --genomeDir
  Base name of the index folder and files.
  
  .. code-block::
  
   --genomeDir $indexdir/$indexname
  

 --outFileNamePrefix
  The path to the directory for the **mapped** output in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   --outFileNamePrefix $outdir/$controlfolder/*yourFastqFile1_*star
  

 --runTreadN
  Number of threads to be used during the computation
  
  .. code-block::
  
   --runTreadN $ncores
  

 --twopassMode
  Basic 2-pass mapping,  with all 1st pass junctions inserted into the genome indices on the fly
  
  .. code-block::
  
   --twopassMode Basic
  

 --outSAMstrandField
  Add strand derived from the intron motif.
  
  .. code-block::
  
   --outSAMstrandField intronMotif
  

 --outSAMattributes 
  Add sam flags to improve compatibility with alternative splicing tools.
  
  .. code-block::
  
   --outSAMattributes  us
  

