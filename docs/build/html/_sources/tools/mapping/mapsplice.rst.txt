.. Links
.. _manual: http://www.netlab.uky.edu/p/bioinfo/MapSplice2UserGuide
.. |tool| replace:: MapSplice 2

MapSplice 2
===========

.. sidebar:: |tool| Factsheet
 
 ============  ======================================================================================================
 **Toolname**  *mapsplice*                                                                                           
 **Version**   *latest conda version*                                                                                
 **License**   `MapSplice Copyright (C) 2012-2013 <https://github.com/merckey/MapSplice2/blob/master/Copyright.txt>`_
 ============  ======================================================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`bowtie_fastadir<bowtie_fastadirMapping>`, :ref:`gtf<gtfMapping>`
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
  * |tool| publication: `MapSplice: Accurate mapping of RNA-seq reads for splice junction discovery <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952873/>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname.rev.2.ebwt` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/mapsplice-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/mapsplice/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -c
  Directory path with chromosome-wise fasta files.
  
  .. code-block::
  
   -c $bowtie_fastadir
  

 -x
  Base name of the index folder and files.
  
  .. code-block::
  
   -x $indexdir/$indexname
  

 --gene-gtf
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --gene-gtf $gtf
  

 -o
  The path to the directory for the **mapped** output in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*mapsplice
  

 -p
  Number of threads to be used during the computation
  
  .. code-block::
  
   -p $ncores
  

 -1
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -1 *yourFastqFile1_*1.fastq
  

 -2
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -2 *yourFastqFile1_*2.fastq
  

