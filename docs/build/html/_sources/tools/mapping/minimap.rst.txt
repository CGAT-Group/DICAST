.. Links
.. _manual: https://lh3.github.io/minimap2/minimap2.html
.. |tool| replace:: Minimap2

Minimap2
========

.. sidebar:: |tool| Factsheet
 
 ============  ========================================================================
 **Toolname**  *minimap*                                                               
 **Version**   *latest conda version*                                                  
 **License**   `MIT Licence <https://github.com/lh3/minimap2/blob/master/LICENSE.txt>`_
 ============  ========================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`fasta<fastaMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/majiq`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `Minimap2: pairwise alignment for nucleotide sequences <https://doi.org/10.1093/bioinformatics/bty191>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/minimap-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/minimap/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -a
  Generate CIGAR and provide output in sam format.
  
  

 -o
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*minimap.sam
  

 -t
  Number of threads to be used during the computation
  
  .. code-block::
  
   -t $ncores
  

 index
  Base name of the index folder and files.
  
  .. code-block::
  
   $indexdir/$indexname
  

 reads
  After all other options call space separated list of file paths to reads in fastq format. One pair of fastq files for paired-end reads.
  
  .. code-block::
  
   *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
  

