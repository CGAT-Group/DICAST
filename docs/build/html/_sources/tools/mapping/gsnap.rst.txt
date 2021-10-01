.. Links
.. _manual: http://research-pub.gene.com/gmap/
.. |tool| replace:: GSNAP

GSNAP
=====

.. sidebar:: |tool| Factsheet
 
 ============  ==================================================================
 **Toolname**  *gsnap*                                                           
 **Version**   *2020-03-12*                                                      
 **License**   `Apache Licence 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_
 ============  ==================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`fasta<fastaMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `Fast and SNP-tolerant detection of complex variants and splicing in short reads <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2844994/>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname/$indexname.contig` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/gsnap-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/gsnap/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --db
  Base name of the index folder and files.
  
  .. code-block::
  
   --db $indexdir/$indexname
  

 --dir
  Base folder of the index files.
  
  .. code-block::
  
   --dir $indexdir
  

 --output-file
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   --output-file $outdir/$controlfolder/*yourFastqFile1_*gsnap.sam
  

 --format
  Define output format (one of sam,  m8).
  
  .. code-block::
  
   --format sam
  

 --force-xs-dir
  Add sam flags to improve compatibility with alternative splicing tools.
  
  .. code-block::
  
   --force-xs-dir us
  

 --nthreads
  Number of threads to be used during the computation
  
  .. code-block::
  
   --nthreads $ncores
  

 reads
  After all other options call space separated list of file paths to reads in fastq format. One pair of fastq files for paired-end reads.
  
  .. code-block::
  
   *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
  

