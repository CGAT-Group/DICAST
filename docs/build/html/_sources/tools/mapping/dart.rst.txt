.. Links
.. _manual: https://github.com/hsinnan75/Dart
.. |tool| replace:: Dart

Dart
====

.. sidebar:: |tool| Factsheet
 
 ============  ===============================================================================================
 **Toolname**  *dart*                                                                                         
 **Version**   *1.4.6*                                                                                        
 **License**   `GNU General Public License version 2 <https://github.com/hsinnan75/Dart/blob/master/LICENSE>`_
 ============  ===============================================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`fasta<fastaMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/aspli`
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/majiq`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `DART: a fast and accurate RNA-seq mapper with a partitioning strategy <https://academic.oup.com/bioinformatics/article/34/2/190/4104410>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname.sa` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/dart-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/dart/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -i
  Base name of the index folder and files.
  
  .. code-block::
  
   -i $indexdir/$indexname
  

 -f
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -f *yourFastqFile1_*1.fastq
  

 -f2
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -f2 *yourFastqFile1_*2.fastq
  

 -o
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*dart.sam
  

