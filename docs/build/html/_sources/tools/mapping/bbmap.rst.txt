.. Links
.. _manual: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. |tool| replace:: BBMap

BBMap
=====

BBMap uses a multi kmer seed and extend strategy for read mapping.

.. sidebar:: |tool| Factsheet
 
 ============  =============================================================================================
 **Toolname**  *bbmap*                                                                                      
 **Version**   *38.43*                                                                                      
 **License**   `BBTools Copyright (c) 2014 <https://github.com/BioInfoTools/BBMap/blob/master/license.txt>`_
 ============  =============================================================================================
 
 **Required Files**
  * *For indexing only:* :ref:`fasta<fastaMapping>`
  * *For mapping:* :ref:`fastq<fastqMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/aspli`
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `Long Read RNA-seq Mapper <http://bib.irb.hr/datoteka/773708.Josip_Maric_diplomski.pdf>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/bbmap-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/bbmap/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -in
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -in *yourFastqFile1_*1.fastq
  

 -in2
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -in2 *yourFastqFile1_*2.fastq
  

 -ref
  Reference genome in fasta format.
  
  .. code-block::
  
   -ref $fasta
  

 -path
  Base name of the index folder and files.
  
  .. code-block::
  
   -path $indexdir/$indexname
  

 -intronlen
  Length of introns.
  
  .. code-block::
  
   -intronlen 20
  

 -xstag
  Add sam flags to improve compatibility with alternative splicing tools.
  
  .. code-block::
  
   -xstag us
  

 -outm
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -outm $outdir/$controlfolder/*yourFastqFile1_*bbmap.sam
  

 -outu
  The path to the **unmapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file.
  
  .. code-block::
  
   -outu $outdir/$controlfolder/*yourFastqFile1_*.bbmap_unmapped.sam
  



Known Issues
^^^^^^^^^^^^

* issue
* another
* another
