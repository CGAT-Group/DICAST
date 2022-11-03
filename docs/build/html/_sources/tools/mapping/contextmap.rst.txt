.. Links
.. _manual: https://www.bio.ifi.lmu.de//files/Software/ContextMap/manual/ContextMap-manual.html
.. |tool| replace:: ContextMap 2.0

ContextMap 2.0
==============

.. warning::
 
 Make sure to put the file :file:`jre-8u241-linux-i586.tar.gz` in the same folder as the ContextMap Dockerfile (:guilabel:`src/contextmap2/`).\nYou can get the file here: https://www.oracle.com/java/technologies/javase/javase8u211-later-archive-downloads.html

.. sidebar:: |tool| Factsheet
 
 ============  ===========================================================================
 **Toolname**  *contextmap*                                                               
 **Version**   *2.7.9*                                                                    
 **License**   `Artistic software License <https://opensource.org/licenses/Artistic-2.0>`_
 ============  ===========================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`bowtie_fastadir<bowtie_fastadirMapping>`
 **Compatible splicing tools**
  * :doc:`../splicing/aspli`
  * :doc:`../splicing/eventpointer`
  * :doc:`../splicing/irfinder`
  * :doc:`../splicing/spladder`
  * :doc:`../splicing/whippet`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `ContextMap 2: fasta and accurate context-based RNA-seq mapping <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0557-5>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/contextmap-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/contextmap/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -reads
  Comma separated list of file paths to reads in fastq format. One pair of fastq files for paired-end mapping
  
  .. code-block::
  
   -reads *yourFastqFile1_*1.fastq,*yourFastqFile1_*2.fastq
  

 -aligner_name
  Used aligner (index tool). We use bowtie2.
  
  .. code-block::
  
   -aligner_name bowtie2
  

 -aligner_bin
  Path to the used aligner. If you use our docker you will not have to wolly about it.
  
  .. code-block::
  
   -aligner_bin /home/biodocker/bin/bowtie2
  

 -indexer_bin
  Path to the indexing tool of the aligner.
  
  .. code-block::
  
   -indexer_bin /home/biodocker/bin/bowtie2-build
  

 -indices
  Comma separated list to your index files base names.
  
  .. code-block::
  
   -indices *IndexChromosome1*,*IndexChromosomes2*,*IndexChromosome3*, . . . 
  

 -genome
  Directory path with chromosome-wise fasta files.
  
  .. code-block::
  
   -genome $bowtie_fastadir
  

