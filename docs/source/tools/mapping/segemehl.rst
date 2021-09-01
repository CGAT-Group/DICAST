.. Links
.. _manual: https://www.bioinf.uni-leipzig.de/Software/segemehl/
.. |tool| replace:: segemehl

segemehl
========

.. sidebar:: |tool| Factsheet
 
 ============  ======================================================================================
 **Toolname**  *segemehl*                                                                            
 **Version**   *0.3.4*                                                                               
 **License**   `GNU General Public License version 3 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_
 ============  ======================================================================================
 
 **Required Files**
  * *For mapping:* :ref:`fastq<fastqMapping>`, :ref:`fasta<fastaMapping>`
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
  * |tool| publication: `Fast Mapping of Short Sequences with Mismatches, Insertions and Deletions Using Index Structures <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000502>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/segemehl-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/segemehl/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -d
  Reference genome in fasta format.
  
  .. code-block::
  
   -d $fasta
  

 -q
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -q *yourFastqFile1_*1.fastq
  

 -q
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -q *yourFastqFile1_*2.fastq
  

 -i
  Base name of the index folder and files.
  
  .. code-block::
  
   -i $indexdir/$indexname
  

 --splits
  Use split reads alignment
  
  

 -o
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*segemehl.sam
  

 -t
  Number of threads to be used during the computation
  
  .. code-block::
  
   -t $ncores
  

