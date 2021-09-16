.. Links
.. _manual: http://crac.gforge.inria.fr/
.. |tool| replace:: CRAC

CRAC
====

.. sidebar:: |tool| Factsheet
 
 ============  ======================================================================================
 **Toolname**  *crac*                                                                                
 **Version**   *2.4.0*                                                                               
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
  * |tool| publication: `CRAC: an integrated approach to the analysis of RNA-seq reads <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-3-r30>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname.ssa` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/crac-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/crac/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -i
  Base name of the index folder and files.
  
  .. code-block::
  
   -i $indexdir/$indexname
  

 -k
  Number of k-mers to be used. 22 is the recommended number for human genome.
  
  .. code-block::
  
   -k 22
  

 -r
  Space separated list of file paths to reads in fastq format. One pair of fastq files for paired-end mapping
  
  .. code-block::
  
   -r *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
  

 -o
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*crac.sam
  

 --detailed-sam
  Return a detailed sam file as output.
  
  

 --stranded
  Reads are from a strand specific RNA-seq protocol.
  
  

