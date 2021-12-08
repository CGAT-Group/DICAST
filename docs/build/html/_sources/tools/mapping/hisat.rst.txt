.. Links
.. _manual: http://daehwankimlab.github.io/hisat2/manual/
.. |tool| replace:: HISAT2

HISAT2
======

.. sidebar:: |tool| Factsheet
 
 ============  ===================================================================================
 **Toolname**  *hisat*                                                                            
 **Version**   *2.2.1*                                                                            
 **License**   `GNU General Public License version 3 <https://www.gnu.org/licenses/gpl-3.0.html>`_
 ============  ===================================================================================
 
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
  * |tool| publication: `Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype <https://www.nature.com/articles/s41587-019-0201-4>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/${indexname}.4.ht2 and $indexdir/${gtfname}_splicesites.txt` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/hisat-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/hisat/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --x
  Base name of the index folder and files.
  
  .. code-block::
  
   --x $indexdir/$indexname
  

 -1
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -1 *yourFastqFile1_*1.fastq
  

 -2
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -2 *yourFastqFile1_*2.fastq
  

 -S
  The path to the **mapped** output file in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -S $outdir/$controlfolder/*yourFastqFile1_*hisat.sam
  

 --known-splicesite-infile
  Provide a list of known splice sites.
  
  .. code-block::
  
   --known-splicesite-infile $indexdir/$indexname/splicesites.txt
  

 -q
  Activate quiet mode so only error messages are printed.
  
  

