.. Links
.. _manual: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
.. |tool| replace:: Subjunc

Subjunc
=======

.. sidebar:: |tool| Factsheet
 
 ============  ======================================================================================
 **Toolname**  *subjunc*                                                                             
 **Version**   *2.0.0*                                                                               
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
  * |tool| publication: `The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote  <https://www.ncbi.nlm.nih.gov/pubmed/23558742>`_
 


Indexing
^^^^^^^^

.. note::
 
 **Indexing might take some time** but only has to be run once per fasta file. Make sure to reuse already computed indices if possible.

DICAST will check if :guilabel:`$indexdir/$indexname.reads` exists. If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
If you want to use your own precomputed index file copy it to :guilabel:`index/subjunc-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.
We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/subjunc/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -i
  Base name of the index folder and files.
  
  .. code-block::
  
   -i $indexdir/$indexname
  

 -r
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -r *yourFastqFile1_*1.fastq
  

 -R
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -R *yourFastqFile1_*2.fastq
  

 -o
  The path to the directory for the **mapped** output in sam format. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir/$controlfolder/*yourFastqFile1_*subjunc
  

 -T
  Number of threads to be used during the computation
  
  .. code-block::
  
   -T $ncores
  

 --SAMoutput
  Return a sam file as output.
  
  

