.. Links
.. _manual: https://bioconductor.org/packages/release/bioc/vignettes/SGSeq/inst/doc/SGSeq.html
.. |tool| replace:: SGSeq

SGSeq
=====

.. sidebar:: |tool| Factsheet
 
 ============  =============================================================================
 **Toolname**  *sgseq*                                                                      
 **Version**   *1.22.0*                                                                     
 **License**   `Artistic software License 2 <https://opensource.org/licenses/Artistic-2.0>`_
 ============  =============================================================================
 
 **Required Files**
  * :ref:`gtf<gtfSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `Prediction and Quantification of Splice Events from RNA-Seq Data <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156132>`_
 

.. note::
 
 |tool| is an R package. Therefore our ENTRYPOINT.sh script for |tool| calls an R script to run the tool. The parameters listed here are the parameters given to the R script.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/sgseq/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --gtf
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --gtf $gtf
 --path_to_bam
  Name of bamfile.
  
  .. code-block::
  
   --path_to_bam $controlfolder/*filename*.bam
 --out
  Output directory. The output will be separated into case and control folder based on the basefolder of the according bam file. If you are running the DICAST pipeline to compare different mapping tools this will include the name of the mapping tool of the used bam file.
  
  .. code-block::
  
   --out $outdir
 --cores
  Number of threads to be used during the computation
  
  .. code-block::
  
   --cores $ncores
