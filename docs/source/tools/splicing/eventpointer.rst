.. Links
.. _manual: https://www.bioconductor.org/packages/release/bioc/vignettes/EventPointer/inst/doc/EventPointer.html
.. |tool| replace:: EventPointer

EventPointer
============

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  =============================================================================
 **Toolname**  *eventpointer*                                                               
 **Version**   *2.4.0*                                                                      
 **License**   `Artistic software License 2 <https://opensource.org/licenses/Artistic-2.0>`_
 ============  =============================================================================
 
 **Required Files**
  * :ref:`gtf<gtfSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `EventPointer: an effective identification of alternative splicing events using junction arrays <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2816-x>`_
 

.. note::
 
 |tool| is an R package. Therefore our ENTRYPOINT.sh script for |tool| calls an R script to run the tool. The parameters listed here are the parameters given to the R script.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/eventpointer/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --gtf
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --gtf $gtf
 --cores
  Number of threads to be used during the computation
  
  .. code-block::
  
   --cores $ncores
 --out
  Output directory. The output will be separated into case and control folder based on the basefolder of the according bam file. If you are running the DICAST pipeline to compare different mapping tools this will include the name of the mapping tool of the used bam file.
  
  .. code-block::
  
   --out $outdir
 --bamfolder
  Location of bam files.
  
  .. code-block::
  
   --bamfolder $controlfolder
 --differential
  1 to run differential analysis, 0 otherwise.
  
  .. code-block::
  
   --differential $differential
