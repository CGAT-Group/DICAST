.. Links
.. _manual: https://bioconductor.org/packages/release/bioc/vignettes/ASpli/inst/doc/ASpli.pdf
.. |tool| replace:: Aspli

Aspli
=====

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  ==============================================================================================
 **Toolname**  *aspli*                                                                                       
 **Version**   *1.12.0*                                                                                      
 **License**   `GNU General Public License version 3 <https://github.com/AlgoLab/galig/blob/master/LICENSE>`_
 ============  ==============================================================================================
 
 **Required Files**
  * :ref:`gtf<gtfSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `ASpli: integrative analysis of splicing landscapes through RNA-Seq assays  <https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab141/6156815>`_
 

.. note::
 
 |tool| is an R package. Therefore our ENTRYPOINT.sh script for |tool| calls an R script to run the tool. The parameters listed here are the parameters given to the R script.


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/aspli/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --gtf
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --gtf $gtf
 --cores
  Number of threads to be used during the computation
  
  .. code-block::
  
   --cores $ncores
 --readLength
  Length of reads.
  
  .. code-block::
  
   --readLength $read_length
 --out
  Output directory. The output will be separated into case and control folder based on the basefolder of the according bam file. If you are running the DICAST pipeline to compare different mapping tools this will include the name of the mapping tool of the used bam file.
  
  .. code-block::
  
   --out $outdir
 --differential
  1 to run differential analysis, 0 otherwise.
  
  .. code-block::
  
   --differential $differential
