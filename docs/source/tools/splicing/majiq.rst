.. Links
.. _manual: https://majiq.biociphers.org/
.. |tool| replace:: MAJIQ

MAJIQ
=====

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  =======================================================================
 **Toolname**  *majiq*                                                                
 **Version**   *2.1-05c9c32*                                                          
 **License**   `Academic License <https://majiq.biociphers.org/Academic_License.txt>`_
 ============  =======================================================================
 
 **Required Files**
  * :ref:`gff<gffSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `A new view of transcriptome complexity and regulation through the lens of local splicing variations <https://elifesciences.org/articles/11752>`_
 


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/majiq/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 build reference
  The path to the gene annotation file in GFF format.
  
  .. code-block::
  
   $gff
 -c
  MAJIQ config file (built based on DICAST config parameters in ENTRYPOINT.sh)
  
  .. code-block::
  
   -c $config
 -j
  Number of threads to be used during the computation
  
  .. code-block::
  
   -j $ncores
 -o
  Output directory with majiq build output. 
  
  .. code-block::
  
   -o $outdir/$outdir_name/build
 maqiq psi
  Run MAJIQ in psi mode with files built from gff as input.
 -j
  Number of threads to be used during the computation
  
  .. code-block::
  
   -j $ncores
 -o
  Output directory with psi output. Used to build splicegraph with voila.
  
  .. code-block::
  
   -o $outdir/$outdir_name/psi
 -n
  Run with bam files as input.
  
  .. code-block::
  
   -n “BAM”
