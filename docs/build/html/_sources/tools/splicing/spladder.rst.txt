.. Links
.. _manual: https://spladder.readthedocs.io/en/latest/
.. |tool| replace:: SplAdder

SplAdder
========

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  ============================================================================
 **Toolname**  *spladder*                                                                  
 **Version**   *2.4.2*                                                                     
 **License**   `BSD License <https://github.com/ratschlab/spladder/blob/master/COPYRIGHT>`_
 ============  ============================================================================
 
 **Required Files**
  * :ref:`gtf<gtfSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `SplAdder: identification, quantification and testing of alternative splicing events from RNA-Seq data <https://pubmed.ncbi.nlm.nih.gov/26873928/>`_
 


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/spladder/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -b
  Name of bamfile.
  
  .. code-block::
  
   -b $controlfolder/*filename*.bam
 -o
  Output directory. The output will be separated into case and control folder based on the basefolder of the according bam file. If you are running the DICAST pipeline to compare different mapping tools this will include the name of the mapping tool of the used bam file.
  
  .. code-block::
  
   -o $outdir
 -a
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   -a $gtf
 --parallel
  Number of threads to be used during the computation
  
  .. code-block::
  
   --parallel $ncores
 -n
  Length of reads.
  
  .. code-block::
  
   -n $read_length
 --output-txt-conf
  Output in txt format.
