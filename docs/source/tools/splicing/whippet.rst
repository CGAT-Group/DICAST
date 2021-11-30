.. Links
.. _manual: https://github.com/timbitz/Whippet.jl
.. |tool| replace:: Whippet

Whippet
=======

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  ==========================================================================
 **Toolname**  *whippet*                                                                 
 **Version**   *0.11.1*                                                                  
 **License**   `MIT Licence <https://github.com/timbitz/Whippet.jl/blob/master/LICENSE>`_
 ============  ==========================================================================
 
 **Required Files**
  * :ref:`fastq<fastqSplicing>` , :ref:`fasta<fastaSplicing>` , :ref:`gtf<gtfSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `Efficient and Accurate Quantitative Profiling of Alternative Splicing Patterns of Any Complexity on a Laptop <https://pubmed.ncbi.nlm.nih.gov/30220560/>`_
 


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/whippet/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --fasta
  Reference genome in fasta format.
  
  .. code-block::
  
   --fasta $fasta
 --gtf
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   --gtf $gtf
 --bam
  Name of bamfile.
  
  .. code-block::
  
   --bam $controlfolder/*filename*.bam
 -x
  Output directory for whippet index.
  
  .. code-block::
  
   -x $outdir/*bamfilename*/graph
 -o
  Output directory. The output will be separated into case and control folder based on the basefolder of the according bam file. If you are running the DICAST pipeline to compare different mapping tools this will include the name of the mapping tool of the used bam file.
  
  .. code-block::
  
   -o $outdir
