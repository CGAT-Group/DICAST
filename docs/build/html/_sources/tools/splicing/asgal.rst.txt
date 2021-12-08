.. Links
.. _manual: https://asgal.algolab.eu/documentation
.. |tool| replace:: ASGAL

ASGAL
=====

.. warning::
 
 |tool| requires the variables ``$fastqpair1suffix`` and ``$fastqpair2suffix`` to be set in the :guilabel:`scripts/asevent_config.sh` file.

.. sidebar:: |tool| Factsheet
 
 ============  ==============================================================================================
 **Toolname**  *asgal*                                                                                       
 **Version**   *1.1.6*                                                                                       
 **License**   `GNU General Public License version 3 <https://github.com/AlgoLab/galig/blob/master/LICENSE>`_
 ============  ==============================================================================================
 
 **Required Files**
  * :ref:`fastq<fastqSplicing>` , :ref:`fasta<fastaSplicing>` , :ref:`gtf<gtfSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `ASGAL: aligning RNA-Seq data to a splicing graph to detect novel alternative splicing events <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2436-3>`_
 


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/asgal/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 --multi
  Set multi option.
 -g
  Reference genome in fasta format.
  
  .. code-block::
  
   -g $fasta
 -a
  The path to the gene annotation file in GTF format for annotation of fusion junctions.
  
  .. code-block::
  
   -a $gtf
 -t
  Transcript file.
  
  .. code-block::
  
   -t $transcript
 -s
  Fastq filename of paired end read 1.
  
  .. code-block::
  
   -s *yourFastqFile1_*1.fastq
 -s2
  Fastq filename of paired end read 2.
  
  .. code-block::
  
   -s2 *yourFastqFile1_*2.fastq
 -o
  Output directory. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir
 -@
  Number of threads to be used during the computation
  
  .. code-block::
  
   -@ $ncores
 --allevents
  Report all events, not only novel ones.
