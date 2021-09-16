.. Links
.. _manual: http://kissplice.prabi.fr/
.. |tool| replace:: kisSplice

kisSplice
=========

.. note::
 
 |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
 If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet
 
 ============  ===========================================================================================================
 **Toolname**  *kissplice*                                                                                                
 **Version**   *2.5.1*                                                                                                    
 **License**   `GNU General Public License version >=2 <https://bioconductor.org/packages/release/bioc/html/kissDE.html>`_
 ============  ===========================================================================================================
 
 **Required Files**
  * :ref:`fastq<fastqSplicing>` , :ref:`gtf<gtfSplicing>` , :ref:`bam<bamSplicing>`
 **Links**
  * |tool| `manual`_
  * |tool| publication: `KIS SPLICE: de-novo calling alternative splicing events from RNA-seq data <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-S6-S5>`_
 


Parameters
^^^^^^^^^^

These are the default parameters set in the :guilabel:`src/kissplice/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.

 -r
  List of fastq files seperated by “ -r “
  
  .. code-block::
  
   -r *yourFastqFile1_*1.fastq -r *yourFastqFile1_*2.fastq
 -o
  Output directory. The output will be separated into case and control folder based on the basefolder of the according fastq file. 
  
  .. code-block::
  
   -o $outdir
