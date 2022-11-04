DICAST Outputs
==================

DICAST provides outputs as you would expect them from each Alternative Splicing tool within :guilabel:`output/<astoolname>-output/<Fastq-filename>_output`.

DICAST also provides a :guilabel:`output/<astoolname>-output/<Fastq-filename>_output_dicast_unified` output format for each tool. This is a simple tsv file that hosts all the events found from each tool. We used this to unify the outputs needed to build each of the plots outputted by DICAST.


.. prompt:: bash $

  output/
  ├── <astoolname>-output
  │   ├── logs
  │   ├── <Fastq-filename>_output
  │   ├── <Fastq-filename>_output_<astoolname>_dicast_unified
  └── plots
      └── <Fastq-filename>
          ├── <mapping tool>-name
          │   ├── A3_compare.png
          │   ├── A5_compare.png
          │   ├── AFE_compare.png
          │   ├── ALE_compare.png
          │   ├── ES_compare.png
          │   ├── IR_compare.png
          │   ├── MEE_compare.png
          │   ├── MES_compare.png
          │   └── overall_compare.png
          └── unmapped
              ├── A3_compare.png
              ├── A5_compare.png
              ├── AFE_compare.png
              ├── ALE_compare.png
              ├── ES_compare.png
              ├── IR_compare.png
              ├── MEE_compare.png
              ├── MES_compare.png
              └── overall_compare.png

------------------------------

DICAST also outputs an UpSet plot for each ``Fastq-filename``-``mapping_tool`` combination.

.. figure:: ../../img/upset_plot.png

This plot shows the events that were found in common by tools and shows you which tools found these events as well.

------------------------------

When run with **ASimulatoR**, DICAST also outputs precision and recall plots for each ``Fastq-filename``-``mapping_tool`` combination.;

for all events

.. figure:: ../../img/overall_compare.png

and for each event

.. figure:: ../../img/ES_compare.png
