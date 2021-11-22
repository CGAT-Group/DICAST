Run one specific tool via Docker
================================

.. toctree::
   :maxdepth: 2

Change config.sh according to your run (see :doc:`How to change your config.sh file <../config/general_parameters>`)

If you have already built the image with <tool>:<tag> (see the :doc:`docker setup <../setup/docker>`) you can run the following command to run the image and start the tool:

.. prompt:: bash $

  docker run -v <your mounted folder>:/MOUNT --user $(id -u):$(id -g) <tool>:<tag>

  # Examples:
  # If you are using our directory structure for your input and are in the dockers directory:
  # Add the --rm flag to remove container, after run.
  
  docker run -v ./:/MOUNT --user $(id -u):$(id -g) --rm gsnap:0.1


Troubleshooting
^^^^^^^^^^^^^^^

  - Check Snakemake Output to see which rule failed.

  - if the rule that failed was named after a tool, check log files under output/<tool>-output/logs/ to see where the error was.
