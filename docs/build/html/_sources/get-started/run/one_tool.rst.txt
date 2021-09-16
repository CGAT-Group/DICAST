Run one specific tool
=====================

.. toctree::
   :maxdepth: 2

Change config.sh according to your run (see :doc:`How to change your config.sh file <../run/config>`)

If you have already built the image with <tool>:<tag> (see the :doc:`docker setup <../setup/docker>`) you can run the following command to run the image and start the tool:

.. prompt:: bash $

  docker run -v <your mounted folder>:/MOUNT --user $(id -u):$(id -g) <tool>:<tag>

  # Examples:
  # If you are using our directory structure for your input and are in the dockers directory:
  docker run -v ./:/MOUNT --user $(id -u):$(id -g) gsnap:0.1

  # Add the --rm flag to clean up after yourself
  docker run -v ./:/MOUNT --user $(id -u):$(id -g) --rm gsnap:0.1
