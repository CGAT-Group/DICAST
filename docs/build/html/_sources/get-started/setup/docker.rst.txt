Setup docker
============

.. toctree::
  :numbered:
  :maxdepth: 1

.. .. sectnum:: depth:1

1. Get docker for your system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Linux
:::::

https://docs.docker.com/engine/install/ubuntu/

MacOS
:::::

https://docs.docker.com/docker-for-mac/install/

Windows
:::::::

https://docs.docker.com/docker-for-windows/install/

2. Build docker images
^^^^^^^^^^^^^^^^^^^^^^

Build all images
::::::::::::::::::

If you intend to use multiple dockers at once you can use our snakemake pipeline, which will take care of building the docker images.

@Amit

Build one image
:::::::::::::::

If you only want to build one specific docker image run the following command:

.. prompt:: bash

  docker build <directory level path to Dockerfile> --tag=<tool>:<tag>

  # Examples:
  # If you are inside of the base directory ("dockers") and want to build gsnap:
  docker build ./gsnap/ --tag=gsnap:0.1

  # If you are already inside the dockers/gsnap directory::
  docker build ./ --tag=gsnap:0.1

3. Run docker images
^^^^^^^^^^^^^^^^^^^^

Use our run script
::::::::::::::::::

If you intend to use multiple dockers at once you can use our run script, which will take care of running the docker images.
@Amit

Run one image
:::::::::::::

If you have already built the image with <tool>:<tag> you can run the following command to run the image and start the tool:

.. prompt:: bash $

  docker run -v <your mounted folder>:/MOUNT --user $(id -u):$(id -g) <tool>:<tag>

  # Examples:
  # If you are using our directory structure for your input and are in the dockers directory:
  docker run -v ./:/MOUNT --user $(id -u):$(id -g) gsnap:0.1

  # Add the --rm flag to clean up after yourself
  docker run -v ./:/MOUNT --user $(id -u):$(id -g) --rm gsnap:0.1


4. Other helpful commands
^^^^^^^^^^^^^^^^^^^^^^^^^
List docker images (for example to get image ids):

.. prompt:: bash $

  docker images

Remove an image:

.. prompt:: bash $

  docker rmi -f <image id>
