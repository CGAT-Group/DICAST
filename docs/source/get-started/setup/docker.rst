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
If you want to build the dockers manually, we provide a ``docker-compose.yml`` file which will let you build them yourself. You can use the command 

.. prompt:: bash

  docker-compose -f path/to/docker-compose.yml build

to build all images. For more information, see the  `docker compose documentation <https://docs.docker.com/compose/>`_.
@Amit

Build one image
:::::::::::::::

If you only want to build one specific docker image, run the following command:

.. prompt:: bash

  docker build <directory level path to Dockerfile> --tag=<tool>:<tag>

  # Examples:
  # If you are inside of the base directory ("dockers") and want to build gsnap:
  docker build ./gsnap/ --tag=gsnap:0.1

  # If you are already inside the dockers/gsnap directory::
  docker build ./ --tag=gsnap:0.1

or use the ``docker-compose.yml`` file and specify which image to build: 

.. prompt:: bash

  docker-compose -f path/to/docker-compose.yml build <tool>

3. Other helpful commands
^^^^^^^^^^^^^^^^^^^^^^^^^
List docker images (for example to get image ids):

.. prompt:: bash $

  docker images

Remove an image:

.. prompt:: bash $

  docker rmi -f <image id>
