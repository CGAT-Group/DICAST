Setup docker
============

.. toctree::
  :numbered:
  :maxdepth: 1

.. .. sectnum:: depth:1

1. Get docker-compose  for your system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This requires you to have administrative rights on your computer, or talk to your system administrator about getting Docker on your system and giving you rights to use the docker user group.

A. Download Docker
:::::::::::::::::::::

Follow the Docker Engine installation manual for getting Docker, your first step. : https://docs.docker.com/engine/install/

B. Post-install Docker steps
::::::::::::::::::::::::::::

.. _post-install: https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user

To run DICAST in user mode entirely, please fulfill this `post-install`_ step to run docker as a non-root user.

C. Install docker-compose
::::::::::::::::::::::::::

This while closely related, docker-compose is the last of DICAST's docker dependencies. Follow the installation manual from docker-compose.
https://docs.docker.com/compose/install/#install-compose


Basic Docker commands
::::::::::::::::::::::
Should you have docker configured on your system, you shouldn't run into a permission error for the following commands.

.. prompt:: bash

  docker images

to list docker images in on your computer.

.. prompt:: bash

  docker ps

to list all running containers only.

.. prompt:: bash

  docker ps -a

to list all running and stopped containers.

.. prompt:: bash

  docker --version


We support Docker version 19 and above.

.. prompt:: bash

  docker-compose --version



2. Build docker images (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While the steps described in this section are handled by DICAST's graphical interface, it can also be accessed via command line, for more control.

Build all images
::::::::::::::::::

If you intend to use multiple dockers at once you can use our snakemake pipeline, which will take care of building the docker images.
If you want to build the dockers manually, we provide a ``docker-compose.yml`` file which will let you build them yourself. You can use the command the following command to build all images.

.. prompt:: bash

  docker-compose -f scripts/Snakemake/docker-compose.yml build

If you'd like to edit DICAST's docker-compose file, see the  `docker-compose Manual <https://docs.docker.com/compose/gettingstarted/>`_.

Build one image
:::::::::::::::

If you only want to build one specific docker image, run the following command to first build some core essential containers:

.. prompt:: bash

  docker-compose -f scripts/Snakemake/docker-compose.yml build base conda bowtie star

And if you want to build any of the other tools, use the following command:

.. prompt:: bash

  docker-compose -f scripts/Snakemake/docker-compose.yml build <tool>

Where <tool> needs to be replaced with one or more of the following tools:

bbmap, contextmap, crac, dart, gsnap, hisat, mapsplice, minimap, segemehl, star, subjunc
asgal, aspli, eventpointer, irfinder, majiq, sgseq, spladder, whippet


3. Pull docker images (Fail safe)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Should the tool you intended to run, not build locally, it's also possible to pull them from DICAST's dockerhub repository at: https://hub.docker.com/repository/docker/dicastproj/dicast


.. prompt:: bash $

  docker pull dicastproj/dicast:tagname

4. Other helpful commands
^^^^^^^^^^^^^^^^^^^^^^^^^
To gracefully stop a running docker container (If perhaps snakemake's process had to be killed):

.. prompt:: bash $

  docker stop <docker-container-name/ID>

Remove an image (to save space, after your analysis):

.. prompt:: bash $

  docker rmi -f <image id>
