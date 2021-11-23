Clone Git Repository
====================

.. toctree::
   :maxdepth: 1

If you don't have git, please install git with the following `Install Git <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.


.. warning::

  Where you choose to install this git may be limited by "where can docker mount".

  Directories that can be mounted as Docker's mounted volumes have permission based limitations. If you don't feel confident about this section, please talk to your administrator.

Docker permission limitation:
:::::::::::::::::::::::::::::

  Docker's mounted volumes must hold the permissions: ``drwxrwxr-x``. If you're on a linux file system this means that all parent folders of your working directory must be more permissive than ``drwxrwsr-x``, except "/".

  This can be ensured with the command for each of the parent directory until '/'.

.. prompt:: bash $

  chmod a+rX,u+w,g+w <directory hosting the DICAST git>

If you have sudo, consider ``/opt/DICAST/`` as your working directory:

Cloning DICAST's git repository:
:::::::::::::::::::::::::::::::::

Clone our project repository to a directory of your choice. This directory will be considered the working directory for most of the commands listed in this documentation.

.. prompt:: bash

	git clone https://github.com/CGAT-Group/DICAST.git

This will give you access to the necessary scripts and the :doc:`directory structure <input>` hosts an example for the inputs in a directory called "sample_input". The directory structure within this git is assumed by DICAST, so please don't modify directory names within this working directory.
