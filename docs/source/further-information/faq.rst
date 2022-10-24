FAQ
===

Here you will find more frequently asked questions soon.

Q: How do I contribute to DICAST?
    | A: The best way to reach us for code updates is via our `github <https://github.com/CGAT-Group/DICAST/issues>`_



.. toctree::
   :maxdepth: 2
    further-information/uninstall-dicast

Q: How do I resolve issue: `docker: Error response from daemon: error while creating mount source path..`
    | A: Check if the folder that docker is trying to mount has the following permissions:`drwxrwsr-x`. Grant them when needed with `chmod a+rX,u+w,g+w `.