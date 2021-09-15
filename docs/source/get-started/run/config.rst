Config files
============

.. warning::

	Something broke while changing the config file? Make sure there is no space between the variable, the equal sign and the value.
	For example:
	
	- Wrong: workdir = "../dockers/"
	- Right: workdir="../dockers/"

To manage the input parameters for multiple tools at once you will have to modify the config files. There is one gerneral config file (config.sh) and one config file for mapping and splicing tools respectively (mapping_config.sh and asevent_config.sh) to handle application specific parameters.

.. note::
	
	Here we will only cover the most important settings. If you want do go further into detail check out our advanced section.


How to change your config.sh file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

	The following explanations assume that you use our directory structure as described in folder structure. I you want to use your own folder structure please look at the information in our advanced section.
	
All config files are in the scripts folder inside of your cloned repository.
If you are using our folder structure for your input you will have to change the following parameters:

-
