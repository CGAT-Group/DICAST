.. Links

.. _manual: *not available*
.. |tool| replace:: Asgal

Asgal
=====

.. warning::

	Asgal requires the variables ``$fastqpair1suffix`` and ``$fastqpair2suffix`` to be set in the :guilabel:`scripts/asevent_config.sh` file.

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *asgal*
	**Version:**   *v1.1.1*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastqdir/*$fastqpair1suffix
		$fastqdir/*$fastqpair2suffix
		$fasta
		$transcript
		$gtf

|tool|


Default parameters:
^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/asevent_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	--multi


		.. code-block:: bash

			--multi

	-g
		Fasta reference file.

		.. code-block:: bash

			-g $fasta

	-a
		Gtf annotation file.

		.. code-block:: bash

			-a $gtf

	-t
		Gene transcripts in fasta format

		.. code-block:: bash

			-t $transcript

	-s
		Fastq filename of paired end read 1.

		.. code-block:: bash

			-f *yourFastqFile1_*1.fastq

	-s2
		Fastq filename of paired end read 2.

		.. code-block:: bash

			-f2 *yourFastqFile1_*2.fastq

	-o
		The path to the output directory for the according fastq file pair. The file will be named after the fastq file basename.

		.. code-block:: bash

			-o $outdir/*yourFastqFile1*-output

	-@
		Set number of threads to be used during the computation.

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			#    $ncores=4

			-@ $ncores

3. Other comments:
^^^^^^^^^^^^^^^^^^


4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication:
