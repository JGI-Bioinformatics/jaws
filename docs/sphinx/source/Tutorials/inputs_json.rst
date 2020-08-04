===========================
How to Create a Inputs File
===========================

An inputs json file has the basic format:

.. code-block:: bash 

   {
    "<workflow name>.<task name>.<variable name>": "<variable type>"
   }


.. note::
	Jaws always expects an inputs json file even if it contains only open and close brackets {}.


You can create an inputs file by scratch or you can build a template based on the WDL using the following command:

.. code-block:: bash 

   jaws run inputs <path to your.wdl>

This command should output a template for the inputs.json file. You can then fill in the vaules of each key.

.. code-block:: bash 

   {
     "bbtools.reads": "File",
     "bbtools.ref": "File"
   }

.. note::
	File paths may be absolute or relative paths, but must contain a slash in the path (e.g. `/tmp/x` or `./x` but not `x`).

To see how to include lists or dictionaries in your input.json files, go to the official wdl docs site for `inputs.json <https://software.broadinstitute.org/wdl/documentation/inputs>`_

