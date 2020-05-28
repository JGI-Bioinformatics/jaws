===========================
How to Create a Inputs File
===========================

An inputs json file has the basic format:

.. code-block:: bash 

   {
    "<workflow name>.<task name>.<variable name>": "<variable type>"
   }


.. note::
	Jaws always expects an inputs json file even if its only open and close brackets.


You can create an inputs file by scratch or you can let the following command build you a template:

.. code-block:: bash 

   jaws wdl inputs <path to your.wdl>

This command should output something like this:

.. code-block:: bash 

   {
     "bbtools.reads": "File",
     "bbtools.ref": "File"
   }

.. note::
	File paths may be absolute or relative paths, but must begin with either a `/` or a `.` (e.g. `/tmp/x` or `./x` but not `x`).

To see more, go to the official wdl docs site for `inputs.json <https://software.broadinstitute.org/wdl/documentation/inputs>`_
There are some examples in this repo `jaws-example-wdl <https://gitlab.com/jfroula/jaws-example-wdl>`_, but you'll need to search through the different example folders to find `inputs.json` files. 

*************
Usage
*************

An example of how you would use the inputs file would be like this:

.. code-block:: bash

   jaws run submit inputs.json my.wdl

