===========================
How to Create a Inputs File
===========================

An inputs json file has the below format.

Typically, you pass a variable to the workflow and within the WDL, 
you can pass the variable to one or more tasks.

.. code-block:: text 

   {
    "<workflow name>.<variable name>": "<value>"
   }

Alternatively, you can pass a variable to an individual task

.. code-block:: text 

   {
    "<workflow name>.<task name>.<variable name>": "<value>"
   }


.. note::
    Jaws always expects an inputs json file even if it contains only open and close brackets {}.


You can create an inputs file by scratch or you can build a template based on the WDL using the following command:

.. code-block:: text 

   jaws inputs <path to your.wdl>

This command should output a template for the inputs.json file. You can then fill in the values of each key.

.. code-block:: text 

   {
     "bbtools.reads": "File",
     "bbtools.ref": "File"
   }

.. warning::
    File paths may be absolute or relative paths, but must start with `./`, `../`, or `/` to be recognized as a path. Globus requires this syntax to ensure the file is copied to the working directory, however, we will changes this so JAWS is in alignment with the practices accepted by openWDL.

To include lists or dictionaries in your input.json files, you would do something like:

.. code-block:: text

   {
      "bbtools.var_list": ["a","b","c"]
      "bbtools.var_map": { "first":  "a",
                           "second": "b",
                           "third":  "c" }
   }


