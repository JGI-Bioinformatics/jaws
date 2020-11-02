===========================
How to Create a Inputs File
===========================

An inputs json file has the below format.

Typically, you pass a variable to the workflow and within the WDL, 
you can pass the variable to one or more tasks.

.. code-block:: bash 

   {
    "<workflow name>.<variable name>": "<value>"
   }

Alternatively, you can pass a variable just to a individual task

.. code-block:: bash 

   {
    "<workflow name>.<task name>.<variable name>": "<value>"
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

.. warning::
    File paths may be absolute or relative paths, but must contain a slash in the path (e.g. `/tmp/x` or `./x` but not `x`). Furthermore, any relative paths that start with ../x will cause an error.

To include lists or dictionaries in your input.json files, you would do something like:

.. code-block:: bash

   {
      "bbtools.var_list": ["a","b","c"]
      "bbtools.var_map": { "first":  "a",
                           "second": "b",
                           "third":  "c" }
   }


