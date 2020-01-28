=====================================
setup your JAWS runtime environment:
=====================================

.. role:: bash(code)
   :language: bash

.. code-block:: bash

   # Do this once to get access to the "jaws" command 
   ln -s /global/project/projectdirs/jaws/prod/cli/ ~/.conda/envs/jaws

   # Do this each time you open a new terminal
   ssh cori
   module load python/2.7-anaconda-2019.07
   source activate jaws

Depending on how you set up your conda, the command may be "conda activate" instead of "source activate".


Test by running :bash:`jaws`

*******************
Testing environment
*******************
To create a testing environment when developing WDLs, follow the instructions in :doc:`The basics </Tutorials/wdl_development>`
