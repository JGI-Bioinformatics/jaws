========================
Current JGI Subworkflows
========================

.. role:: bash(code)
   :language: bash

A subworkflow is just a WDL that can be imported into a main WDL and then called like a task. 

Here is a more detailed description of subworkflows on `Cromwell's docs <https://Cromwell.readthedocs.io/en/stable/SubWorkflows>`_

####################
Why use subworkflows
####################

* re-use code
* use code that has been expertly vetted

#######################
Example of Subworkflow
#######################

**Main WDL: main.wdl**

.. figure:: /Figures/main.png

**Subworkflow: sub.wdl**

.. figure:: /Figures/sub.png

################################
How to Run Subworkflows in JAWS
################################

You can try running this example

.. code-block:: bash

    # activate the environment you set up above
    source ~/jaws-prod.sh

    # clone the example code 
    git clone https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories/jaws-tutorial-examples.git

    cd jaws-tutorial-examples/subworkflow

    # run jaws run submit <workflow> <inputs> <site>
    jaws run submit main.wdl inputs.json cori
