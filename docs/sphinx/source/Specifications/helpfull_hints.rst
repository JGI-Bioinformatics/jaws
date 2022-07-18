###################################
Helpfull hints when developing WDLs
###################################

Here is a link to the developer's `wiki <https://bitbucket.org/berkeleylab/jgi-workflows/wiki/Home>`_  page.  It has valuable gotha's and advice and general helpful hints when developing a WDL.

Please feel free to add problems that you have solved that will be helpful to others. It is editable by anyone with write permissions to the "jgi-workflows" repository.  If you would like permissions to edit the wiki please email jaws-support@lbl.gov or you can just send an email with suggestions.


# Using JAWS from your own service

If you have your own service which needs to interact with JAWS, there are a few options.

## Method 1: Use the jaws-client provided by the jaws team.

If you want to run the client as a command-line tool, you can simply add the client to your $PATH.  This would be appropriate for services not written in python.

.. code-block:: text

    export PATH=/global/cfs/projectdirs/jaws/jaws-prod:$PATH
    export JAWS_CLIENT_CONFIG=/global/cfs/projectdirs/jaws/jaws-prod/jaws-prod.conf
    export JAWS_USER_CONFIG=~/jaws.conf

This requires that your service has access to the same file system (e.g. on a cori node, not a VM).

## Method 2: Install the jaws-client into your own python virtual environment

.. code-block:: text

    source /your/venv/bin/activate
    pip install wheel
    git clone https://code.jgi.doe.gov/advanced-analysis/jaws.git
    cd ./jaws/client
    python python setup.py bdist_wheel
    pip install pip install ./dist/jaws_client*whl

You must also set the environment variables (including `$PATH`) as shown in Method 1.

## Method 3: Use the REST interface

The REST interface at `jaws.lbl.gov` may be used to query runs but does not support submissions at this time.
