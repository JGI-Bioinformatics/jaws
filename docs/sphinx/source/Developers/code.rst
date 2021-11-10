.. _code.rst:

===============
Developer Guide
===============

This guide is intended for people who wish to contribute to the JAWS
core itself.


--------------
Overview
--------------
JAWS is designed to run workflows by leveraging Workflow Domain Language(WDL) specifications
and Cromwell. It is meant to help users run workflows across many different systems.

--------------
Code Structure
--------------

This section is designed to give users an overview of the different code bases.
JAWS consists of different services that all exist in the same code base. They consist
of the JAWS central codebase, JAWS site code and client code. 

-----------
JAWS Client
-----------

JAWS client contains the client code needed to interact with the Central API service. It
implements a command line interface

^^^^
CLI
^^^^
.. automodule:: client.jaws_client.cli
   :members:

^^^^^
Workflow helpers
^^^^^

.. automodule:: client.jaws_client.workflow
   :members:

^^^^^
Copy progress
^^^^^
.. automodule:: client.jaws_client.copy_progress
   :members:


^^^^
WDL runtime validator
^^^^
:mod:`wdl_runtime_validator.py` contains all the logic for validating the
runtime section of WDLs. 

.. automodule:: client.jaws_client.wdl_runtime_validator
   :members:

------------
JAWS Central
------------

^^^^^^^^
API
^^^^^^^^
The `jaws_central.analysis` file contains the API for the central application. It uses RPC style calls
for processing requests.

.. automodule:: central.jaws_central.analysis
   :members:


^^^^^^
Globus
^^^^^^

:class:`jaws_central.globus.GlobusService` contains the class used for transferring globus transfers using the 
Globus SDK.

.. automodule:: central.jaws_central.globus
   :members:

^^^^^^ 
Auth
^^^^^^

.. automodule:: central.jaws_central.auth
   :members:
