================================
Best Practices for Creating WDLs
================================

This link has helpful examples

`WDL best practices (scott frazer's gist) <https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211>`_


Docker Images
-----------------------------------
* save you docker images at hub.docker.com. Jaws will pull images from there by default.

* include a "Dockerfile" in your git repository

* keep older version of the docker images so older workflows can be reproduced

* use conda or venv environments instead of docker for testing

* use Cromwell.jar instead of jaws for testing


Workflows/Tasks:
----------------
* Workflows should be comprised of distinct tasks than can be executed outside the context of the workflow

* tasks should take file paths as inputs, not folders

* separate tasks which use varying compute resources whenever possible

* grouping tasks into sub-workflows makes them more understandable, reuseable, and more easily maintained

* if multiple commands are used in a "command" stanza, they should be chained using `&&` or `set -euo pipefail` so subsequent commands are not executed after the first failure


There are opportunities to participate in code reviews with other WDL developers; there is also a #jaws-developers slack channel (see below link).

The `ContactUs <contact_us.html>`_ page offers ways for getting help.

