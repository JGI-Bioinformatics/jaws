######################
Known Issues with JAWS
######################

.. role:: red
.. role:: bash(code)
  :language: bash

.. raw:: html

	<span class="red">Fixed in v2.3</span>

* JAWS v2.1. There is sometimes a serious lag (e.g. 30min+) for the commands :bash:`jaws task-log` and :bash:`jaws run task-status`. These are issues that will be dealt with with a re-designing of how they get their information and should be fixed in JAWS v2.2.

.. raw:: html

	<span class="red">Fixed in v2.3</span>

* Sometimes JAWS runs will sit in the "queued" state for a long time.  If one of the various daemons that are watching for new runs goes down, we need to do a manual re-start.  For now, we have a test run every hour to catch down services and need to restart daemons manually.  Automated re-starts may also be part of JAWS 2.2.
Fixed in v2.3
