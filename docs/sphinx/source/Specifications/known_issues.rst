######################
Known Issues with JAWS
######################

.. role:: bash(code)
  :language: bash

* JAWS v2.1. There is sometimes a serious lag (e.g. 30min+) for the commands :bash:`jaws run task-log` and :bash:`jaws run task-status`. These are issues that will be dealt with with a re-designing of how they get their information and should be fixed in JAWS v2.2.

|

* Sometimes JAWS runs will sit in the "queued" state for a long time.  If one of the various daemons that are watching for new runs goes down, we need to do a manual re-start.  For now, we have a test jaws run occur every hour to catch down services and need to restart daemons manually.  Automated re-starts may also be part of JAWS 2.2.
