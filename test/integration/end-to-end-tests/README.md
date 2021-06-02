# Directory Contents

**`deployment-tests`**
This fold contains pytests that are run everytime we deploy to the `staging` branch. They represent `integration-tests` that are run after these other stages (found in the .gitlab-ci.yml file:

stages:
 - unit-test
 - package
 - deploy-jaws
 (- integration-test)

See `deploy-jaws-central-staging-integration-tests` in the .gitlab-ci.yml.


**`nightly-tests`**
These are tests that run "end-to-end" JAWS runs that make sure all the JAWS systems are functioning properly (differentiating them from unit-tests). These pytests are initiated by the native gitlab scheduling system (akin to cronjobs). These scheduled jobs are run by the gitlab runner with the tag of `central` and therefore, run on the machine where `central` is installed.

See `nightly-end-to-end-tests` in the .gitlab-ci.yml.


**`manually_initiated_tests`**
These are tests that are not part of the previous categories and are initiated by manually running pytest. They include things like stress tests, etc.

-------------------
## Useful Commands for the Developer

**pytest --fixture**
This command will give you some documentation about which fixtures are available to use in the test scripts.

**pytest --capture=no**
Use this if you want to have your print statements within the test scripts to be printed out to stdout, otherwise, all stdout is captured except for the pytest generated output.

#### pytest logging
In your test scripts, you can use python's logging method to control which logs you want printed to stdout.  For instance, if you are debugging tests, you can run 

`pytest -v --log-file-level debug` 

which will print out any debug, info, warning, error, or critical log levels.

You should have something like this code in your test script.

```
import logging
logging.basicConfig(level=logging.DEBUG)
mylogger = logging.getLogger()

def test_setup():
    mylogger.info('Inside Setup')

def test_setup_module():
    mylogger.debug('Inside Setup')
```

#### Log Levels
DEBUG
Detailed information, typically of interest only when diagnosing problems.

INFO
Confirmation that things are working as expected.

WARNING
An indication that something unexpected happened, or indicative of some problem in the near future (e.g. ‘disk space low’). The software is still working as expected.

ERROR
Due to a more serious problem, the software has not been able to perform some function.

CRITICAL
A serious error, indicating that the program itself may be unable to continue running.
