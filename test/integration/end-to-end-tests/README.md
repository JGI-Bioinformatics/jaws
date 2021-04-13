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

