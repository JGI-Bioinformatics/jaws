# Changelog

# Inbound changes not assigned to version # yet
- Add function to chmod files/folders transferred via Globus (#1739)
- Increase wall clock time to 14 days for dori site (#1754))
- change NMDC site's refdata dir path (#1756)
- add function to chmod files/folders transferred via Globus (#1739)
- try to address db update failures after to long transfer operation (#1746)
- Cromwell run submission error message shall be saved in the run log (#1568)

### External Facing Changes

### Internal Facing Changes

## 2.0.2  Summary [10/23/2023 - staging // 11/09/2023 - Prod]
### External Facing Changes
- added JOB_ID (HTCondor_ID) to tasks report (#1742)
- increase allowed length of wdl and json filenames (#1744)
- tasks report shall report "runtime_minutes" (#1745 jaws-support#141)
- fix workflow_root not defined bug (#1747)

## 2.0.1  Summary [09/19/2023 - staging // 10/09/2023 - Prod]
- update the access times for copied input files to ensure the files are not purged and so call caching works (#1735)
- Add image_version.yml file to keep track of container deployments (#1733)
- cancelling runs shall set "result" field too (#1724)
- fixed bug where if WDL name started with compute-site name, unexpected path substitution occurred (#1713)
- Change dori site's max ram to 1500
- Add "cpu_hours" to run status info and task log after run completes (#1719)
- task-summary has been deprecated as it's been merged with the task-log (#1719); command was renamed from "task-log" to simply "tasks"
- cached tasks are now added to the task-log after a run completes (#1727)
- cancelling a run now updates the task-log (#1664)
- task-log format has changed (!1520)
- task-summary shall use timestamps from the task log instead of Cromwell metadata (#1611)
- ensure infiles exist upon (re)submit before sending to Cromwell; i.e. haven't been purged (jaws-support#110)
- update atime for infiles when (re)submit; to avoid purging files prematurely (#1689)
- outputs.json shall contain relpaths instead of abspaths (#1652)
- runs' output folders shall once again contain a copy of the WDL, inputs-JSON, and subworkflows-ZIP (#1710)
- fixed bug when Cromwell submission fails during input processing, before a folder is created, and was not recognized by JAWS (#1711)
- corrected outputs folder permissions (#1712)
- fixed transfer error when file was named pipe instead of regular file (#1725)

## 2.0.0 Summary [08/23/2023]
### Internal Facing Changes
- use Cromwell hogGroups to prevent a single user from hogging entire cluster
- Use `parallel_sync` Python package for file transfers
- Add close() in update_status to force use of a fresh connection on the next update
- chmod's the output files using concurrent.futures
- get the Pagurus filewriter tests working again, also disentangle filewriter class from script.
- bugfix cancel where compute resources were not being released
- bugfix in the scrontab.sh script that prevented deployment to perlmutter

## 1.1.1 Summary [05/10/2023]
- Update Dockerfile base image to 3.11 and replace setup.py with pyproject.toml (#1635)
- dori now is now a viable site because globus configuration has been fixed

## 1.1.0 Summary [04/12/2023]
- improved support for generating metadata for large runs (>10k tasks)
- added 'jaws resubmit <run_id>' command
- Enabled setfacl for dori
- You can now use nested-array datastructures in your outputs{} section in the wdl (#1620)

## 1.1.0 Summary [03/24/2023 --> Staging]
- ci/cd now makes the end-to-end test manual and allowed to fail or not run
- task-log is now much faster and supports >10k tasks; it no longer uses cromwell metadata and so there are no longer records for cached tasks which did not actually execute
- fixed performance_metrics path to match htcondor's path
- Added job in gitlab-ci.yml to deploy to dori in addition to cori (triggered by MR is created)

## 1.0.13b Summary [03/01/2023]
- AMS CI/CD 
- A merge request will trigger the end-to-end pytest in jaws-tests 

## 1.0.13 Summary [02/24/2023 --> Staging]
- Update max_cpu check

## 1.0.12 Summary [02/08/2023 --> Staging][02/23/2023 --> PROD]
- accurately set task start time for cancelled tasks too
- add max_cpu check, similar to max_ram check so user doesn't exceed available resources on site. This only works for hard-coded values in the runtime.

### Internal
- Changed cromwell_run_id to run_id in dispatch table of rpc_operations.py for cancel_run

## 1.0.11 Summary [01/27/2023]
### User Facing

- throttle user runs to only run 10 cromwell submissions at a time
- fixed bug that included queue time with run time

## 1.0.10 Summary
### User Facing

- add new jaws-nmdc site
- site names are now lower-cased everywhere for consistency
- add support for Pair and Map outputs


### Internal
- Remove main.py from ticket #1380
- pagurus now part of site repo
- remove build.sh jobs from ci/cd pipeline

## 1.0.5 (2022-9-26) Summary
- Copy fixes from jaws_central related to environment variables to control configuration.

## 3.0 (2022-7-14)Summary
### User Facing Changes
#### Changed
- Condor and AWS complete beta-testing
- sort task-log and task-summary by "queued" time
- "done" is now always the final state of a run.
- errors report shall include .log files if they exist
- new "running-tasks" command will show logs of active tasks
- fixed shifterimg pull via tag bug

## 2.8.6 (2022-6-21) Summary
### User Facing Changes
#### Changed

- MR#1254 Add child class of ConfigParser that allows properly named environment variables to override config file settings.
- Restore /refdata (NB: not implemented for AWS)
- Created user group to restrict access to private compute site.
- Fixed bug where stderr file contents were not included in errors report.
- Fixed bug where AWS upload failure did not cause Run to change state to "upload failed" (stayed in "uploading")
- Fixed bug where very long WDLs could not be saved in Db and would fail.
- Cancelled runs will now be downloaded and available via the "get" command (useful for debugging).
- Include the stdout.submit in the errors report (previously only included stderr.submit file)
- Include transfer failure error message in run log
- Always send email (previously didn't after transfer failed or cancelled)

## 2.8.5 (2022-5-24) Summary
### User Facing Changes
#### Changed
- Install performance metrics script pagurus in site venv during CICD deployment
- add beta deployment at AWS
- remove "--default-container" option.  If a container is not specified for a task, "ubuntu:latest" is used by default; if a particular container is desired, it should be specified in the WDL, not on the command-line.
- add user-supplied "tag" to Run completion notification email
- add "--webhook" option for users to supply an URL, to which JAWS will post Run info upon completion
- deprecate /refdata dir -- specify reference files in your inputs json like any other input file and the JAWS file caching system will allow the file to be reused between multiple runs and will not be deleted until the files have not been accessed for some period of time (e.g. 14d).
- if the compute node has fast local disk (e.g. SSD) or large local disk (e.g. HD), then the /fast_scratch and/or /big_scratch volumes shall be mounted to the task's container
- run submissions are failed if any infile is missing (was warning); if your inputs.json refers to paths inside a container then use String instead of File for these
- make copy run inputs progress bar in human-readable units
- remove unnecessary file copy operations when compute-site is same as input-site (i.e. faster)
- `task-log` format has changed (was previously generated from JTM messages).  This also provides task status, in a more compact format.
- `task-status` was deprecated; use `task-log` instead.
- `task-summary` has been added, which provides more information about each task.

## 2.7.7 (2022-4-12) Summary
### User Facing Changes
#### Changed
- As before, the jaws-get command saves the wdl, json and zip files used for the run; however, the names have changed so all files are prefixed with "run_<run_id>"
- All jaws commands now use local time in the output, whereas before, some of the commands showed time in UTC.
- Ensure sharded tasks have unique names for task-log/task-status (#1064)
- Globus transfers will skip files when error when NFS problems instead of quit

#### Deprecated
- Deprecate wfcopy command (#1040)

#### Fixed
- Fix Cromwell-caching (#1059)
- WDL validation will now allow variable for memory (#1098)
- WDL validation will now allow sub workflow http import (#1097)
- Fix “get” command for list of output files (#1132)
#### Added
- Cromwell tasks echo HOSTNAME to stderr (#1044)
- Add stdout output to “errors” report (#1101)
- Allow users to submit subworkflows zip (#1089)
- Support for shifter container version tags (#1125)


### Internal Changes
#### Fixed
- Fixed end-to-end tests for jaws-get
- Fixed JTM to reuse workers that are already running
- Fixed the RabbitMQ settings in the CI/CD deployment for JTM service monitoring from TAHOMA

#### Added
- Reviewed and updated the documentation for dockerization of JAWS, and how to setup a local dev instance.
- Added initial dockerization implementation:
- Build.sh script to build the docker images
- Push images into gitlab repository
- Started collecting performance data with [pagurus](https://github.com/tylern4/pagurus) on cori.
- Updated site daemon to add performance metrics to elasticsearch
- Separated elasticsearch metrics into dev,staging and prod envs
- Updated site daemon to add run information to elasticsearch
- Updated JTM to access  high-memory nodes (1.5TB) on TAHOMA

#### Changed
- Refactored data transfer to use factory class and  plugin architecture
- Added a prefix, “slurm_jid=” to the SLURM job id in the `reason` filed of the `task-log` command
- Added HTCondor backed configuration for Cromwell and updated CI/CD deployment for HTCondor
- Jaws-site cromwell submit file handles instead of paths, for AWS
- Update site daemon to add performance metrics to elasticsearch
- Separate elasticsearch metrics into dev,staging and prod envs
- Factor out the rpc and base image into a separate docker image to speed up builds
- Update CI/CD to support builds (optional currently)

## 2.7.6 (2022-4-1) Summary

### User Facing Changes
#### Fixed
- Fixed WDL comments in the runtime section causing errors.
- The documentation (https://jaws-docs.readthedocs.io/en/latest/Specifications/configuringJTM_in_wdls.html) matches the jaws command `jaws list-sites`.
#### Added
- Default docker image is ubuntu:20.04 (WDL command always runs inside a container)
#### Changed
- The source file `jaws-prod.sh` can be used in a script with `set -u`.

### Internal Changes
#### Fixed
- Workflows were crashing on tahoma because of bad permissions on the `inputs` and `refdata` dirs.
### Added
- Created a status bar-chart for tasks as part of the jaws dashboard
- LaunchDarkley implemented for jaws dashboard
#### Changed
- WDL validation now allows variables to be used in the runtime section (i.e memory, time, etc.).

## 2.5.0 (2021-10-21) Summary
This release includes changes necessary to deploy to an additional computing site, tahoma.  Also included are several bug fixes and new features to improve the user experience.

### Deployment
- A new jaws-site, "tahoma" is deployed at PNNL (!947, !948, !958, !968, !974, !975, !979)
- JAWS-Sites should now restart automatically after scheduled maintenance (#155).

### End-to-end Tests
- Improved integration tests (!968, !973, !984)

### Bug Fixes
- Getting the metadata or errors report for a run with many (e.g. 1000) subworkflows resulted in a timeout error because jaws-site would request the metadata for each subworkflow in a separate request to the Cromwell server, which in rare cases took longer than the jaws-central REST server was willing to wait for a response.  This has been fixed as Cromwell offers an option to retrieve all subworkflow metadata with the main workflow; however, this changes the structure of the metadata JSON document slightly: rather than tasks having a "subWorkflowId", they now contain "subWorkflowMetadata" and the metadata is nested within the task.  The errors report, which is a filtered metadata document has also changed slightly as a result (!980, !982, !986).
- A workflow task which created a tmpfile could not be transferred by Globus due to file permissions, causing the results of such a run to not be downloadable.  This is resolved by skipping such tmpfiles (!981).

### New Features
- The client `jaws get` command now retrieves only files tagged as "outputs" in the WDL, by default, skipping intermediate files (e.g. stdout, stderr, rc, etc.) that are often of no use to the user.  The complete cromwell-executions output can still be retrieved by using the `--complete` option (!949, !953, !954, !955, !956, !966).  We also now generate an "outputs.json" file for each run which contains the relative paths to each output file as well as any non-path (e.g. string, integer) outputs (!964, !965, !968, !970, !971, !976, !977).
- The task-log and task-status tables have been improved for readability.  The cromwell_run_id column (UUID string) that was previously used to distinguish between tasks of the same name in different subworkflow calls has been dropped and, instead, subworkflow names (which are human-readable labels) now prefix the task names.  Additionally the "attempts" column has been dropped as we report only the last attempt of each task and we have actually never seen cromwell make multiple attempts at running a task (!983, !991).
- A progress bar is now displayed when input and output files are copied by the client, so the user can see the copy progress and see if there is a problem with the file system (!963(
- Modification to jaws-auth for supporting keycloak and integration with the upcoming jaws-dashboard (!950, !951).

## 2.4.0 (2021-08-30) Summary
This release contains various bugfixes and deployment changes. It also
addresses a security issue associated with Cromwell deployments at HPC sites.
Admin commands were also added to view runs and queues of all users.

### Deployment
- JAWS readthedocs is deployed using the Gitlab CI/CD pipeline (!894)
- JAWS software is no longer installed in `/tmp` but on paths located on cluster parallel filesystems (!881, !886, !912)
- Client now deploys bash autocomplete for autocompletion of jaws cmds (!871, !913)
- Cromwell configuration changes to allow for URLs in inputs JSON file (!877)
- Cromwell and womtool JARS are retrieved using wget (!935)
- Cromwell binds to localhost to address security concerns (!909)
- Tahoma deployment initiated. To be completed at a later release (!905, !923)
- Versions of JAWS dependencies updated (!895, !897, !891)

### Bug fixes
- Fixed a bug in client for a function that converts UTC to localtime (!904)
- JAWS site catches and recognizes CromwellException (!896)
- JAWS site would create a new subprocess for every client message. This would create hundreds of
  processes. JAWS was switched to use a single RPC client. (!889, !890, !911, !916, !917)
- RMQ queue property changed to remove temporary queues (!884)
- State descriptions updated (!914)

### WDL and Inputs validation
- JAWS now validates a WDLs runtime (!782)
- JAWS can intelligently parse and validate filepaths specified in the inputs.json file (!899)

### JAWS Admin commands
- Admins can now list all active jobs and see the run history and queue of all users (!901, !910)

## 2.3.0 (2021-05-17) Summary
This is a major release that simplifies the command structure for JAWS.  All the subcommands under `jaws wdl` are eleminated. This means there is no interaction with the jaws catalog. Instead, WDLs are stored in two gitlab repositories:

1) [public WDLs](https://code.jgi.doe.gov/official-jgi-workflows/wdl-specific-repositories), maintained by the public.
2) [sub-WDLs](https://code.jgi.doe.gov/official-jgi-workflows/jgi-wdl-tasks) maintained by JAWS staff.

Also `jaws run` is simplified to `jaws`.  For example, `jaws run submit` becomes `jaws submit`. See jaws --help.

Output from JAWS runs are now retrieved by running the `jaws get` command which copies results to your local folder or by running `jaws status --verbose` which will show the path to saved results on the JAWS scratch (you should copy files to your local folder).

Also, you can turn off caching (jaws submit --no-cache), the jaws errors command is more comprehensive in what it covers, and there are various flags added to the subcommands (see below).

### Major CLI changes
- `jaws get` which copies Run output now excludes Cromwell inputs and tmp folders
- `jaws submit` will print a warning, but not fail, if a path-like input is in the inputs.json file but doesn't exist on the filesystem (i.e. if you had /opt/img/data in the inputs.json and it only exists in the docker container).
- `jaws errors` command now includes a) backend stderr, b) task runtime parameters, c) task-log error messages; stdout was removed
- `jaws errors` now includes errors from subworkflows
- `jaws errors` format was changed to multi-level dictionary for improved readability.
- add `--no-cache` option to `jaws submit` to disable Cromwell result-caching for a Run
- add `--tag` option to `jaws submit` so users may add metadata (e.g. title, external ID, description) to runs; these are included in the Run records (e.g. `jaws status`, `jaws history`)
- add `--site` filter to `jaws queue` and `jaws history` commands
- add `--result` filter to `jaws history` command
- simplified `jaws status`, `jaws queue`, `jaws history` output
- `jaws task-status` now provides real-time results
- add `cromwell_run_id` and `cromwell_job_id` columns to `task-status` output so as not to collapse status messages when task-names are not unique (e.g. in scatter or subworkflow)
- add `jaws cancel-all` command
- fix formatting of jaws cli-client help screens


### Minor core changes
- include user's original WDL and JSON file paths in the Run records (e.g. `jaws status`, `jaws history`)
- fixed run input validation bug where URLs were mistaken for file paths
- faster file staging when compute-site is same as submission-site for a run
- fixed bug where transitioning Run status was delayed when there were many Task Logs to process
- corrected `jaws cancel` return code (now 0 for success)


## 2.2.0 (2021-04-09)
Minor release includes Globus client authentication support. Users no longer need to activate an endpoint
which required users to have an account at a compute site. This now uses a shared endpoint and service account to
do Globus transfers on behalf of the user.

Addresses error messaging and bugs in the jaws_rpc module

Addresses the stuck in "uploading" issue (#797, #781, #790, #603, #407). Jobs were stuck in
uploading due to slow processing of the messages from RMQ. MR !702 adds improvement in message processing.

When running and observing the job using `jaws queue` the run would show the
status as `queued` even though running `jaws run task-status` would show that the
some tasks are running. MR !749 fixes this issue

Includes fixes the inconsistent status between run and task-status.

When running and observing the job using `jaws queue` the run would show the
status as `queued` even though running `jaws run task-status` would show that the
some tasks are running. MR !749 fixes this issue

### Major Core Changes
- Globus endpoints are now set to shared endpoints for NERSC, LRC and EMSL (!560, !586)
- Users no longer have to authenticate or activate a Globus endpoint. Instead we use app client credentials [model](https://globus-sdk-python.readthedocs.io/en/stable/examples/client_credentials.html?highlight=secret). (!610)
- Properly sets the transfer paths for shared endpoints. There is a mechanism where non-root globus paths require virtual relative paths. (!621, !627)
- No longer require Globus accounts/tokens for transfers (!664)
- WDLS and JSON sent to output directory for reproducibility (!667, !669)
- Input files are no longer symlinked to output directory but are copied over (!679)
- Cromwell caching is turned on (!569)

### Minor core changes
- RPC log messages were not written to their log files, fix includes adding logger to jaws_rpc (!725)
- Unbound local error with jaws_rpc fix (!726)
- Improve messaging of error with inputs. No longer rejects bad paths in inputs.json but warns the user (!712)

Changes include fixes to regressions that took place after updates to staging (!717, !716)
 - Change install path location of Central back to /opt/jaws (!717)
 - Add log rotation back to JTM  (!716)
 - Rsync uses flags `-rLtq` and also allows `chmod` of the input files (!707)
 - Moves origin unmodified JSON to the output directory (!710)
 - Check permissions of wdl and inputs json in cromwell output directory (!709)

### Deployment changes
- Staging directory is now created by deployment scripts and has setgid set for genome group. Output directories are also created by deployment script (!562)
- Change software installation of JGI and Central away from `/tmp` (!642, commit hash `6f00ca1f`)
- Use the same directory as Globus uploads. Gets rid of "staging" diretory (commit hash `ddc65562343ef6d568917f4033c35486e1d47df9`)

### Database schema changes
- Database schema includes wdl and json file paths (!643)
- User no longer requires globus tokens. No longer necessary in schema (!664)
- Label column added so users can specify "label" their runs (!643)

### CLI changes
- User no longer specifies the output directory, instead there is a shared output directory where JAWS will transfer
data files upon completion. (!562)
- Users can retrieve files using the new `jaws run get` command. This allows users to set the permissions to the original
owner who submitted rather than owned by the `jaws` service account user. (!638)


## v2.1.1 (2020-11-02)
- Hotfix change to fix JAWS production deployment to production

## v2.1 (2020-11-02)
Major feature release that includes JGI Cloud support

### Major Core Changes
 - No longer reformat output automatically; instead users may use "wfcopy" command to reformat the output after their run completes.  Users may find it useful for comparison to runs performed by Cromwell outside of JAWS.
 - The output of tasks are returned as the tasks complete, rather than waiting for the end.  As the output includes the stdout/stderr files, the "run outputs" command has been deprecated and an new "run errors" command summarizes any errors (replaces "run output --failed").
 - "run task-log" and "run task-status" are now real-time; previously were out of sync by up to 10 seconds due to update interval. (!495)
 - Simplified user config file, provided a .sh file to source for activating jaws which simplifies use of multiple jaws deployments, and provide wheel file if you wish to install your own client or use it in your own python software.
 - Add "info" command to provide jaws version and deployment user client is using, as well as the link to the documentation appropriate for that release.
 - "run list-sites" now includes the maximum requestable RAM available at each Site
 - "run metadata" now returns any subworkflows' metadata too.
 - assorted minor bugfixes and improvements to usability.

### Deployment changes
- JAWS deploys to three different environments - dev, staging, and prod. It utilizes bash scripts to automate deployment.
- Simplify .gitlab-ci.yml by using environment variables to set configuration for different environments
- Due to instability with Spin, Central, central-RabbitMQ and central-MySQL were moved off of Spin and onto LBL IT SVM.

### Documentation
- Added instructions on activating lrc endpoint via Globus
- Added Known Issues section for users
- Let users know of task-log issue where there is a delay.


## v2.0 (2020-06-30)
Major feature release with new Gitlab integration and central gitlab repository

### Major core features
- Merging of JAWS services into a single repository
- Testing added (!6, !53, !73, )
- Integration with gitlab CI/CD deployment pipeline and usage of gitlab-runners on Cori, JGI, and AWS (!123, !121, !126, !127, !137)
- Packaging of JAWS services using Python wheels for deployment
- Automated deployment using shims and supervisord
- Shared RPC module is a python package and dependency to Site, Central, JTM (!191)
- Site and JTM communication path is done via RPC (!197)
- Central and Site no longer share a database (!202)


### Refactoring
- Central refactored to use setup.py for installation, follow PEP8 standard, logging and file rotation, mysql connection bugfixes, deprecated functions removed (!29)
- Auth to use logging, config class singleton, and to be setuptools friendly (!31)

### Documentation
- Modification to how JAWS works section (!46)
- Major changes to quickstart guide to reflect 2.0 usage (!49)
- Update of JAWS commands in documentation (!136)

### Controlled Vocabulary
All changes should be under one of these categories

`Added` for new features.
`Changed` for changes in existing functionality.
`Deprecated` for soon-to-be removed features.
`Removed` for now removed features.
`Fixed` for any bug fixes.
`Security` in case of vulnerabilities.
