# Changelog

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
