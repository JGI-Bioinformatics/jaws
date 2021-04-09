# Changelog

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
- Users can retrieve files using the new `jaws get` command. This allows users to set the permissions to the original
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