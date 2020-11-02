# Changelog

## v2.1 (2020-11-02)
Major feature release that includes JGI Cloud support 

### Major Core Changes
 - No longer reformat output automatically; instead users may use "wfcopy" command to reformat the output after their run completes.  Users may find it useful for runs performed by Cromwell, outside of JAWS.
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