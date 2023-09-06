# JAWS Template files
This directory contains template files that are used to create
shims/configs for JAWS site services. These files are created using
`envsubst` and environment variables that set by the `configs/` directory,
`setup_environment.sh` and the `all.sh` scripts.

## Templates
- **all.sh**: Script that creates all the compound paths, like the inputs/outputs directory for JAWS data files.
- **jaws-perlmutter-gitlab-runner.sh**: The shim created to run a gitlab-runner on Perlmutter using scrontab.
- **jaws-site.service.templ**: Template used for creating systemd service files
- **jaws-site.sh**: Template for running jaws-site in a scronjob 
- **setup_facl.sh**: Template for setfacl script 
- **starter.sh**: Template for script to submit JAWS site to scronjob.
- **site.env.templ**: Template for environment file to be used with Apptainer. Contains sensitive information like passwords.

### Virtual environment shims
These shims use a Python virtual environment to run the JAWS services. The virtual environment
is usually installed by the CI/CD pipeline.

- **perf-metrics.sh**: Shim for running performance metrics to elasticsearch
- **rpc-server.sh**: Shim for running JAWS RPC server
- **runs.sh**: Shim for running the run-daemon 
- **task-log.sh**: Shim for running the task-log daemon
- **transfer.sh**: Shim for running the transfer-daemon

### Supervisor templates
- **supervisor.conf**: General configuration file for supervisord.
- **supervisor.site.conf**: Template for running JAWS Site services.
- **supervisorctl.sh**: Shim for running supervisorctl.
- **supervisord.sh**: Shim for running supervisor daemon.

### Container shims and templates
The `container_runtime_templates` contains templates for shims to run the specified container runtime (e.g. apptainer).
These files are shims that are used to run services with apptainer. It takes care of the
paths needed for mounts as well as any arguments to run JAWS in its respective environment (eg LOG_LEVEL=DEBUG)
- **apptainer-instance.sh**: Shim for running JAWS services using `apptainer instance`.
- **apptainer-run.sh**: Shim for running JAWS services using `apptainer run`