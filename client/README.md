# JAWS Client

Command-line client for JGI Analysis Workflows Service (JAWS)

## Terms

- workflow : analysis pipeline specification in WDL (or CWL) format
- task : a workflow is comprised of many tasks
- job : a reservation on computing resources for performing a task
- inputs file : JSON file with parameters (e.g. paths) required by the workflow
- analysis run (AKA analysis AKA run) : execution of a workflow with specific inputs file
- jaws-central : the REST server the CLI interacts with
- jaws-site(s) : RPC servers at each computing cluster
- workflow manager (AKA Cromwell) : perform the run by submitting tasks to the task manager 
- task manager (AKA JTM) : middleware between the workflow manager and the scheduler (e.g. SLURM) associates tasks to jobs
- scheduler (e.g. SLURM) : controls the allocation of shared computing resources to run jobs
- catalog : a managed repository of shared workflows

## Developer Documentation

### Packaging

The resulting wheel is located in the "dist/" folder. The version of the wheel via "git describe" (see setup.py).
This requires the python package "wheel" to be installed.
A new version can be created by tagging the source tree accordingly:

```
  git tag 1.0.0
  make dist
```

To make development easier, setuptools provides a development mode, which obviates the need to install a package,
when working out of the git repo directly. Changes to the source tree will be immediately active.

```
  python setup.py develop
```

### Versioning

This project follows the [semantic versioning](https://semver.org/) spec.
