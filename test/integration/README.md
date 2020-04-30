# Integration Testing

## Architecture (high level)

    +---------+          +------+             +-----------------+
    | central | <------- | site | <---------  | compute cluster |
    +---------+ <--+     |      | --------->  |                 |
                   |     +------+             +-----------------+
                   |        .
                   |        .
                   |        .
                   |        .
                   |     +------+
                   |     | site |
                   +---- |      |
                         +------+

## Site Preparation Guide

This section describes how to deploy JAWS at a new site. This is not exclusively applicable to
integration testing and can be used to prepare a future JAWS site. Site is an overloaded term
and is context dependent: organizationally it is the site that JAWS runs at (eg LBNL), the
whole of the site service hosted at that organization and the application part of the
site service.

### Site Overview

     to central
       ^
       |                         +-----------> to worker filesystem
       |                         |
    +----------------------+-----+------+---------+
    | Site  |              |            |         |
    |  App  |<------------>| Cromwell   |         |
    |-------+              |            |         |
    |                      +------------+         |
    |                            ^        +-------+
    |                            |  +---->|  RMQ  |----> to worker
    |                            v  v     |       |
    |                       +----+------+ +-------+
    |                       |   JTM     |         |
    |                       |           |         |
    +-----------------------------------+---------|
                                 |
                                 |
                                 +---> to scheduler

[Cromwell](https://github.com/broadinstitute/cromwell) is a workflow management engine developed by the Broad Institute, with focus on bioinformatics workflows.

### Network Connectivity

#### Site to External

The JAWS-site service needs to be able to communicate with the JAWS-central service. All connections
are outgoing and are initiated by the site to the target system:

- site to RabbitMQ broker for RPC between site and central (eg site to rmq.lbl.gov, default port 5672)
- site to central database for additional RPC between site and central (eg site to jaws-db.lbl.gov, default port 3306)

JAWS is deployed in an automated fashion, using the Gitlab CI system, using a "gitlab-runner". The current
implementation requires all components of JAWS-site being deployed by the same runner.

- runner (on site host) to code.jgi.doe.gov, port 443

#### Site to Internal

Site submits workflow runs to Cromwell, which submits tasks to JTM-Manager (JGI Task Manager, Manager process) via RabbitMQ message.  The Manager submits jobs to a compute cluster but those jobs do not contain specific tasks to run, but JTM-Workers instead.  The Workers consume tasks from the appropriate RabbitMQ queue, with queues for different sized workers (with respect to RAM, CPU resources).  The RabbitMQ instance can be hosted beside the JAWS-site service, or standalone. In any case the site-service and the worker running on the compute cluster need to be able to communicate with the rabbitmq instance.

- site to rabbitmq broker (eg. jtm-rmq.newsite.org)
- cromwell to rabbitmq broker (eg. jtm-rmq.newsite.org)
- worker to rabbitmq broker (eg.  jtm-rmq.newsite.org)

### Filesystem Access

Cromwell, JTM, and Site on need to be able to read/write the same filesystem as the worker nodes (ie. cromwell-executions folder).  Cromwell prepares the directory structure for a task, JTM-Worker(s) generate the task output files, and after all tasks complete, Site reorganizes the datafiles and writes them to a results folder, which is accessible from the Site's Globus endpoint.

### Globus Endpoint Requirements

The JAWS application is registered as a Globus app and has a unique client ID.  Users must have their own Globus accounts and authorize the JAWS app to peform transfers on their behalf.  The storage resources where input files are stored and output files are desired, must have a Globus endpoint capable of reading and writing, respectively.  Additionally, the computing-site must have a Globus endpoint.  The basedir of the Globus endpoint may not be root, and must be specified in the JAWS config, as must the staging and results subdirectories, which must be under the endpoint basedir.  The cromwell-executions folder does not have to be under the Globus endpoint basedir as Site copies the files in/out of the staging and results directories.

### Accounts

For privilege separation between the user-supplied code and JAWS services, two Linux user accounts are necessary:
an account running the services (eg. jaws-usr) and another one executing user codes on the cluster (eg. jtm-usr).  The jaws-usr belongs to the jaws and jtm groups, but the jtm-usr belongs only to the jtm group.

    +------------------------------------------------------+
    | Service (eg. jaws-svc)  |    User (eg. jaws-usr)     |
    |-------------------------|----------------------------+
    |  rabbitmq               |  jaws-usr                  |
    |  database               |  jaws-usr                  |
    |  jaws-site              |  jaws-usr                  |
    |  jtm-manager            |  jtm-usr                   |
    |  jtm-worker             |  jtm-usr  >----------------+--------> Slurm
    +------------------------------------------------------|

The jtm-usr needs to be able to access Slurm, to submit jobs, monitor their status, etc.
The jaws-usr needs to be able to access Globus, to transfer files.

### Compute Node Requirements


TODO:

- Container Runtime (Shifter, Singularity)
- Python3 interpreter

