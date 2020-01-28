# JAWS Central

Central server for JGI Analysis Workflows Service (JAWS).  Command-line clients (see: JAWS Client) interact (only) with this server via REST and JSON requests.

## Services Provided:

This server provides two services:

### Workflows Catalog

Workflow specifications (i.e. WDL files and associated markdown README files) are stored in a relational database.  Workflows are versioned and are owned by one or more users.  Workflow versions which have be tagged as "released" (i.e. production releases) are immutable, but unreleased (i.e. development versions) may be updated by owners.

### Analysis Management

This server allows users to add, query, and delete analysis runs.  Consequently, the server manages the analysis database, which maps analysis-runs to the computing-site to which they have been submitted, so it can relay requests to the appopriate site via JSON-RPC2 (see: JAWS Site).  Besides routing (e.g. submit analysis to site-X, query status of analysis Y), it may provide some basic aggregation of results (e.g. list all running analyses).

## Related JAWS Services

### JAWS Dashboard

All webpages are provided by the JAWS Dashboard Flask server, which runs independently of this service.

### JAWS Auth

There is also a separate JAWS OAuth2 service which supports authentication for all JAWS services, including this one.

### JAWS Site

Each computing facility has a JAWS Site server for managing the analyses running there.  This is an AmqpStorm-RPC service using JSON-RPC2 messages and provides an interface to the Cromwell servers running there.  Cromwell, in turn, submits tasks to the JGI Task Manager (JTM) running at each site, which manages the jobs on the compute cluster, which is usually managed by a scheduler (e.g. SLURM, SGE, UGE) but not always (e.g. AWS, Google Cloud).
