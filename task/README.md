# jaws-task service

This service stores and serves all task-related data.

 * Cromwell creates Tasks via the jaws-backend RPC client
 * jaws-task submits to job scheduler (e.g. JTM, Slurm, Google Cloud)
 * jaws-task may receive updates on job status from job scheduler
 * jaws-central may request job log or status
