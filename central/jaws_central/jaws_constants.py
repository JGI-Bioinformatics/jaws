# states should use concise names, all lowercase, and uses spaces (no camelCase, CAPS, or under_scores)
run_status_msg = {
    "uploading": "Your run inputs are being sent to the compute site via Globus.",
    "upload failed": "The Globus transfer of your run to the compute-site failed.  Possible cause may be due to inactive globus enpoint.  Please go to https://app.globus.org/file-manager, on the left side of the page, select ENDPOINTS, click the > icon to the right of the NERSC DTN endpoint, then click Activate.",  # noqa
    "upload inactive": "Globus transfer stalled; try reactivating the endpoint. Please go to https://app.globus.org/file-manager, on the left side of the page, select ENDPOINTS, click the > icon to the right of the NERSC DTN endpoint, then click Activate.",  # noqa",
    "upload complete": "Your run inputs have been transferred and are ready to submit to Cromwell",
    "missing input": "The run was uploaded but some of the required files were missing..",
    "submitted": "The run has been submitted to Cromwell and tasks should start to queue within moments",
    "submission failed": "The run was submitted to Cromwell but rejected due to invalid input.",
    "queued": "At least one task has requested resources but no tasks have started running yet",
    "running": "The run is being executed; you can check `task-status` for more detail",
    "succeeded": "The run has completed successfully and is waiting for the output to be downloaded",
    "failed": "The run has failed; see: `output --failed` and `metadata` for more detail",
    "aborting": "Your run is in the process of being canceled",
    "aborted": "The run was cancelled",
    "downloading": "The run output is being sent via Globus",
    "download failed": "Globus failed to return the results to the user",
    "download inactive": "Globus transfer stalled; try reactivating the endpoint. Please go to https://app.globus.org/file-manager, on the left side of the page, select ENDPOINTS, click the > icon to the right of the NERSC DTN endpoint, then click Activate.",  # noqa",
    "download complete": "The run output (suceeded or failed) has been returned to the user.",
}
