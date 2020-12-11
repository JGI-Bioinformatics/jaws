task_status_msg = {
    "ready": "The job has been prepared by Cromwell and submitted to JTM",
    "queued": "The job was received by JTM-manager and sent to JTM-worker",
    "pending": "The job was receive by JTM-worker and is awaiting resources",
    "running": "The job is currently executing",
    "success": "The job completed successfully",
    "failed": "The job has failed",
    "outofresource": "The job exceeeded the reserved RAM; increase the amount in the WDL and rerun",
    "terminated": "The run was terminated",
    "invalidtask": "The task definition in the message from jtm_submit is not valid",
    "timeout": "The worker timed out",
    "lostconnection": "The manager lost connection with the worker",
}
