from jaws_site import config
from jaws_site.cromwell import Cromwell
from jaws_site.tasks import get_task_log_error_messages


def _get_cromwell_errors_report(cromwell_run_id):
    """Get errors report from Cromwell metadata.
    :param cromwell_run_id: Cromwell's UUID for the run
    :type cromwell_run_id: str
    :return: errors for each failed task
    :rtype: dict
    """
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    metadata = cromwell.get_metadata(cromwell_run_id)
    return metadata.errors()


def get_errors(session, cromwell_run_id):
    """Retrieve error messages from Cromwell Metadata and TaskLog.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: error report
    :rtype: dict
    """
    errors_report = _get_cromwell_errors_report(cromwell_run_id)
    for task_name in errors_report:
        if task_name == cromwell_run_id:
            continue
        cromwell_job_id = errors_report[task_name]["cromwell_job_id"]
        task_log_err_msgs = get_task_log_error_messages(session, cromwell_job_id)
        errors_report[task_name]["task-log"] = "\n".join(task_log_err_msgs)
    return errors_report
