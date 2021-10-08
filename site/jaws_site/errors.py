import logging
from jaws_site import config
from jaws_site.cromwell import Cromwell, CromwellError
from jaws_site.tasks import get_task_log_error_messages


logger = logging.getLogger(__package__)


def _get_cromwell_errors_report(cromwell_run_id):
    """Get errors report from Cromwell metadata.
    :param cromwell_run_id: Cromwell's UUID for the run
    :type cromwell_run_id: str
    :return: errors for each failed task
    :rtype: dict
    """
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    return cromwell.get_all_errors(cromwell_run_id)


def get_errors(session, cromwell_run_id):
    """Retrieve error messages from Cromwell Metadata and TaskLog.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: error report
    :rtype: dict
    """
    try:
        full_errors_report = _get_cromwell_errors_report(cromwell_run_id)
    except CromwellError as error:
        logger.error(f"Failed to retrieve Cromwell errors report for {cromwell_run_id}: {error}")
        return {}
    for a_cromwell_run_id, errors_report in full_errors_report.items():
        if "calls" in errors_report:
            for task_name in errors_report["calls"]:
                for item in errors_report["calls"][task_name]:
                    if "jobId" in item:
                        cromwell_job_id = item["jobId"]
                        item["taskLog"] = get_task_log_error_messages(session, cromwell_job_id)
    return full_errors_report
