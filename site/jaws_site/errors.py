import logging
import copy
from jaws_site import config
from jaws_site.cromwell import Cromwell, CromwellError
from jaws_site.tasks import get_task_log_error_messages


logger = logging.getLogger(__package__)


def _get_cromwell_errors_report(cromwell_run_id):
    """
    Get errors report from Cromwell metadata.

    :param cromwell_run_id: Cromwell's UUID for the run
    :type cromwell_run_id: str
    :return: errors for each failed task
    :rtype: dict
    """
    cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
    metadata = cromwell.get_metadata(cromwell_run_id)
    return metadata.errors()


def _add_task_log_errors(session, report):
    """
    Recursive function to parse "calls" section of errors struct and add
    error messages from task log.

    :param session: db session handle
    :type session: sqlalchmey.session
    :param report: report dictionary with "calls"
    :type report: dict
    :return: new report dictionary with "calls" and task log errors added
    :rtype: dict
    """
    calls = report["calls"]
    new_calls = {}
    for task_name, task_calls in calls.items():
        # the errors report only contains the last attempt (i.e. exactly 1 call)
        new_calls[task_name] = []
        call = task_calls[0]
        new_call = {}
        if "subWorkflowMetadata" in call:
            # this is a subworkflow, so recurse
            new_call["subWorkflowMetadata"] = _add_task_log_errors(
                session, call["subWorkflowMetadata"]
            )
        else:
            # this is a regular task
            new_call = copy.deepcopy(call)
            if "jobId" in call:
                # the task was submitted to the back-end so check task-log for errors
                new_call["taskLog"] = get_task_log_error_messages(
                    session, call["jobId"]
                )
        new_calls[task_name].append(new_call)
    new_report = {"calls": new_calls}
    return new_report


def get_errors(session, cromwell_run_id):
    """
    Retrieve error messages from Cromwell Metadata and TaskLog.

    :param cromwell_run_id: Cromwell run ID
    :type params: dict
    :return: error report
    :rtype: dict
    """
    try:
        errors_report = _get_cromwell_errors_report(cromwell_run_id)
    except CromwellError as error:
        logger.error(
            f"Failed to retrieve Cromwell errors report for {cromwell_run_id}: {error}"
        )
        return {}

    # if there are no errors, there will not be a "calls" section
    if "calls" not in errors_report:
        return errors_report

    # generate new errors report with any task-log error messages added
    new_report = _add_task_log_errors(session, errors_report)
    return new_report
