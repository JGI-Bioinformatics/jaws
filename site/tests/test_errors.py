from jaws_site import errors, tasks
from deepdiff import DeepDiff


def test_get_errors(monkeypatch):

    example_cromwell_run_id = "AAAA-BBBB-CCCC"
    example_cromwell_errors_report = {
        "calls": {
            "example_sub": [
                {
                    "subWorkflowMetadata": {
                        "calls": {
                            "example_task": [
                                {
                                    "jobId": "300",
                                }
                            ],
                        }
                    }
                }
            ],
        }
    }
    example_task_log_errors = {
        "300": ["Out of time!"],
    }
    expected_errors_report = {
        "calls": {
            "example_sub": [
                {
                    "subWorkflowMetadata": {
                        "calls": {
                            "example_task": [
                                {
                                    "jobId": "300",
                                    "taskLog": ["Out of time!"],
                                }
                            ],
                        }
                    }
                }
            ],
        }
    }

    def mock_get_cromwell_errors_report(cromwell_run_id):
        return example_cromwell_errors_report

    monkeypatch.setattr(
        errors, "_get_cromwell_errors_report", mock_get_cromwell_errors_report
    )

    def mock_select_task_log_error_messages(session, cromwell_job_id):
        if cromwell_job_id in example_task_log_errors:
            return example_task_log_errors[cromwell_job_id]
        else:
            return []

    monkeypatch.setattr(
        tasks, "_select_task_log_error_messages", mock_select_task_log_error_messages
    )

    mock_session = None
    actual_errors_report = errors.get_errors(mock_session, example_cromwell_run_id)
    assert (
        bool(DeepDiff(actual_errors_report, expected_errors_report, ignore_order=True))
        is False
    )
