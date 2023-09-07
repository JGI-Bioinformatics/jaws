from tests.conftest import MockSession, MockLogger
from datetime import datetime
from jaws_site.tasks import TaskLog


def test_did_run_start(monkeypatch):
    # test 1: yes
    def mock__select_rows(self):
        return [
            [
                "call-count_seqs",
                "2023-04-03 17:20:35",
                "2023-04-03 17:20:52",
                "2023-04-03 17:20:55",
                0,
                None,
                "0:0:17",
                "0:0:03",
            ],
        ]

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is True

    # test 2: No
    def mock__select_rows(self):
        return [
            [
                "call-count_seqs",
                "2023-04-03 17:20:35",
                None,
                None,
                None,
                None,
                None,
                None,
            ],
        ]

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is False


def test_table(monkeypatch):
    def mock__select_rows(self):
        return [
            [
                "call-do_something",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
                0,
                None,
                "0:01:00",
                "0:02:00",
            ],
        ]

    expected = {
        "header": [
            "TASK",
            "STATUS",
            "QUEUED",
            "RUNNING",
            "COMPLETED",
            "RC",
            "CANCELLED",
            "QUEUE_DUR",
            "RUN_DUR",
        ],
        "data": [
            [
                "call-do_something",
                "completed",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
                0,
                None,
                "0:01:00",
                "0:02:00",
            ]
        ],
    }

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    actual = task_log.table(local_tz="GMT")
    assert actual == expected


def test__utc_to_local(monkeypatch):
    def mock__select_rows(self):
        return []

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    actual = task_log._utc_to_local("2023-04-24 08:00:00", "US/Pacific")
    expected = "2023-04-24 01:00:00"
    assert actual == expected
