from tests.conftest import MockSession, MockLogger
from datetime import datetime
from jaws_site.tasks import TaskLog


def test_did_run_start(monkeypatch):
    # test 1: yes
    def mock__select_rows(self):
        self.data = [
            ["call-count_seqs", "queued", "2023-04-03 17:20:35"],
            ["call-count_seqs", "running", "2023-04-03 17:20:52"],
            ["call-count_seqs", "succeeded", "2023-04-03 17:20:55"],
        ]

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is True

    # test 2: No
    def mock__select_rows(self):
        self.data = [
            [
                "call-count_seqs",
                "queued",
                datetime.strptime("2023-04-03 17:20:35", "%Y-%m-%d %H:%M:%S"),
            ]
        ]

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is False


def test_table(monkeypatch):
    def mock__select_rows(self):
        self.data = [
            [
                "call-do_something",
                "queued",
                "2023-04-24 11:00:00",
            ],
            ["call-do_something", "running", "2023-04-24 11:01:00"],
            ["call-do_something", "succeeded", "2023-04-24 11:03:00"],
        ]

    expected = {
        "header": [
            "TASK",
            "STATUS",
            "QUEUED",
            "RUNNING",
            "FINISHED",
            "QUEUE_DUR",
            "RUN_DUR",
        ],
        "data": [
            [
                "call-do_something",
                "succeeded",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
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
    actual = task_log.table()
    assert actual == expected


def test__utc_to_local(monkeypatch):
    def mock__select_rows(self):
        self.data = []

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    actual = task_log._utc_to_local("2023-04-24 08:00:00", "US/Pacific")
    expected = "2023-04-24 01:00:00"
    assert actual == expected
