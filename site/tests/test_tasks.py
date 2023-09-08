from tests.conftest import MockSession, MockLogger
from jaws_site.tasks import TaskLog
from datetime import datetime


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"


def test_did_run_start(monkeypatch):
    # test 1: yes
    def mock__select_rows(self):
        return [
            [
                "call-count_seqs",
                "done",
                datetime.strptime("2023-04-03 17:20:35", DATETIME_FMT),
                datetime.strptime("2023-04-03 17:20:52", DATETIME_FMT),
                datetime.strptime("2023-04-03 17:20:55", DATETIME_FMT),
                0,
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
                "queued",
                datetime.strptime("2023-04-03 17:20:35", DATETIME_FMT),
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
                "done",
                datetime.strptime("2023-04-24 11:00:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:01:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:03:00", DATETIME_FMT),
                0,
                "0:01:00",
                "0:02:00",
            ],
        ]

    expected = {
        "header": [
            "TASK",
            "STATUS",
            "QUEUE_START",
            "RUN_START",
            "RUN_END",
            "RC",
            "QUEUE_DUR",
            "RUN_DUR",
        ],
        "data": [
            [
                "call-do_something",
                "done",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
                0,
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
    actual = task_log.table(local_tz="UTC")
    assert actual == expected


def test__utc_to_local_str(monkeypatch):
    def mock__select_rows(self):
        return []

    monkeypatch.setattr(TaskLog, "_select_rows", mock__select_rows)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger, local_tz="Africa/Algiers")

    timestamp = datetime.strptime("2023-04-24 08:00:00", "%Y-%m-%d %H:%M:%S")
    print(type(timestamp))  # ECCE
    actual = task_log._utc_to_local_str(timestamp)
    expected = "2023-04-24 09:00:00"
    assert actual == expected
