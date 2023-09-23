from tests.conftest import MockSession, MockLogger
from jaws_site.tasks import TaskLog
from datetime import datetime


DATETIME_FMT = "%Y-%m-%d %H:%M:%S"


def test_did_run_start(monkeypatch):
    # test 1: yes
    def mock__select_num_started(self):
        return 1

    monkeypatch.setattr(TaskLog, "_select_num_started", mock__select_num_started)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is True

    # test 2: No
    def mock__select_num_started(self):
        return 0

    monkeypatch.setattr(TaskLog, "_select_num_started", mock__select_num_started)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    assert task_log.did_run_start() is False


def test_table(monkeypatch):
    def mock__select_all_rows(self):
        return [
            [
                123,
                "ABCD-EFGH-IJKL-MNOP",
                "2421",
                "call-do_something",
                "done",
                datetime.strptime("2023-04-24 11:00:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:01:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:03:00", DATETIME_FMT),
                1,
                2,
                0,
                None,
                None,
                None,
                None,
                None,
            ],
        ]

    expected = {
        "header": [
            "TASK_DIR",
            "STATUS",
            "QUEUE_START",
            "RUN_START",
            "RUN_END",
            "RC",
            "QUEUE_MINUTES",
            "RUN_MINUTES",
            "CACHED",
            "TASK_NAME",
            "REQ_CPU",
            "REQ_GB",
            "REQ_MINUTES",
        ],
        "data": [
            [
                "call-do_something",
                "done",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
                0,
                1,
                2,
                None,
                None,
                None,
                None,
                None,
            ]
        ],
    }

    monkeypatch.setattr(TaskLog, "_select_all_rows", mock__select_all_rows)
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
    task_log = TaskLog(
        mock_session, mock_cromwell_run_id, mock_logger, local_tz="Africa/Algiers"
    )

    timestamp = datetime.strptime("2023-04-24 08:00:00", "%Y-%m-%d %H:%M:%S")
    actual = task_log._utc_to_local_str(timestamp)
    expected = "2023-04-24 09:00:00"
    assert actual == expected


def test_time_minutes():
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    minutes = task_log.time_minutes("01:02:00")
    assert minutes == 62


def test_memory_gb():
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    gb = task_log.memory_gb("1 TB")
    assert gb == 1024
