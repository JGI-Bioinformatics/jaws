import pytest
from datetime import datetime, timezone

from jaws_site.tasks import TaskLog

from tests.conftest import MockLogger, MockSession

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
    def mock_select(self):
        return [
            [
                9999,
                "call-do_something",
                "14999",
                "done",
                datetime.strptime("2023-04-24 11:00:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:01:00", DATETIME_FMT),
                datetime.strptime("2023-04-24 11:03:00", DATETIME_FMT),
                1,
                2,
                None,
                None,
                None,
                None,
                None,
                0,
            ],
        ]

    expected = {
        "header": [
            "TASK_DIR",
            "JOB_ID",
            "STATUS",
            "QUEUE_START",
            "RUN_START",
            "RUN_END",
            "QUEUE_MIN",
            "RUN_MIN",
            "CACHED",
            "TASK_NAME",
            "REQ_CPU",
            "REQ_GB",
            "REQ_MIN",
            "CPU_HRS",
            "RETURN_CODE",
        ],
        "data": [
            [
                "call-do_something",
                "14999",
                "done",
                "2023-04-24 11:00:00",
                "2023-04-24 11:01:00",
                "2023-04-24 11:03:00",
                1,
                2,
                None,
                None,
                None,
                None,
                None,
                None,
                0,
            ]
        ],
    }

    monkeypatch.setattr(TaskLog, "select", mock_select)
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    actual = task_log.table(local_tz="UTC")
    assert actual == expected


def test__utc_to_local_str(monkeypatch):
    def mock_select(self):
        return []

    monkeypatch.setattr(TaskLog, "select", mock_select)
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


def test_memory_gb():
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    gb = task_log.memory_gb("1 TB")
    assert gb == 1024


def test_cpu_hours(monkeypatch) -> float:
    def mock__select_all_cpu_minutes(self):
        return [
            [1, 60],
            [2, 100],
        ]

    monkeypatch.setattr(
        TaskLog, "_select_all_cpu_minutes", mock__select_all_cpu_minutes
    )

    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)
    cpu_hrs = task_log.cpu_hours()
    assert cpu_hrs == 4.33


def test_int_or_none():
    mock_session = MockSession()
    mock_logger = MockLogger()
    mock_cromwell_run_id = "ABCD-EFGH-IJKL-MNOP"
    task_log = TaskLog(mock_session, mock_cromwell_run_id, mock_logger)

    assert task_log.int_or_none(None) is None
    assert task_log.int_or_none([]) is None
    assert task_log.int_or_none({}) is None

    actual = task_log.int_or_none(int(1))
    assert type(actual) is int and actual == 1

    actual = task_log.int_or_none("2")
    assert type(actual) is int and actual == 2

    actual = task_log.int_or_none(3.1)
    assert type(actual) is int and actual == 3


@pytest.fixture
def task_log():
    return TaskLog(session=None, cromwell_run_id="test_id")


@pytest.mark.parametrize(
    "input_time,expected_output",
    [
        (datetime(2023, 1, 1, 12, 0, 0, tzinfo=timezone.utc), "2023-01-01 04:00:00"),
        (datetime(2023, 6, 1, 12, 0, 0, tzinfo=timezone.utc), "2023-06-01 05:00:00"),
        ("2023-04-24 08:00:00", "2023-04-24 01:00:00"),
        (None, None),
    ],
)
def test_utc_to_local_str(task_log, input_time, expected_output):
    result = task_log._utc_to_local_str(input_time)
    assert result == expected_output


def test_utc_to_local_str_custom_timezone(task_log):
    task_log.set_local_tz("Europe/London")
    input_time = datetime(2023, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    expected_output = "2023-01-01 12:00:00"
    result = task_log._utc_to_local_str(input_time)
    assert result == expected_output


def test_utc_to_local_str_invalid_string():
    task_log = TaskLog(session=None, cromwell_run_id="test_id")
    with pytest.raises(ValueError):
        task_log._utc_to_local_str("invalid datetime string")


def test_utc_to_local_str_invalid_type():
    task_log = TaskLog(session=None, cromwell_run_id="test_id")
    with pytest.raises(AttributeError):
        task_log._utc_to_local_str(123)  # Integer instead of datetime or string


def test_utc_to_local_str_daylight_savings():
    task_log = TaskLog(session=None, cromwell_run_id="test_id")
    task_log.set_local_tz("America/New_York")
    winter_time = datetime(2023, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    summer_time = datetime(2023, 7, 1, 12, 0, 0, tzinfo=timezone.utc)
    assert task_log._utc_to_local_str(winter_time) == "2023-01-01 07:00:00"
    assert task_log._utc_to_local_str(summer_time) == "2023-07-01 08:00:00"
