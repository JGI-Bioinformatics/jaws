import pytest
from datetime import datetime
from unittest.mock import Mock, patch
from pathlib import Path

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
                100,
                100,
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
            "INPUT_SIZE",
            "OUTPUT_SIZE",
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
                100,
                100,
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


# Test for _get_return_code
@pytest.fixture
def task_logger():
    mock_session = Mock()
    mock_cromwell_run_id = "test_run_id"
    mock_logger = Mock()
    return TaskLog(mock_session, mock_cromwell_run_id, logger=mock_logger)


@pytest.fixture
def temp_dir(tmp_path):
    return tmp_path


def test_get_return_code_valid(task_logger, temp_dir):
    rc_file = temp_dir / "execution" / "rc"
    rc_file.parent.mkdir(parents=True)
    rc_file.write_text("0")

    result = task_logger._get_return_code(str(temp_dir))
    assert result == 0


def test_get_return_code_non_zero(task_logger, temp_dir):
    rc_file = temp_dir / "execution" / "rc"
    rc_file.parent.mkdir(parents=True)
    rc_file.write_text("1")

    result = task_logger._get_return_code(str(temp_dir))
    assert result == 1


def test_get_return_code_large_number(task_logger, temp_dir):
    rc_file = temp_dir / "execution" / "rc"
    rc_file.parent.mkdir(parents=True)
    rc_file.write_text("1000000")  # A large number that should still be valid

    result = task_logger._get_return_code(str(temp_dir))
    assert result == 1000000


def test_get_return_code_negative_number(task_logger, temp_dir):
    rc_file = temp_dir / "execution" / "rc"
    rc_file.parent.mkdir(parents=True)
    rc_file.write_text("-1")

    result = task_logger._get_return_code(str(temp_dir))
    assert result == -1
