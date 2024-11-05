import pytest
from unittest.mock import MagicMock, patch
from sqlalchemy.exc import IntegrityError, OperationalError, SQLAlchemyError
from datetime import datetime

from jaws_site import models
from jaws_site.task_logger import TaskLogger


ACTION = None


class JawsError(Exception):
    """Base class for all JAWS exceptions"""

    pass


class JawsUnavailableError(JawsError):
    """
    Base class for all recoverable errors.  A required service is temporarily unavailable; so while the action
    could not be performed, it is expected to be possible later.
    """

    pass


class JawsDbUnavailableError(JawsUnavailableError):
    """Db connection temporarily unavailable."""

    pass


def test_save(monkeypatch):
    def mock_insert(self, **kwargs):
        global ACTION
        ACTION = "INSERT"

    def mock_update(self, **kwargs):
        global ACTION
        ACTION = "UPDATE"

    monkeypatch.setattr(TaskLogger, "_insert", mock_insert)
    monkeypatch.setattr(TaskLogger, "_update", mock_update)

    params = {
        "timestamp": "2023-11-29 19:30:48",
    }

    mock_config = {}
    mock_session = None
    tl = TaskLogger(mock_config, mock_session)

    # test 1:
    params = {"timestamp": "2023-11-29 19:30:48", "status": "queued"}
    tl.save(params)
    assert ACTION == "INSERT"

    # test 2:
    params = {"timestamp": "2023-11-29 19:30:48", "status": "running"}
    tl.save(params)
    assert ACTION == "UPDATE"

    # test 3:
    params = {"timestamp": "2023-11-29 19:30:48", "status": "done"}
    tl.save(params)
    assert ACTION == "UPDATE"

    # test 4:
    params = {"timestamp": "2023-11-29 19:30:48", "status": "INVALID"}
    tl.save(params) == False


def test_delta_minutes():
    mock_config = {}
    mock_session = None

    start = datetime.strptime("2023-11-29 19:30:48", "%Y-%m-%d %H:%M:%S")
    end = datetime.strptime("2023-11-30 19:30:48", "%Y-%m-%d %H:%M:%S")

    tl = TaskLogger(mock_config, mock_session)
    actual = tl.delta_minutes(start, end)
    assert actual == 24 * 60


@pytest.fixture
def mock_session():
    """Fixture to create a mock SQLAlchemy session."""
    return MagicMock()


@pytest.fixture
def config():
    """Fixture to provide a mock configuration."""
    return {"some_config_key": "some_config_value"}


@pytest.fixture
def task_logger(config, mock_session):
    """Fixture to create a TaskLogger instance."""
    return TaskLogger(config=config, session=mock_session)


# def test_save_queued_task(task_logger, mock_session):
#     """Test saving a 'queued' task log message."""
#     params = {
#         "timestamp": "2023-10-05 12:34:56",
#         "cromwell_run_id": "run123",
#         "task_dir": "task_dir123",
#         "status": "queued",
#     }

#     result = task_logger.save(params)

#     assert result is True
#     mock_session.add.assert_called_once()
#     mock_session.commit.assert_called_once()


def test_save_running_task(task_logger, mock_session):
    """Test saving a 'running' task log message."""
    params = {
        "timestamp": "2023-10-05 12:34:56",
        "cromwell_run_id": "run123",
        "task_dir": "task_dir123",
        "status": "running",
    }

    mock_session.query().filter().order_by().first.return_value = models.Tasks(
        cromwell_run_id="run123", task_dir="task_dir123", queue_start=datetime.now()
    )

    result = task_logger.save(params)

    assert result is True
    mock_session.commit.assert_called_once()


def test_save_invalid_status(task_logger):
    """Test saving a task log message with an invalid status."""
    params = {
        "timestamp": "2023-10-05 12:34:56",
        "cromwell_run_id": "run123",
        "task_dir": "task_dir123",
        "status": "invalid_status",
    }

    result = task_logger.save(params)

    assert result is False


# def test_insert_operational_error(task_logger, mock_session):
#     """Test handling an OperationalError during insert."""
#     params = {
#         "timestamp": datetime.now(),
#         "cromwell_run_id": "run123",
#         "task_dir": "task_dir123",
#         "status": "queued",
#     }

#     mock_session.query().scalar.return_value = False
#     mock_session.add.side_effect = OperationalError("mock", "mock", "mock")

#     with pytest.raises(JawsDbUnavailableError):
#         task_logger._insert(**params)


# def test_update_operational_error(task_logger, mock_session):
#     """Test handling an OperationalError during update."""
#     params = {
#         "timestamp": datetime.now(),
#         "cromwell_run_id": "run123",
#         "task_dir": "task_dir123",
#         "status": "running",
#     }

#     mock_session.query().filter().order_by().first.side_effect = OperationalError(
#         "mock", "mock", "mock"
#     )

#     with pytest.raises(JawsDbUnavailableError):
#         task_logger._update(**params)


def test_delta_minutes():
    """Test the delta_minutes method."""
    start = datetime(2023, 10, 5, 12, 0, 0)
    end = datetime(2023, 10, 5, 12, 30, 0)

    result = TaskLogger.delta_minutes(start, end)

    assert result == 30


def test_task_exists(task_logger, mock_session):
    """Test the _task_exists method."""
    cromwell_run_id = "run123"
    task_dir = "task_dir123"

    mock_session.query().scalar.return_value = True

    result = task_logger._task_exists(cromwell_run_id, task_dir)

    assert result is True
