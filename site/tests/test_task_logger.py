from datetime import datetime

import pytest
from jaws_site.task_logger import TaskLogger

ACTION = None


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
    with pytest.raises(ValueError):
        tl.save(params)


def test_delta_minutes():
    mock_config = {}
    mock_session = None

    start = datetime.strptime("2023-11-29 19:30:48", "%Y-%m-%d %H:%M:%S")
    end = datetime.strptime("2023-11-30 19:30:48", "%Y-%m-%d %H:%M:%S")

    tl = TaskLogger(mock_config, mock_session)
    actual = tl.delta_minutes(start, end)
    assert actual == 24 * 60
