from tests.conftest import MockSession, MockLogger
from datetime import datetime
from jaws_site.tasks import TaskLog


def test_did_run_start(monkeypatch):
    # test 1: yes
    def mock__select_rows(self):
        self.data = [
            [
                "call-count_seqs",
                "queued",
                datetime.strptime("2023-04-03 17:20:35", "%Y-%m-%d %H:%M:%S"),
            ],
            [
                "call-count_seqs",
                "running",
                datetime.strptime("2023-04-03 17:20:52", "%Y-%m-%d %H:%M:%S"),
            ],
            [
                "call-count_seqs",
                "succeeded",
                datetime.strptime("2023-04-03 17:20:55", "%Y-%m-%d %H:%M:%S"),
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
