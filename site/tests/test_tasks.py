from jaws_site import tasks
from jaws_site.tasks import TaskLog
from deepdiff import DeepDiff


def test_save_job_log(monkeypatch):
    example_log = [
        "AAAA-BBBB-CCCC",
        2345,
        "queued",
        "running",
        "2021-03-17 12:22:22",
        None,
    ]

    def mock_save_job_log(self, job_log):
        pass

    monkeypatch.setattr(TaskLog, "_save_job_log", mock_save_job_log)
    mock_session = None

    tasks = TaskLog(mock_session)
    tasks.save_job_log(*example_log)


def test_get_job_logs(monkeypatch):

    test_cromwell_run_id = "AAAA-BBBB-CCCC"
    test_cromwell_job_id = "2345"
    test_job_logs = [
        [
            test_cromwell_run_id,
            test_cromwell_job_id,
            "created",
            "ready",
            "2021-03-17 12:00:00",
            None,
        ],
        [
            test_cromwell_run_id,
            test_cromwell_job_id,
            "ready",
            "queued",
            "2021-03-17 12:11:11",
            None,
        ],
        [
            test_cromwell_run_id,
            test_cromwell_job_id,
            "queued",
            "running",
            "2021-03-17 12:22:22",
            None,
        ],
    ]

    def mock_get_job_logs(self, cromwell_run_ids):
        return test_job_logs

    monkeypatch.setattr(TaskLog, "_get_job_logs", mock_get_job_logs)
    mock_session = None

    tasks = TaskLog(mock_session)
    job_logs = tasks.get_job_logs([test_cromwell_run_id])
    assert test_cromwell_job_id in job_logs
    assert len(job_logs[test_cromwell_job_id]) == len(test_job_logs)
    assert job_logs[test_cromwell_job_id][0][0] == "created"


def test_task_status(monkeypatch):
    def mock_get_task_log(self, run_id):
        example_log = [
            [
                "EX_RUN_ID",
                "main.ex_task_1",
                1,
                "2222",
                "created",
                "ready",
                "2020-03-22 12:47:10",
                None,
            ],
            [
                "EX_RUN_ID",
                "main.ex_task_1",
                1,
                "2222",
                "ready",
                "queued",
                "2020-03-22 12:47:20",
                None,
            ],
            [
                "EX_RUN_ID",
                "main.ex_task_1",
                1,
                "2222",
                "queued",
                "pending",
                "2020-03-22 12:47:25",
                None,
            ],
            [
                "EX_RUN_ID",
                "main.ex_task_1",
                1,
                "2222",
                "pending",
                "running",
                "2020-03-22 12:48:01",
                None,
            ],
            [
                "EX_RUN_ID",
                "main.ex_task_2",
                1,
                "2223",
                "created",
                "ready",
                "2020-03-22 12:48:11",
                None,
            ],
            [
                "EX_RUN_ID",
                "main.ex_task_2",
                1,
                "2223",
                "ready",
                "queued",
                "2020-03-22 12:48:16",
                None,
            ],
        ]
        return example_log

    example_run_id = 1
    expected_status = [
        [
            "EX_RUN_ID",
            "main.ex_task_1",
            1,
            "2222",
            "pending",
            "running",
            "2020-03-22 12:48:01",
            None,
        ],
        [
            "EX_RUN_ID",
            "main.ex_task_2",
            1,
            "2223",
            "ready",
            "queued",
            "2020-03-22 12:48:16",
            None,
        ],
    ]

    monkeypatch.setattr(TaskLog, "get_task_log", mock_get_task_log)
    mock_session = None

    tasks = TaskLog(mock_session)

    task_status = tasks.get_task_status(example_run_id)
    assert bool(DeepDiff(task_status, expected_status, ignore_order=True)) is False


def test_get_run_status(monkeypatch):
    def mock_get_task_status(session, run_id):
        example_status = []
        if run_id == 102:
            example_status = [
                [
                    "EX_RUN_ID",
                    "main.ex_task_1",
                    1,
                    "2222",
                    "created",
                    "ready",
                    "2020-03-22 12:47:10",
                    None,
                ],
            ]
        elif run_id == 103:
            example_status = [
                [
                    "EX_RUN_ID",
                    "main.ex_task_2",
                    1,
                    "2222",
                    "ready",
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ]
            ]
        elif run_id == 104:
            example_status = [
                [
                    "EX_RUN_ID",
                    "main.ex_task_2",
                    1,
                    "2222",
                    "ready",
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ],
                [
                    "EX_RUN_ID",
                    "main.ex_task_3",
                    1,
                    "2222",
                    "queued",
                    "pending",
                    "2020-03-22 12:47:25",
                    None,
                ],
            ]
        elif run_id == 105:
            example_status = [
                [
                    "EX_RUN_ID",
                    "main.ex_task_4",
                    1,
                    "2222",
                    "pending",
                    "running",
                    "2020-03-22 12:48:01",
                    None,
                ],
                [
                    "EX_RUN_ID",
                    "main.ex_task_5",
                    1,
                    "2223",
                    "created",
                    "ready",
                    "2020-03-22 12:48:11",
                    None,
                ],
            ]
        elif run_id == 106:
            example_status = [
                [
                    "EX_RUN_ID",
                    "main.ex_task_6",
                    1,
                    "2223",
                    "running",
                    "success",
                    "2020-03-22 12:48:16",
                    None,
                ],
            ]
        return example_status

    run_id_and_expected = {
        101: None,
        102: "queued",
        103: "queued",
        104: "running",
        105: "running",
        106: "running",
    }

    monkeypatch.setattr(TaskLog, "get_task_status", mock_get_task_status)
    mock_session = None

    for run_id, expected in run_id_and_expected.items():
        assert tasks.get_run_status(mock_session, run_id) == expected
