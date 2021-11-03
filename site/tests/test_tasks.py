from jaws_site import tasks
from jaws_site.tasks import TaskLog
from deepdiff import DeepDiff


def test_save_job_log(monkeypatch):
    example_log = [
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
            test_cromwell_job_id,
            "created",
            "ready",
            "2021-03-17 12:00:00",
            None,
        ],
        [
            test_cromwell_job_id,
            "ready",
            "queued",
            "2021-03-17 12:11:11",
            None,
        ],
        [
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
                "main.ex_task_1",
                "2222",
                "created",
                "ready",
                "2020-03-22 12:47:10",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                "ready",
                "queued",
                "2020-03-22 12:47:20",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                "queued",
                "pending",
                "2020-03-22 12:47:25",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                "pending",
                "running",
                "2020-03-22 12:48:01",
                None,
            ],
            [
                "main.ex_task_2",
                "2223",
                "created",
                "ready",
                "2020-03-22 12:48:11",
                None,
            ],
            [
                "main.ex_task_2",
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
            "main.ex_task_1",
            "2222",
            "running",
            "2020-03-22 12:48:01",
            None,
        ],
        [
            "main.ex_task_2",
            "2223",
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
                    "main.ex_task_1",
                    "2222",
                    "ready",
                    "2020-03-22 12:47:10",
                    None,
                ],
            ]
        elif run_id == 103:
            example_status = [
                [
                    "main.ex_task_2",
                    "2222",
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ]
            ]
        elif run_id == 104:
            example_status = [
                [
                    "main.ex_task_2",
                    "2222",
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ],
                [
                    "main.ex_task_3",
                    "2222",
                    "pending",
                    "2020-03-22 12:47:25",
                    None,
                ],
            ]
        elif run_id == 105:
            example_status = [
                [
                    "main.ex_task_4",
                    "2222",
                    "running",
                    "2020-03-22 12:48:01",
                    None,
                ],
                [
                    "main.ex_task_5",
                    "2223",
                    "ready",
                    "2020-03-22 12:48:11",
                    None,
                ],
            ]
        elif run_id == 106:
            example_status = [
                [
                    "main.ex_task_6",
                    "2223",
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
        104: "queued",
        105: "running",
        106: "running",
    }

    monkeypatch.setattr(TaskLog, "get_task_status", mock_get_task_status)
    mock_session = None

    for run_id, expected in run_id_and_expected.items():
        assert tasks.get_run_status(mock_session, run_id) == expected


def test_get_job_metadata(monkeypatch):
    example_task_summary = [
        ["main_workflow.goodbye", "12129", False, "0:00:44.787000"],
        ["main_workflow.hello", "12130", False, "0:00:45.297000"],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "12134",
            False,
            "0:00:40.851000",
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "12133",
            False,
            "0:00:40.203000",
        ],
        [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "12131",
            False,
            "0:00:40.282000",
        ],
        [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "12132",
            False,
            "0:00:41.661000",
        ],
    ]

    expected_task_info = {
        "12129": ["main_workflow.goodbye", "0:00:44.787000"],
        "12130": ["main_workflow.hello", "0:00:45.297000"],
        "12134": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "0:00:40.851000",
        ],
        "12133": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "0:00:40.203000",
        ],
        "12131": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "0:00:40.282000",
        ],
        "12132": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "0:00:41.661000",
        ],
    }

    mock_session = None
    tasks = TaskLog(mock_session)
    task_info = tasks.get_task_info(example_task_summary)
    assert bool(DeepDiff(task_info, expected_task_info, ignore_order=True)) is False


def test_get_task_log(monkeypatch):

    example_run_id = 9

    def mock_get_cromwell_run_id(self, run_id):
        assert run_id == example_run_id
        return "AAAA"

    monkeypatch.setattr(TaskLog, "_get_cromwell_run_id", mock_get_cromwell_run_id)

    def mock_get_task_summary(self, cromwell_run_id):
        example_task_summary = [
            ["main_workflow.goodbye", "5480", False, "0:00:44.787000"],
            ["main_workflow.hello", "5481", False, None],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
                "5482",
                False,
                None,
            ],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "5484",
                False,
                None,
            ],
            ["main_workflow.hello_world", None, True, None],
        ]
        return example_task_summary

    monkeypatch.setattr(TaskLog, "_get_task_summary", mock_get_task_summary)

    def mock_get_job_logs(self, cromwell_run_id):
        example_logs = {
            "5480": [
                ["created", "ready", "2021-04-15 12:42:08", None],
                ["ready", "queued", "2021-04-15 12:42:28", None],
                ["queued", "pending", "2021-04-15 12:42:29", None],
                ["pending", "running", "2021-04-15 12:43:44", None],
                ["running", "succeeded", "2021-04-15 12:45:01", None],
            ],
            "5481": [["created", "queued", "2021-04-15 12:49:95", None]],
            "5482": [["created", "queued", "2021-04-15 01:01:04", None]],
            "5484": [["created", "queued", "2021-04-15 01:11:52", None]],
        }
        return example_logs

    monkeypatch.setattr(TaskLog, "get_job_logs", mock_get_job_logs)

    expected = [
        ["main_workflow.hello_world", None, None, None, None, "Cached call"],
        [
            "main_workflow.goodbye",
            "5480",
            "created",
            "ready",
            "2021-04-15 12:42:08",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            "ready",
            "queued",
            "2021-04-15 12:42:28",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            "queued",
            "pending",
            "2021-04-15 12:42:29",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            "pending",
            "running",
            "2021-04-15 12:43:44",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            "running",
            "succeeded",
            "2021-04-15 12:45:01",
            "run_time=0:00:44.787000",
        ],
        [
            "main_workflow.hello",
            "5481",
            "created",
            "queued",
            "2021-04-15 12:49:95",
            None,
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "5482",
            "created",
            "queued",
            "2021-04-15 01:01:04",
            None,
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "5484",
            "created",
            "queued",
            "2021-04-15 01:11:52",
            None,
        ],
    ]

    mock_session = None
    tasks = TaskLog(mock_session)
    actual = tasks.get_task_log(example_run_id)
    print(actual)  # DEBUG
    print(expected)  # DEBUG
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False
