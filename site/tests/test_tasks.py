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

    tasks = TaskLog(mock_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    tasks.save_job_log(*example_log)


def test_job_logs(monkeypatch):

    test_cromwell_run_id = "AAAA-BBBB-CCCC"
    test_cromwell_job_id = "2345"

    def mock_select_job_logs(self):
        self._job_logs = [
            [
                test_cromwell_job_id,
                "queued",
                "running",
                "2021-03-17 12:22:22",
                None,
            ],
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
        ]
        return self._job_logs

    monkeypatch.setattr(TaskLog, "_select_job_logs", mock_select_job_logs)
    mock_session = None

    tasks = TaskLog(mock_session, cromwell_run_id=test_cromwell_run_id)
    job_logs = tasks.job_logs()
    assert job_logs[test_cromwell_job_id][0][0] == "created"


def test_task_status(monkeypatch):
    def mock_task_log(self):
        example_log = [
            [
                "main.ex_task_1",
                "2222",
                False,
                "created",
                "ready",
                "2020-03-22 12:47:10",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "ready",
                "queued",
                "2020-03-22 12:47:20",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "queued",
                "pending",
                "2020-03-22 12:47:25",
                None,
            ],
            [
                "main.ex_task_1",
                "2222",
                False,
                "pending",
                "running",
                "2020-03-22 12:48:01",
                None,
            ],
            [
                "main.ex_task_2",
                "2223",
                False,
                "created",
                "ready",
                "2020-03-22 12:48:11",
                None,
            ],
            [
                "main.ex_task_2",
                "2223",
                False,
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
            False,
            "running",
            "2020-03-22 12:48:01",
            None,
        ],
        [
            "main.ex_task_2",
            "2223",
            False,
            "queued",
            "2020-03-22 12:48:16",
            None,
        ],
    ]

    monkeypatch.setattr(TaskLog, "task_log", mock_task_log)
    mock_session = None

    tasks = TaskLog(mock_session, run_id=example_run_id)
    task_status = tasks.task_status()
    assert bool(DeepDiff(task_status, expected_status, ignore_order=True)) is False


def test_run_status(monkeypatch):
    def mock_task_status(session):
        example_status = []
        if run_id == 102:
            example_status = [
                [
                    "main.ex_task_1",
                    "2222",
                    False,
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
                    False,
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
                    False,
                    "queued",
                    "2020-03-22 12:47:20",
                    None,
                ],
                [
                    "main.ex_task_3",
                    "2222",
                    False,
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
                    False,
                    "running",
                    "2020-03-22 12:48:01",
                    None,
                ],
                [
                    "main.ex_task_5",
                    "2223",
                    False,
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
                    False,
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

    monkeypatch.setattr(TaskLog, "task_status", mock_task_status)
    mock_session = None

    for run_id, expected in run_id_and_expected.items():
        assert tasks.get_run_status(mock_session, run_id) == expected


def test_cromwell_job_summary(monkeypatch):
    example_run_id = 99

    def mock_cromwell_task_summary(self):
        example_task_summary = [
            ["main_workflow.goodbye", "12129", False, "0:10:00"],
            ["main_workflow.hello", "12130", False, "0:10:00"],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
                "12134",
                False,
                "0:10:00",
            ],
            [
                "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "12133",
                False,
                "0:10:00",
            ],
            [
                "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
                "12131",
                False,
                "0:10:00",
            ],
            [
                "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
                "12132",
                False,
                "0:10:00",
            ],
        ]
        return example_task_summary

    monkeypatch.setattr(TaskLog, "cromwell_task_summary", mock_cromwell_task_summary)

    expected_task_info = {
        "12129": ["main_workflow.goodbye", "0:10:00"],
        "12130": ["main_workflow.hello", "0:10:00"],
        "12134": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "0:10:00",
        ],
        "12133": [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "0:10:00",
        ],
        "12131": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "0:10:00",
        ],
        "12132": [
            "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "0:10:00",
        ],
    }

    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    task_info = tasks.cromwell_job_summary()
    assert bool(DeepDiff(task_info, expected_task_info, ignore_order=True)) is False


def test_task_log(monkeypatch):

    example_run_id = 9

    def mock_get_cromwell_run_id(self):
        self._cromwell_run_id = "AAAA"
        return self._cromwell_run_id

    monkeypatch.setattr(TaskLog, "_get_cromwell_run_id", mock_get_cromwell_run_id)

    def mock_cromwell_task_summary(self):
        example_task_summary = [
            ["main_workflow.goodbye", "5480", False, "0:10:00"],
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

    monkeypatch.setattr(TaskLog, "cromwell_task_summary", mock_cromwell_task_summary)

    def mock_job_logs(self):
        example_logs = {
            "5480": [
                ["created", "ready", "2021-04-15 12:42:08", None],
                ["ready", "queued", "2021-04-15 12:42:28", None],
                ["queued", "pending", "2021-04-15 12:42:29", None],
                ["pending", "running", "2021-04-15 12:43:44", None],
                ["running", "success", "2021-04-15 12:45:01", None],
            ],
            "5481": [["created", "queued", "2021-04-15 12:49:95", None]],
            "5482": [["created", "queued", "2021-04-15 01:01:04", None]],
            "5484": [["created", "queued", "2021-04-15 01:11:52", None]],
        }
        return example_logs

    monkeypatch.setattr(TaskLog, "job_logs", mock_job_logs)

    expected = [
        ["main_workflow.hello_world", None, True, None, None, None, None],
        [
            "main_workflow.goodbye",
            "5480",
            False,
            "created",
            "ready",
            "2021-04-15 12:42:08",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            False,
            "ready",
            "queued",
            "2021-04-15 12:42:28",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            False,
            "queued",
            "pending",
            "2021-04-15 12:42:29",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            False,
            "pending",
            "running",
            "2021-04-15 12:43:44",
            None,
        ],
        [
            "main_workflow.goodbye",
            "5480",
            False,
            "running",
            "success",
            "2021-04-15 12:45:01",
            None,
        ],
        [
            "main_workflow.hello",
            "5481",
            False,
            "created",
            "queued",
            "2021-04-15 12:49:95",
            None,
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "5482",
            False,
            "created",
            "queued",
            "2021-04-15 01:01:04",
            None,
        ],
        [
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "5484",
            False,
            "created",
            "queued",
            "2021-04-15 01:11:52",
            None,
        ],
    ]

    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    actual = tasks.task_log()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_task_summary(monkeypatch):
    def mock_task_log(self):
        self._task_log = [
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "created",
                "ready",
                "2021-12-07 20:39:09",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "ready",
                "queued",
                "2021-12-07 20:39:09",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "queued",
                "pending",
                "2021-12-07 20:39:09",
                "slurm_jid=45352308",
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "pending",
                "running",
                "2021-12-07 20:39:16",
                None,
            ],
            [
                "fq_count.count_seqs",
                "8919",
                False,
                "running",
                "success",
                "2021-12-07 20:39:16",
                None,
            ],
        ]

    monkeypatch.setattr(TaskLog, "task_log", mock_task_log)

    expected = [
        [
            "fq_count.count_seqs",
            "8919",
            False,
            "success",
            "2021-12-08 04:39:09",
            "0:00:07",
            "0:00:00",
            "00:10:00",
        ]
    ]

    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    actual = tasks.task_summary()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False
