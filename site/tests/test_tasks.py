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
    test_job_logs = [
        {
            "cromwell_job_id": test_cromwell_job_id,
            "status_from": "created",
            "status_to": "ready",
            "timestamp": "2021-03-17 12:00:00",
            "comment": None,
        },
        {
            "cromwell_job_id": test_cromwell_job_id,
            "status_from": "ready",
            "status_to": "queued",
            "timestamp": "2021-03-17 12:11:11",
            "comment": None,
        },
        {
            "cromwell_job_id": test_cromwell_job_id,
            "status_from": "queued",
            "status_to": "running",
            "timestamp": "2021-03-17 12:22:22",
            "comment": None,
        },
    ]

    def mock_get_job_logs(self):
        self._job_logs = test_job_logs

    monkeypatch.setattr(TaskLog, "_get_job_logs", mock_get_job_logs)
    mock_session = None

    tasks = TaskLog(mock_session, cromwell_run_id=test_cromwell_run_id)
    job_logs = tasks.job_logs()
    assert job_logs[0]["status_from"] == "created"


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
            {
                "name": "main_workflow.goodbye",
                "jobId": "12129",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello",
                "jobId": "12130",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
                "jobId": "12134",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "jobId": "12133",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
                "jobId": "12131",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
                "jobId": "12132",
                "cached": False,
                "maxTime": "0:10:00",
            },
        ]
        return example_task_summary

    monkeypatch.setattr(TaskLog, "cromwell_task_summary", mock_cromwell_task_summary)

    expected_task_info = {
        "12129": {"name": "main_workflow.goodbye", "max_time": "0:10:00"},
        "12130": {"name": "main_workflow.hello", "max_time": "0:10:00"},
        "12134": {
            "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "max_time": "0:10:00",
        },
        "12133": {
            "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "max_time": "0:10:00",
        },
        "12131": {
            "name": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.goodbye",
            "max_time": "0:10:00",
        },
        "12132": {
            "name": "main_workflow.hello_and_goodbye_2:hello_and_goodbye.hello",
            "max_time": "0:10:00",
        },
    }

    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    task_info = tasks.cromwell_job_summary()
    assert bool(DeepDiff(task_info, expected_task_info, ignore_order=True)) is False


def test_task_log(monkeypatch):

    example_run_id = 9

    def mock_get_cromwell_run_id(self):
        self._cromwell_run_id = "AAAA"

    monkeypatch.setattr(TaskLog, "_get_cromwell_run_id", mock_get_cromwell_run_id)

    def mock_cromwell_task_summary(self):
        self._cromwell_task_summary = [
            {
                "name": "main_workflow.goodbye",
                "jobId": "5480",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello",
                "jobId": "5481",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
                "jobId": "5482",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
                "jobId": "5484",
                "cached": False,
                "maxTime": "0:10:00",
            },
            {
                "name": "main_workflow.hello_world",
                "jobId": None,
                "cached": True,
                "maxTime": "0:10:00",
            },
        ]
        return self._cromwell_task_summary

    monkeypatch.setattr(TaskLog, "cromwell_task_summary", mock_cromwell_task_summary)

    def mock_job_logs(self):
        self._job_logs = {
            "5480": [
                {
                    "status_from": "created",
                    "status_to": "ready",
                    "timestamp": "2021-04-15 12:42:08",
                    "comment": None,
                },
                {
                    "status_from": "ready",
                    "status_to": "queued",
                    "timestamp": "2021-04-15 12:42:28",
                    "comment": None,
                },
                {
                    "status_from": "queued",
                    "status_to": "pending",
                    "timestamp": "2021-04-15 12:42:29",
                    "comment": None,
                },
                {
                    "status_from": "pending",
                    "status_to": "running",
                    "timestamp": "2021-04-15 12:43:44",
                    "comment": None,
                },
                {
                    "status_from": "running",
                    "status_to": "success",
                    "timestamp": "2021-04-15 12:45:01",
                    "comment": None,
                },
            ],
            "5481": [
                {
                    "status_from": "created",
                    "status_to": "queued",
                    "timestamp": "2021-04-15 12:49:95",
                    "comment": None,
                }
            ],
            "5482": [
                {
                    "status_from": "created",
                    "status_to": "queued",
                    "timestamp": "2021-04-15 01:01:04",
                    "comment": None,
                }
            ],
            "5484": [
                {
                    "status_from": "created",
                    "status_to": "queued",
                    "timestamp": "2021-04-15 01:11:52",
                    "comment": None,
                }
            ],
        }
        return self._job_logs

    monkeypatch.setattr(TaskLog, "job_logs", mock_job_logs)

    expected = [
        {
            "name": "main_workflow.hello_world",
            "cromwell_job_id": None,
            "cached": True,
            "status_from": None,
            "status_to": None,
            "timestamp": None,
            "comment": None,
        },
        {
            "name": "main_workflow.goodbye",
            "cromwell_job_id": "5480",
            "cached": False,
            "status_from": "created",
            "status_to": "ready",
            "timestamp": "2021-04-15 12:42:08",
            "comment": None,
        },
        {
            "name": "main_workflow.goodbye",
            "cromwell_job_id": "5480",
            "cached": False,
            "status_from": "ready",
            "status_to": "queued",
            "timestamp": "2021-04-15 12:42:28",
            "comment": None,
        },
        {
            "name": "main_workflow.goodbye",
            "cromwell_job_id": "5480",
            "cached": False,
            "status_from": "queued",
            "status_to": "pending",
            "timestamp": "2021-04-15 12:42:29",
            "comment": None,
        },
        {
            "name": "main_workflow.goodbye",
            "cromwell_job_id": "5480",
            "cached": False,
            "status_from": "pending",
            "status_to": "running",
            "timestamp": "2021-04-15 12:43:44",
            "comment": None,
        },
        {
            "name": "main_workflow.goodbye",
            "cromwell_job_id": "5480",
            "cached": False,
            "status_from": "running",
            "status_to": "success",
            "timestamp": "2021-04-15 12:45:01",
            "comment": None,
        },
        {
            "name": "main_workflow.hello",
            "cromwell_job_id": "5481",
            "cached": False,
            "status_from": "created",
            "status_to": "queued",
            "timestamp": "2021-04-15 12:49:95",
            "comment": None,
        },
        {
            "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.goodbye",
            "cromwell_job_id": "5482",
            "cached": False,
            "status_from": "created",
            "status_to": "queued",
            "timestamp": "2021-04-15 01:01:04",
            "comment": None,
        },
        {
            "name": "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
            "cromwell_job_id": "5484",
            "cached": False,
            "status_from": "created",
            "status_to": "queued",
            "timestamp": "2021-04-15 01:11:52",
            "comment": None,
        },
    ]

    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    actual = tasks.task_log()
    print(actual)  # DEBUG
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False
