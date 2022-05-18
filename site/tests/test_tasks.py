from jaws_site import tasks
from jaws_site import cromwell
from jaws_site.tasks import TaskLog
from deepdiff import DeepDiff
from tests.conftest import mock_task_status_table, mock_task_summary_table, this_date


def test_task_status_table(monkeypatch):
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
    task_status = tasks.task_status_table()
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


def test_task_summary_table(monkeypatch):
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
        return self._task_log

    monkeypatch.setattr(TaskLog, "task_log", mock_task_log)

    def mock_get_cromwell_run_id(self):
        self._cromwell_run_id = "AAAA"
        return self._cromwell_run_id

    monkeypatch.setattr(TaskLog, "_get_cromwell_run_id", mock_get_cromwell_run_id)

    def mock_cromwell_job_summary(self):
        self._cromwell_job_summary = {"8919": ["fq_count.count_seqs", "00:10:00"]}
        return self._cromwell_job_summary

    monkeypatch.setattr(TaskLog, "cromwell_job_summary", mock_cromwell_job_summary)

    expected = [
        [
            "fq_count.count_seqs",
            "8919",
            False,
            "success",
            "2021-12-07 20:39:09",
            "0:00:07",
            "0:00:00",
            "00:10:00",
        ]
    ]

    example_run_id = 1
    mock_session = None
    tasks = TaskLog(mock_session, run_id=example_run_id)
    actual = tasks.task_summary_table()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_get_task_cromwell_dir_mapping(monkeypatch):
    class TaskMetadata:
        def __init__(self, task_data):
            self.tasks = task_data

        def summary(self):
            return self.tasks

    def mock_cromwell(*args, **kwargs):
        class CromMetadata:
            task1 = [
                [
                    "align.stats",
                    1,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-stats/execution",
                ],
            ]
            task2 = [
                [
                    "align.shard_wf:shard_wf.indexing",
                    2,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-indexing/execution",
                ],
                [
                    "align.shard_wf:shard_wf.map[0]",
                    3,
                    False,
                    "01:00:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-map/shard-0/execution",
                ],
                [
                    "align.shard_wf:shard_wf.merge",
                    4,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-merge/execution",
                ],
                [
                    "align.shard_wf:shard_wf.shard",
                    5,
                    False,
                    "00:30:00",
                    "/scratch/cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-shard/execution",
                ],
            ]
            task3 = [
                [
                    "no_cromwell_dir",
                    1,
                    False,
                    "00:30:00",
                    "",
                ],
            ]

            def __init__(self):
                self.tasks = {
                    "align.stats": TaskMetadata(CromMetadata.task1),
                    "align.bbmap_shard_wf": TaskMetadata(CromMetadata.task2),
                    "no_cromwell_dir": TaskMetadata(CromMetadata.task3),
                }

        return CromMetadata()

    monkeypatch.setattr(cromwell.Cromwell, "get_metadata", mock_cromwell)

    expected = {
        "cromwell-executions/align/C1/call-stats/execution": "align.stats",
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-indexing/execution": "align.shard_wf:shard_wf.indexing",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-map/shard-0/execution": "align.shard_wf:shard_wf.map[0]",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-merge/execution": "align.shard_wf:shard_wf.merge",  # noqa
        "cromwell-executions/align/C1/call-shard_wf/align.shard_wf/C2/call-shard/execution": "align.shard_wf:shard_wf.shard",  # noqa
    }
    mock_session = None
    tasks = TaskLog(mock_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    actual = tasks.get_task_cromwell_dir_mapping()
    assert bool(DeepDiff(actual, expected, ignore_order=False)) is False


def test_task_summary(mock_db_session, monkeypatch):
    monkeypatch.setattr(TaskLog, "task_summary_table", mock_task_summary_table)

    exp_results = {
        "task_abcd": {
            "cached": False,
            "cromwell_job_id": "cromwell_abcd",
            "max_time": "03:00:00",
            "queue_wait": "01:00:00",
            "queued": this_date,
            "result": "success",
            "run_time": "02:00:00",
        },
        "task_efgh": {
            "cached": True,
            "cromwell_job_id": "cromwell_efgh",
            "max_time": "06:00:00",
            "queue_wait": "04:00:00",
            "queued": this_date,
            "result": "success",
            "run_time": "05:00:00",
        },
    }

    tasks = TaskLog(mock_db_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    obs_results = tasks.task_summary()
    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False


def test_task_status(mock_db_session, monkeypatch):
    monkeypatch.setattr(TaskLog, "task_status_table", mock_task_status_table)

    exp_results = {
        "task_abcd": {
            "cromwell_job_id": "cromwell_abcd",
            "reason": "reason_abcd",
            "status": "success",
            "timestamp": this_date,
            "cached": False,
        },
        "task_efgh": {
            "cromwell_job_id": "cromwell_efgh",
            "reason": "reason_efgh",
            "status": "success",
            "timestamp": this_date,
            "cached": True,
        },
    }

    tasks = TaskLog(mock_db_session, cromwell_run_id="EXAMPLE-CROMWELL-RUN-ID")
    obs_results = tasks.task_status()
    assert bool(DeepDiff(obs_results, exp_results, ignore_order=True)) is False
