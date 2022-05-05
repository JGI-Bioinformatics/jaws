import os
import shutil
import pytest
from dataclasses import dataclass
from jaws_site import perf_metrics_es, runs, config
from pprint import pprint
from deepdiff import DeepDiff


@pytest.mark.parametrize(
    "dir_name, exp_output",
    [
        ("/a/b/c/cromwell-executions/d/e/f", "cromwell-executions/d/e/f"),
        ("/cromwell-executions/a", "cromwell-executions/a"),
        ("a/not-a-cromwell-dir/b/c", "None")
    ]
)
def test_remove_beginning_path(dir_name, exp_output):
    obs_output = perf_metrics_es.remove_beginning_path(dir_name)
    assert obs_output == exp_output


@pytest.mark.parametrize(
    "working_dir, expect",
    [
        (
            "/var/udiMount/global/cscratch1/sd/jaws_jtm/\
jaws-prod/cromwell-executions/rna_count_wrapper/a304bf3c-f042-4af4-8cc5-20bfba245dc9\
/call-rna_count_DGE/DGE.rna_count_DGE/\
d04e8fb6-4a37-4f87-8e4e-8a7107523df9/call-diffGeneExp/execution/DGE_files",
            "a304bf3c-f042-4af4-8cc5-20bfba245dc9"
        ),
        (
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/\
cromwell-executions/jgi_dap_leo/f92649f2-31f9-4fc4-b1fd-18ed4de94c12\
/call-trimAlign_expt/shard-37/execution",
            "f92649f2-31f9-4fc4-b1fd-18ed4de94c12"
        ),
        (
            "cromwell-executions/\
                fq_count/8b1669c3-9bca-4047-99b4-e3556b71c4f4/\
                    call-count_seqs/execution",
            "8b1669c3-9bca-4047-99b4-e3556b71c4f4"
        ),
        (
            "/global/cfs/cdirs/jaws/jaws-install/jaws-supervisord-prod",
            "naw"
        ),
        (
            "None",
            "naw"
        ),
        (  # Hopefully tests that lru_cache is working
            "/global/cfs/cdirs/jaws/jaws-install/jaws-supervisord-prod",
            "naw"
        ),
        (
            "None",
            "naw"
        )
    ],
)
def test_extract_cromwell_es(working_dir, expect):
    out = perf_metrics_es.extract_jaws_info(working_dir=working_dir)
    assert out == expect


@pytest.fixture
def mock_runs(monkeypatch):
    class RunDbError(Exception):
        pass

    class MockRuns:
        def get_run_id_from_cromwell_id(*args):
            if data_obj.raises:
                raise RunDbError()
            else:
                return data_obj.run_id

    @dataclass
    class Data():
        raises = False
        run_id = None

    def mock_runs(*args, **kwargs):
        return MockRuns

    monkeypatch.setattr(runs, 'Run', mock_runs)

    data_obj = Data()
    return data_obj


@pytest.fixture
def mock_csv_file(tmp_path):
    csv_file = tmp_path / "metrics.csv"
    content = """
@timestamp,mem_rss,mem_vms,num_fds,current_dir
04-25-2022 01:35:03.801668,1,10,100,/a/cromwell-executions/aa/aaa
04-26-2022 01:35:03.801668,2,20,200,/b/cromwell-executions/bb/bbb
04-27-2022 01:35:03.801668,3,30,300,/c/cromwell-executions/cc/ccc
"""
    csv_file.write_text(content)
    return csv_file.as_posix()


@pytest.mark.parametrize(
    "cromwell_id, exp_run_id, raise_exception",
    [
        (1, 10, False),
        (2, "None", True)
    ]

)
def test_get_run_id(cromwell_id, exp_run_id, raise_exception, mock_runs):
    mock_runs.run_id = exp_run_id
    obs_run_id = perf_metrics_es.Metrics(None, None).get_run_id(cromwell_id)
    assert obs_run_id == exp_run_id


def test_process_csv(monkeypatch, mock_csv_file, mock_runs):
    exp_output = [
        {
            '@timestamp': '2022-04-25T01:35:03.801668',
            'cromwell_id': 'aaa',
            'current_dir': 'cromwell-executions/aa/aaa',
            'jaws_run_id': 1,
            'mem_rss': 1,
            'mem_total': 11,
            'mem_vms': 10,
            'num_fds': 100
        },
        {
            '@timestamp': '2022-04-26T01:35:03.801668',
            'cromwell_id': 'bbb',
            'current_dir': 'cromwell-executions/bb/bbb',
            'jaws_run_id': 1,
            'mem_rss': 2,
            'mem_total': 22,
            'mem_vms': 20,
            'num_fds': 200
        },
        {
            '@timestamp': '2022-04-27T01:35:03.801668',
            'cromwell_id': 'ccc',
            'current_dir': 'cromwell-executions/cc/ccc',
            'jaws_run_id': 1,
            'mem_rss': 3,
            'mem_total': 33,
            'mem_vms': 30,
            'num_fds': 300
        }
    ]

    mock_runs.run_id = 1
    obs_output = perf_metrics_es.Metrics(None, None).process_csv(mock_csv_file)
    assert obs_output == exp_output


def test_process_metrics(monkeypatch, mock_csv_file, mock_runs, config_file):

    from jaws_site import runs_es
    obs_output = []

    def mock_add_taskname_mapping(self, cromwell_id):
        mapping = {
            "cromwell-executions/aa/aaa": "task_a",
            "cromwell-executions/bb/bbb": "task_b",
            "cromwell-executions/cc/ccc": "task_c",
        }
        self.tasks[cromwell_id] = mapping

    def mock_get_run_id(_, cromwell_id):
        run_ids = {
            "aaa": 1,
            "bbb": 2,
            "ccc": 3,
        }
        return run_ids.get(cromwell_id)

    def mock_send_rpc_run_metadata(_, doc):
        obs_output.append(doc)

    monkeypatch.setattr(perf_metrics_es.Metrics, "add_taskname_mapping", mock_add_taskname_mapping)
    monkeypatch.setattr(perf_metrics_es.Metrics, "get_run_id", mock_get_run_id)
    monkeypatch.setattr(runs_es, "send_rpc_run_metadata", mock_send_rpc_run_metadata)

    done_dir = config.conf.get("PERFORMANCE_METRICS", "done_dir")
    proc_dir = config.conf.get("PERFORMANCE_METRICS", "processed_dir")
    proc_file = os.path.join(proc_dir, os.path.basename(mock_csv_file))
    done_file = os.path.join(done_dir, os.path.basename(mock_csv_file))
    exp_output = [
        {
            '@timestamp': '2022-04-25T01:35:03.801668',
            'cromwell_id': 'aaa',
            'current_dir': 'cromwell-executions/aa/aaa',
            'jaws_run_id': 1,
            'mem_rss': 1,
            'mem_total': 11,
            'mem_vms': 10,
            'num_fds': 100,
            'task_name': 'task_a'
        },
        {
            '@timestamp': '2022-04-26T01:35:03.801668',
            'cromwell_id': 'bbb',
            'current_dir': 'cromwell-executions/bb/bbb',
            'jaws_run_id': 2,
            'mem_rss': 2,
            'mem_total': 22,
            'mem_vms': 20,
            'num_fds': 200,
            'task_name': 'task_b'
        },
        {
            '@timestamp': '2022-04-27T01:35:03.801668',
            'cromwell_id': 'ccc',
            'current_dir': 'cromwell-executions/cc/ccc',
            'jaws_run_id': 3,
            'mem_rss': 3,
            'mem_total': 33,
            'mem_vms': 30,
            'num_fds': 300,
            'task_name': 'task_c'
        }
    ]

    assert done_dir
    assert proc_dir

    if not done_dir or not proc_dir:
        return

    os.makedirs(done_dir, exist_ok=True)
    os.makedirs(proc_dir, exist_ok=True)

    shutil.copyfile(mock_csv_file, done_file)

    perf_metrics_es.Metrics(None, None).process_metrics()

    assert bool(DeepDiff(obs_output, exp_output, ignore_order=True)) is False
    assert os.path.isfile(proc_file)

    shutil.rmtree(done_dir)
    shutil.rmtree(proc_dir)
