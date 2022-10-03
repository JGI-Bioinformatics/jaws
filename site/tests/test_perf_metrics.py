import pytest
import os
from jaws_site import perf_metrics
from jaws_site.perf_metrics import PerformanceMetrics
from tests.conftest import (
    MockSession,
    MockRunModel,
)
from jaws_site import runs
from jaws_site.runs import RunNotFoundError


tests_dir = os.path.dirname(os.path.abspath(__file__))


class MockRpcClient:
    def __init__(self, params=None, logger=None):
        pass

    def request(self, method, params={}):
        response = {"result": None}
        return response


@pytest.fixture
def mock_rpc_client(run):
    return MockRpcClient()


def test__init(mock_sqlalchemy_session, monkeypatch):
    a_rpc_client = MockRpcClient()
    PerformanceMetrics(mock_sqlalchemy_session, a_rpc_client)


def test_get_run_id(monkeypatch, mock_sqlalchemy_session):
    def mock_from_cromwell_run_id(session, cromwell_run_id):
        session = MockSession()
        data = MockRunModel(cromwell_run_id="xxx-xxx")
        run = runs.Run(session, data)
        return run

    monkeypatch.setattr(
        perf_metrics.runs.Run, "from_cromwell_run_id", mock_from_cromwell_run_id
    )

    a_rpc_client = MockRpcClient()
    pm = PerformanceMetrics(mock_sqlalchemy_session, a_rpc_client)
    ret = pm.get_run_id("xxx-xxx")
    assert ret == "99"

    # Test Exception
    def mock_from_cromwell_run_id(session, cromwell_run_id):
        raise RunNotFoundError

    monkeypatch.setattr(
        perf_metrics.runs.Run, "from_cromwell_run_id", mock_from_cromwell_run_id
    )

    a_rpc_client = MockRpcClient()
    pm = PerformanceMetrics(mock_sqlalchemy_session, a_rpc_client)
    ret = pm.get_run_id("xxx-xxx")
    assert ret == 0


def test_extract_jaws_info():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "cda3cb3f-535c-400d-ab61-2e41aeb35a80",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "e7f02164-2d3d-4cfb-828a-f3da23c43280",
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.extract_jaws_info(task_dir)
        assert result == expected


def test_remove_beginning_path():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.remove_beginning_path(task_dir)
        assert result == expected


def test_parse_perf_metrics_task_dir():
    test_data = [
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/jgi_dap_leo/cda3cb3f-535c-400d-ab61-2e41aeb35a80/call-trimAlign_expt/shard-9/execution",  # noqa
            "jgi_dap_leo.trimAlign_expt[9]",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-dev/cromwell-executions/main_workflow/e7f02164-2d3d-4cfb-828a-f3da23c43280/call-hello_and_goodbye_1/sub.hello_and_goodbye/3327f701-769a-49fe-b407-eb4be3a4a373/call-hello/execution",  # noqa
            "main_workflow.hello_and_goodbye_1:hello_and_goodbye.hello",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-bbmap_indexing/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.bbmap_indexing",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-shard/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.shard",
        ],
        [
            "cromwell-executions/main_align_wdl/8d357504-8eda-4d1c-b3c4-f87a0bdd29f8/call-bbmap_shard_wf/align.bbmap_shard_wf/27536cf7-4b8f-486c-9085-872625afaef1/call-alignment/shard-0/execution",  # noqa
            "main_align_wdl.bbmap_shard_wf:bbmap_shard_wf.alignment[0]",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/nmdc_metag/9969e560-7f3e-4305-b66e-324d199d3b33/call-annotation/awf.annotation/8ea94de6-e56f-4b4a-9015-25eb998e68cc/call-f_annotate/shard-0/fa.f_annotate/0e16887c-11ec-48d3-a343-27a9876a7c70/call-smart",  # noqa
            "nmdc_metag.annotation:annotation.f_annotate:f_annotate.smart",
        ],
        [
            "s3://jaws-site-prod/cromwell-execution/jgi_meta/bbfad8f4-f5de-43c2-94ef-4bd43f1de4d3/call-bbcms",  # noqa
            "jgi_meta.bbcms",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/inhomo/8cc6f043-13b1-4e18-b685-b9533e6704cf/call-processTaxon/shard-5",  # noqa
            "inhomo.processTaxon[5]",
        ],
        [
            "/tahoma/mscjgi/scratch/jaws-prod/cromwell-executions/nmdc_metag/c16846ff-3a7b-444e-a26b-ce484eb205b5/call-annotation/awf.annotation/e0910a3c-6ba1-43e3-8b4b-d275fb0601fb/call-s_annotate/shard-0/sa.s_annotate/1671df94-89d9-4418-a949-737038f458a0/call-fasta_merge",  # noqa
            "nmdc_metag.annotation:annotation.s_annotate:s_annotate.fasta_merge",
        ],
        [
            "/global/cscratch1/sd/jaws_jtm/jaws-prod/cromwell-executions/bbmap_shard_wf/be518d2a-6232-4f50-b6cf-7e1a3a995ad3/call-alignment/shard-0",  # noqa
            "bbmap_shard_wf.alignment[0]",
        ],
    ]

    for task_dir, expected in test_data:
        result = perf_metrics.parse_cromwell_task_dir_name(task_dir)
        assert result == expected


def test_compute_rate():
    import pandas as pd
    import numpy as np
    from pandas import Timestamp

    df_time = pd.DataFrame(
        {"B": [0, 1, 2, np.nan, 4]},
        index=[
            pd.Timestamp("20130101 09:00:00"),
            pd.Timestamp("20130101 09:00:02"),
            pd.Timestamp("20130101 09:00:03"),
            pd.Timestamp("20130101 09:00:05"),
            pd.Timestamp("20130101 09:00:06"),
        ],
    )
    result = perf_metrics.compute_rates(df_time)
    """expected result
    2013-01-01 09:00:00   0.0
    2013-01-01 09:00:02  10.0
    2013-01-01 09:00:03  20.0
    2013-01-01 09:00:05   NaN
    2013-01-01 09:00:06  40.0"""
    assert result.to_dict()["B"][Timestamp("2013-01-01 09:00:00")] == 0.0
    assert result.to_dict()["B"][Timestamp("2013-01-01 09:00:06")] == 40.0
