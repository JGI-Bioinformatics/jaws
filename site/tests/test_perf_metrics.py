import os

import pytest
from jaws_site import perf_metrics, runs
from jaws_site.perf_metrics import PerformanceMetrics
from jaws_site.runs import RunNotFoundError

from tests.conftest import MockRunModel, MockSession

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


def test_compute_rate():
    import numpy as np
    import pandas as pd
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
