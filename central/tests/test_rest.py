import pytest
from datetime import datetime
import jaws_central.rest
import jaws_central.runs
import jaws_central.models
import jaws_central.config
from tests.conftest import MockSession, MockRunModel


class MockUser:
    @property
    def transfer_refresh_token(self):
        return "abcdefghijklmnopqrstuvwxyz"


class MockQuery:
    @staticmethod
    def get(user):
        return MockUser()


class MockResponse:
    def __init__(self, status_code):
        self.status_code = status_code


class MockGetInactiveUploadRun:
    def __init__(self, *args, **kwargs):
        self.cromwell_run_id = "90b333ed-8095-474e-9ced-df86f3c99241"
        self.download_id = None
        self.id = 123
        self.input_site_id = "NERSC"
        self.compute_site_id = "NERSC"
        self.status = "upload inactive"
        self.result = None
        self.submission_id = "2ba6eb76-22e0-4d49-8eeb-b0c6683dfa30"
        self.submitted = datetime.strptime("2020-05-14 23:08:50", "%Y-%m-%d %H:%M:%S")
        self.updated = datetime.strptime("2020-05-14 23:27:15", "%Y-%m-%d %H:%M:%S")
        self.upload_id = 50
        self.user_id = "dduck"
        self.tag = "example"
        self.wdl_file = "example.wdl"
        self.json_file = "example.json"
        self.manifest_json = "{}"


class MockRunWithId:
    def __init__(self, *args, **kwargs):
        self.id = 123


@pytest.fixture()
def mock_database(monkeypatch):
    monkeypatch.setattr(jaws_central.models.db, "session", MockSession)


def test_list_sites(configuration):
    user = "test_user"
    result, code = jaws_central.rest.list_sites(user)
    expected_sites = ["JGI", "NERSC", "AWS"]
    assert isinstance(result, list)
    for site in result:
        assert isinstance(site, dict)
        assert "site_id" in site
        assert site["site_id"] in expected_sites
        assert "max_ram_gb" in site


def test_run_metadata(monkeypatch):
    def mock_get_run(user_id, run_id):
        return MockRunModel()

    def mock_abort_if_pre_cromwell(run):
        return

    def mock_rpc_call(user_id, run, method, params={}):
        assert isinstance(user_id, str)
        assert isinstance(run, MockRunModel)
        assert method == "run_metadata"

    monkeypatch.setattr(jaws_central.rest, "_get_run", mock_get_run)
    monkeypatch.setattr(
        jaws_central.rest, "_abort_if_pre_cromwell", mock_abort_if_pre_cromwell
    )
    monkeypatch.setattr(jaws_central.rest, "rpc_call", mock_rpc_call)
    jaws_central.rest.run_metadata("test_user", 123)


def test_get_performance_metrics(monkeypatch):
    """
    test if the response is parsed correctly
    """

    def mock_get_run(user_id, run_id):
        return MockRunWithId()

    def mock_search_es(host, port, api_key, index, query, aggregations=None):
        assert isinstance(host, str)
        assert isinstance(port, str)
        assert isinstance(api_key, str)
        assert isinstance(index, str)
        assert isinstance(query, dict)

        response = {
            "took": 2,
            "_shards": {},
            "hits": {
                "max_score": 1.0,
                "hits": [
                    {
                        "_index": "dummyperfmetrics",
                        "_type": "_doc",
                        "_id": "47JAY38BP88SDbn3qnsw",
                        "_score": 1.0,
                        "_source": {
                            "read_chars": 2668002692,
                            "pid": 31390,
                            "write_chars": 410010798,
                            "cmdline": "",
                            "memory_percent": 0.026769763565699004,
                            "write_count": 4142,
                            "current_dir": "",
                            "cpu_percent": 0,
                            "jaws_run_id": 123,
                            "name": "jtm",
                            "num_fds": 17,
                            "mem_vms": 361418752,
                        },
                    }
                ],
            },
        }
        return response

    monkeypatch.setattr(jaws_central.rest, "_get_run", mock_get_run)

    monkeypatch.setattr(jaws_central.rest, "_search_elastic_search", mock_search_es)
    metrics = jaws_central.rest.get_performance_metrics("test_user", 123)
    assert len(metrics) == 1
    assert "name" in metrics[0]
    assert "cpu_percent" in metrics[0]
    assert type(metrics[0]["num_fds"]) == int
