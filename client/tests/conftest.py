import pytest
import os
import jaws_client.user
import jaws_client.config
import json


@pytest.fixture
def configuration():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config = jaws_client.config.JawsConfig(f'{base_dir}/jaws_client.ini')
    return config


class MockUser:

    def __init__(self):
        pass

    def header(self):
        return {"Authorization": "Bearer ABCDEFGHIJKLMNOP"}


class MockResult:

    def __init__(self, response, status_code):
        self.response = response
        self.status_code = status_code

    def json(self):
        return json.dumps(self.response)


def get_queue(url, headers=None):
    return MockResult({}, 200)


def post_history(url, data=None, headers=None):
    return MockResult({}, 200)


@pytest.fixture
def mock_user(monkeypatch):
    monkeypatch.setattr(jaws_client.user, "User", MockUser)


@pytest.fixture
def wdl_path(tmp_path):
    wdls = tmp_path / "wdls"
    wdls.mkdir()
    return wdls


@pytest.fixture()
def input_file(wdl_path):
    wdl_dir = wdl_path
    inputs = wdl_dir / "inputs.json"
    path = wdl_dir.as_posix()
    contents = """
{
    "file1": "%s/test.wdl"
}"""
    inputs.write_text(contents % path)
    test_wdl = wdl_dir / "test.wdl"
    test_wdl.write_text("""
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
}
    """)
    return inputs


@pytest.fixture
def no_runtime_file(wdl_path):
    wdl_dir = wdl_path
    no_runtime_wdl = wdl_dir / "test.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
}"""
    no_runtime_wdl.write_text(contents)
    return no_runtime_wdl


@pytest.fixture
def good_runtime_file(wdl_path):
    wdl_dir = wdl_path
    good_runtime = wdl_dir / "test.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
  runtime {
    cpu: "1"
  }
}"""
    good_runtime.write_text(contents)
    return good_runtime


@pytest.fixture
def bad_runtime_file(wdl_path):
    wdl_dir = wdl_path
    bad_runtime = wdl_dir / "test.wdl"
    contents = """
task hello_world {

  String hi = "Hello world"

  command {
    echo ${hi}
  }
  runtime {
    cpus: "1"
  }
}"""
    bad_runtime.write_text(contents)
    return bad_runtime
