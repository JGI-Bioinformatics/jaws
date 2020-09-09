import pytest
import jaws_parsl.send
from click.testing import CliRunner

@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "script.sh"
    content = """#!/bin/bash
    echo 'hello'
    """

    cfg.write_text(content)
    return cfg.as_posix()


def test_send(config_file):
    runner = CliRunner()
    result = runner.invoke(jaws_parsl.send.send, ['-cl', 'cori', '-cmd', config_file])
    assert result.exit_code == 1 # should fail w/connection error right now

