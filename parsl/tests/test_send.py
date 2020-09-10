import pytest
import jaws_parsl.send
from click.testing import CliRunner
from unittest.mock import Mock
from tests.__mocks__ import pika

@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "script.sh"
    content = """#!/bin/bash
    echo 'hello'
    """

    cfg.write_text(content)
    return cfg.as_posix()

def test_send(config_file, monkeypatch):
    mocked_pika = Mock()
    mocked_pika.BlockingConnection.return_value = pika.Connection()
    monkeypatch.setattr('jaws_parsl.send.pika', mocked_pika)
    runner = CliRunner()
    result = runner.invoke(jaws_parsl.send.send, ['-cl', 'cori', '-cmd', config_file])
    assert result.exit_code == 0

