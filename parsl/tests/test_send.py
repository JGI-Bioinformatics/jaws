import jaws_parsl.send
from click.testing import CliRunner
from unittest.mock import Mock
from tests.__mocks__ import pika


def test_send(monkeypatch, script):
    mocked_pika = Mock()
    mocked_pika.BlockingConnection.return_value = pika.Connection()
    monkeypatch.setattr('jaws_parsl.send.pika', mocked_pika)
    runner = CliRunner()
    result = runner.invoke(jaws_parsl.send.send, ['-c', 4, '-m', '2GB', '-cmd', script])
    assert result.exit_code == 0
