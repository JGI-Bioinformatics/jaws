import jaws_parsl.send
from click.testing import CliRunner
from unittest.mock import Mock
from tests.__mocks__ import multiprocessing


def test_send(monkeypatch, script):
    mocked_mp = Mock()
    mocked_mp.Client.return_value = multiprocessing.Client()
    monkeypatch.setattr('jaws_parsl.send.Client', mocked_mp)
    runner = CliRunner()
    result = runner.invoke(jaws_parsl.send.send, ['-c', 4, '-m', '2G', '-s', 'cori', '-cmd', script])
    assert result.exit_code == 0
