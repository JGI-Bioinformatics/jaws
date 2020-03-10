import click.testing
import pytest
import requests
from .conftest import get_queue, post_history

from jaws_client.analysis import run


def test_cli_queue(mock_user, monkeypatch, configuration):
    monkeypatch.setattr(requests, 'get', get_queue)

    runner = click.testing.CliRunner()
    result = runner.invoke(run, ['queue'])
    assert result.exit_code == 0


def test_cli_history(mock_user, monkeypatch, configuration):
    monkeypatch.setattr(requests, "post", post_history)
    runner = click.testing.CliRunner()
    result = runner.invoke(run, ["history"])
    assert result.exit_code == 0


@pytest.mark.xfail
def test_cli_status():
    pass


@pytest.mark.xfail
def test_cli_tasks():
    pass


@pytest.mark.xfail
def test_cli_metadata():
    pass


@pytest.mark.xfail
def test_cli_log():
    pass


@pytest.mark.xfail
def test_cli_cancel():
    pass


@pytest.mark.xfail
def test_cli_delete():
    pass


@pytest.mark.xfail
def test_cli_submit():
    pass
