import jaws_parsl.updates_listener
from unittest.mock import Mock
from tests.__mocks__ import pika


def test_start_file_logger(tmp_path):
    tmp_log = tmp_path / "test.txt"
    jaws_parsl.updates_listener.start_file_logger(tmp_log)


def test_on_message_callback(caplog):
    ch = pika.Channel()
    method = pika.spec.Basic.Deliver()
    properties = pika.spec.BasicProperties()
    body = b'message body'
    jaws_parsl.updates_listener.on_message_callback(ch, method, properties, body)
    for record in caplog.records:
        assert record.levelname == "DEBUG"
    assert "Received message script." in caplog.text
    assert "message body" in caplog.text


def test_listen(monkeypatch, caplog):
    mocked_pika = Mock()
    mocked_pika.BlockingConnection.return_value = pika.Connection()
    monkeypatch.setattr('jaws_parsl.updates_listener.pika', mocked_pika)
    updates_channel = jaws_parsl.updates_listener.UpdatesChannel("localhost", "queue")
    updates_channel.listen()
    for record in caplog.records:
        assert record.levelname == "INFO"
    assert "Waiting for messages." in caplog.text
