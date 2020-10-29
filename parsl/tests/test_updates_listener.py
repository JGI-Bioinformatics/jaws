import jaws_parsl.updates_listener
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
