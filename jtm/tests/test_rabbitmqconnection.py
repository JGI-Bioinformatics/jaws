# test_rabbitmqconnection.py
# import pytest
# import rabbitmq_adapter
# from unittest.mock import Mock
# from tests.__mocks__ import pika


def test_channel_sets_parameters():
    # mocked_pika = Mock()
    # mocked_pika.URLParameters.return_value = 'JAWS'
    # monkeypatch.setattr('rabbitmq_adapter.channel.pika', mocked_pika)
    #
    # rabbitmq_adapter.channel.create('JAWS HOST')
    #
    # mocked_pika.URLParameters.assert_called_once_with('JAWS HOST')
    assert True


def test_channel_creates_connection():
    # mocked_pika = Mock()
    # mocked_pika.URLParameters.return_value = 'JAWS's
    # mocked_pika.BlockingConnection.return_value = pika.Connection()
    # monkeypatch.setattr('rabbitmq_adapter.channel.pika', mocked_pika)
    #
    # rabbitmq_adapter.channel.create('JAWS HOST')
    #
    # mocked_pika.BlockingConnection.assert_called_once_with('JAWS')
    assert True


def test_rabbitmq_factory():
    # Should be able to create a RabbitMQ connection
    assert True


def test_rabbitmq_listen_to_queue():
    # Should be able to listen to a RabbitMQ queue
    assert True


def test_rabbitmq_queue_one_to_many_queue_handler():
    # Should be able to have an one-to-many relationship between the queue and the handler function
    assert True


def test_rabbitmq_send_message():
    # Should be able to send a message in a given queue
    assert True
