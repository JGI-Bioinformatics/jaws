#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
#
# RabbitMQ connection and communication wrapper classes.
#
import pika
import sys
from jaws_jtm.common import logger
from jaws_jtm.lib.msgcompress import zdumps
from amqpstorm import Connection


class RmqConnectionAmqpstorm(object):
    __connection = None

    def __init__(self, config):
        assert config
        self.config = config
        RMQ_HOST = self.config.configparser.get("RMQ", "host")
        RMQ_USER = self.config.configparser.get("RMQ", "user")
        RMQ_PORT = self.config.configparser.getint("RMQ", "port")
        RMQ_PASS = self.config.configparser.get("RMQ", "password")
        RMQ_VHOST = self.config.configparser.get("RMQ", "vhost")
        self.__connection = Connection(RMQ_HOST, RMQ_USER, RMQ_PASS,
                                       port=RMQ_PORT,
                                       virtual_host=RMQ_VHOST,
                                       heartbeat=120,
                                       timeout=180,)

    def open(self):
        return self.__connection

    def close(self):
        self.__connection.close()


class RmqConnectionHB(object):
    """
    Class for RMQ connection with hb option
    # With heartbeat_interval=0
    # ref) http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq,
    # http://stackoverflow.com/questions/34721178/pika-blockingconnection-rabbitmq-connection-closed
    #
    # WIth heartbeat_interval=5
    # ref) https://github.com/pika/pika/blob/master/examples/basic_consumer_threaded.py
    # for functool and threading method to prevent connection lost
    """
    __connection = None

    def __init__(self, config=None):
        if not self.__connection:
            if config:
                self.config = config
                RMQ_USER = self.config.configparser.get("RMQ", "user")
                RMQ_HOST = self.config.configparser.get("RMQ", "host")
                RMQ_PORT = self.config.configparser.getint("RMQ", "port")
                RMQ_PASS = self.config.configparser.get("RMQ", "password")
                RMQ_VHOST = self.config.configparser.get("RMQ", "vhost")

            creds = pika.PlainCredentials(RMQ_USER, RMQ_PASS)
            params = pika.ConnectionParameters(credentials=creds,
                                               host=RMQ_HOST,
                                               virtual_host=RMQ_VHOST,
                                               port=RMQ_PORT,
                                               heartbeat=5,  # for functools method
                                               )
            self.__connection = pika.BlockingConnection(params)

        else:
            print("Already connected.")

    def open(self):
        return self.__connection

    def close(self):
        if self.__connection:
            self.__connection.close()

    def is_open(self):
        return self.__connection.is_open()


def rmq_close(rmq_conn_obj):
    rmq_conn_obj.close()


def send_msg_callback(ch, method, props, msg, exch=None, queue=None, delivery_mode=1):
    exch_name = exch if exch is not None else ''
    routing_key = queue if queue is not None else props.reply_to
    assert routing_key

    try:

        ch.basic_publish(exchange=exch_name,
                         routing_key=props.reply_to,  # use the queue which the client created
                         properties=pika.BasicProperties(delivery_mode=delivery_mode,  # make message persistent
                                                         correlation_id=props.correlation_id),
                         body=zdumps(str(msg)))
        ch.basic_ack(delivery_tag=method.delivery_tag)
    except Exception as detail:
        logger.exception("Exception: Failed to send a msg to {} with {}".format(exch_name, routing_key))
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)
