#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
#
# RabbitMQ connection and communication wrapper classes.
#
import pika
import amqpstorm
import time

from jaws_jtm.common import logger
from jaws_jtm.lib.msgcompress import zdumps


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
        self.__connection = amqpstorm.Connection(
            RMQ_HOST,
            RMQ_USER,
            RMQ_PASS,
            port=RMQ_PORT,
            virtual_host=RMQ_VHOST,
            heartbeat=120,
            timeout=180,
        )

    def open(self):
        return self.__connection

    def close(self):
        self.__connection.close()

    def __exit__(self):
        self.__connection.close()


class JtmAmqpstormBase(object):
    def __init__(self, config, ppid=None, max_retries=None):
        self.config = config
        self.ppid = ppid
        self.max_retries = max_retries
        self.connection = None
        self.jgi_jtm_main_exch = config.configparser.get("JTM", "jgi_jtm_main_exch")
        self.jtm_task_request_q = config.configparser.get("JTM", "jtm_task_request_q")
        self.jtm_task_result_q = config.configparser.get("JTM", "jtm_task_result_q")
        self.jtm_status_request_q = config.configparser.get("JTM", "jtm_status_request_q")
        self.jtm_status_result_q = config.configparser.get("JTM", "jtm_status_result_q")
        self.jtm_inner_main_exch = config.configparser.get("JTM", "jtm_inner_main_exch")
        self.inner_result_queue_name = config.configparser.get(
            "JTM", "jtm_inner_result_q"
        )
        self.inner_request_queue_name = config.configparser.get(
            "JTM", "jtm_inner_request_q"
        )
        self.jtm_task_kill_exch = config.configparser.get("JTM", "jtm_task_kill_exch")
        self.jtm_task_kill_q = config.configparser.get("JTM", "jtm_task_kill_q")
        self.task_status = config.constants.TASK_STATUS
        self.done_flag = config.constants.DONE_FLAGS

    def create_connection(self):
        """Create a connection.
        :return:
        """
        attempts = 0
        while True:
            attempts += 1
            try:
                self.connection = RmqConnectionAmqpstorm(config=self.config).open()
                break
            except amqpstorm.AMQPError as why:
                logger.exception(why)
                if self.max_retries and attempts > self.max_retries:
                    break
                time.sleep(min(attempts * 2, 30))
            except KeyboardInterrupt:
                break

    def send_reply(self, message, mode, reply_msg):
        properties = {"correlation_id": message.correlation_id}
        try:
            response = amqpstorm.Message.create(
                message.channel, zdumps(str(reply_msg)), properties
            )
            response.publish(message.reply_to)
            message.ack()
        except amqpstorm.AMQPError as why:
            logger.exception(why)
            raise
        else:
            logger.info("Send reply back: {}, {}".format(mode, reply_msg))


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
            params = pika.ConnectionParameters(
                credentials=creds, host=RMQ_HOST, virtual_host=RMQ_VHOST, port=RMQ_PORT
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
