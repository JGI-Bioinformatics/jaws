#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
import sys

from jaws_jtm.config import pika, RMQ_HOST, RMQ_USER, RMQ_PASS, RMQ_VHOST, RMQ_PORT, JGI_JTM_MAIN_EXCH
from jaws_jtm.common import zdumps, logger

"""
Created on Apr 11 2013

RabbitMQ connection and communication wrapper classes.

:author Seung-Jin Sul
"""


class RmqConnection(object):
    """
    Class for RabbitMQ connection management
    """

    __connection = None

    def __init__(self):  # return a pika connection
        if not self.__connection:
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
        self.__connection.close()

    def is_open(self):
        return self.__connection.is_open()


class RmqConnectionHB(object):
    """
    Class for RMQ connection with hb option
    # With heartbeat_interval=0
    # ref) http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq,
    # http://stackoverflow.com/questions/34721178/pika-blockingconnection-rabbitmq-connection-closed
    """

    __connection = None

    def __init__(self):
        if not self.__connection:
            creds = pika.PlainCredentials(RMQ_USER, RMQ_PASS)
            params = pika.ConnectionParameters(
                credentials=creds,
                host=RMQ_HOST,
                virtual_host=RMQ_VHOST,
                port=RMQ_PORT,
                heartbeat=0,
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
    exch_name = exch if exch is not None else ""
    routing_key = queue if queue is not None else props.reply_to

    try:

        ch.basic_publish(
            exchange=JGI_JTM_MAIN_EXCH,
            routing_key=props.reply_to,  # use the queue which the client created
            properties=pika.BasicProperties(
                delivery_mode=delivery_mode,  # make message persistent
                correlation_id=props.correlation_id,
            ),
            body=zdumps(str(msg)),
        )
        ch.basic_ack(delivery_tag=method.delivery_tag)
    except Exception as detail:
        logger.exception(
            "Exception: Failed to send a msg to {} with {}".format(
                exch_name, routing_key
            )
        )
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)
