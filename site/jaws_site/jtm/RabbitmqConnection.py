#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

from .Config import *
from .Common import *

"""
Created on Apr 11 2013

RabbitMQ connection and communication wrapper classes.

:author Seung-Jin Sul
"""

import atexit

class RmqConnection(object):
    """
    Class for RabbitMQ connection management
    """
    __connection = None

    def __init__(self): # return a pika connection
        atexit.register(rmq_close, rmqConnObj=self)
        if not self.__connection:
            creds = pika.PlainCredentials(RMQ_USER, RMQ_PASS)
            params = pika.ConnectionParameters(credentials=creds,
                                               host=RMQ_HOST,
                                               virtual_host=RMQ_VHOST,
                                               port=RMQ_PORT
                                               # heartbeat_interval=0
                                               # heartbeat_interval=0,  # https://github.com/pika/pika/issues/636, https://pika.readthedocs.io/en/latest/modules/parameters.html
                                               # connection_attempts=100,
                                               # socket_timeout=1
                                               )
            self.__connection = pika.BlockingConnection(params)

            ## Tip?
            ##credentials = pika.PlainCredentials(RABBITMQ_USER, RABBITMQ_PASS)
            ##connection = pika.BlockingConnection(pika.ConnectionParameters(
            ##    credentials=credentials,
            ##        host=RABBITMQ_HOST,
            ##         socket_timeout=300))
            ##channel = connection.channel()
            ##
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

    def __init__(self): # return a pika connection
        atexit.register(rmq_close, rmqConnObj=self)
        if not self.__connection:
            creds = pika.PlainCredentials(RMQ_USER, RMQ_PASS)
            params = pika.ConnectionParameters(credentials=creds,
                                               host=RMQ_HOST,
                                               virtual_host=RMQ_VHOST,
                                               # heartbeat_interval=hb,  # http://stackoverflow.com/questions/14572020/handling-long-running-tasks-in-pika-rabbitmq
                                               port=RMQ_PORT,
                                               heartbeat=600,
                                               blocked_connection_timeout=600,
                                               # connection_attempts=100,
                                               # socket_timeout=1
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


def rmq_close(rmqConnObj):
    # print "Running 'atexit()' handler"
    rmqConnObj.close()


def send_msg_callback(ch, method, props, msg, exch=None, queue=None, deliveryMode=1):
    exchName = exch if exch is not None else ''
    rKey = queue if queue is not None else props.reply_to
    # logger.debug("exch, props.reply_to, queue = {} {} {}".format(exch, props.reply_to, queue))

    try:
        # ch.basic_publish(exchange=exchName,
        #                  routing_key=rKey,
        #                  properties=pika.BasicProperties(
        #                      delivery_mode=deliveryMode,  # default =1 = non-persistent
        #                      correlation_id=props.correlation_id),
        #                  body=zdumps(str(msg)))
        # ch.queue_declare(queue=props.reply_to,
        #                  durable=True,
        #                  exclusive=False,
        #                  auto_delete=True)
        # ch.queue_bind(exchange=JGI_JTM_MAIN_EXCH,
        #                         queue=props.reply_to,
        #                         routing_key=props.reply_to)
        ch.basic_publish(exchange=JGI_JTM_MAIN_EXCH,
                         routing_key=props.reply_to,  # use the queue which the client created
                         properties=pika.BasicProperties(delivery_mode=deliveryMode,  # make message persistent
                                                         correlation_id=props.correlation_id),
                         body=zdumps(str(msg)))
        ch.basic_ack(delivery_tag=method.delivery_tag)
    except Exception as detail:
        logger.exception("Exception: Failed to send a msg to {} with {}".format(exchName, rKey))
        logger.exception("Detail: %s", str(detail))
        sys.exit(1)


# def send_msg_normal(exch, queue='', deliveryMode=1):
#     try:
#         while 1:
#             ch.basic_publish(exchange=exchName, routing_key='', body=message)
#             logger.debug("Send HB to workers, %r" % (message,))
#             time.sleep(hbInterval)
#     except Exception as e:
#         logger.critical("Something wrong in send_hb_to_workers_thread(): %s", e)
#         raise
