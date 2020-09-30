import argparse
import datetime
import logging
import time
import os
import json
from functools import partial

import pika
import parsl
from parsl import bash_app, AUTO_LOGNAME
from parsl.config import Config

logger = None

G_TASK_COUNTER = 0
G_TASK_TABLE = {}
G_UPDATES_CHANNEL = None


def on_task_callback(task_id, future):
    logger.debug(f"[Task:{task_id}] Received callback")
    global G_UPDATES_CHANNEL
    
    try:
        result = future.result()
    except Exception as e:    
        logger.exception(f"[Task:{task_id}] failed with exception : {e}")
        G_UPDATES_CHANNEL.push_task_status(task_id, 'FAILED')
    else:
        logger.info(f"[Task:{task_id}] completed successfully")
        G_UPDATES_CHANNEL.push_task_status(task_id, 'COMPLETED')
    

def on_message_callback(ch, method, properties, body):
    
    logger.debug("Received message script.")
    logger.debug(f" Body : {body}")
    global G_UPDATES_CHANNEL
    
    try:
        message = json.loads(body)
    except:
        logger.exception(f"Failed to decode message: {message}")
        return

    global G_TASK_TABLE

    executor = message['cluster']    
    future = run_script(message['command'])
    task_id = future.tid
    logger.debug(f"[Task:{task_id} Launched")
    G_UPDATES_CHANNEL.push_task_status(task_id, 'LAUNCHED')
    
    G_TASK_TABLE[task_id] = {'future' : future,
                             'received_at' : time.time(),
                             'completed_at' : None}
    
    future.add_done_callback(partial(on_task_callback, task_id))

    logger.debug("Done")


@bash_app(executors=['cori'])
def run_script(script, stdout=AUTO_LOGNAME, stderr=AUTO_LOGNAME):
    cmd = 'bash ' + script
    return cmd

class UpdatesChannel():

    def __init__(self, address, qname):
        self.address = address
        self.qname = qname
        self.connection = pika.BlockingConnection(pika.ConnectionParameters(address))
        self.channel = self.connection.channel()
        self.channel.queue_declare(queue=qname)
        logger.debug(f"Task updates on {address}:{qname}")

    def _publish(self, message):
        logger.debug(f"Sending message {message}")
        self.channel.basic_publish(exchange='',
                                   routing_key=self.qname,
                                   body=message)

    def push_task_status(self, task_id, status):
        status_message = {'task_id': task_id,
                          'timestamp': str(datetime.datetime.now()),
                          'status': status}
        self._publish(json.dumps(status_message))
        

class TasksChannel():

    def __init__(self, address, qname):

        self.address = address
        self.qname = qname
        connection = pika.BlockingConnection(pika.ConnectionParameters(address))
        channel = connection.channel()
        channel.queue_declare(queue=qname)
        channel.basic_consume(queue=qname,
                              auto_ack=True,
                              on_message_callback=on_message_callback)
        self.channel = channel
        
    def listen(self):
        logger.info(' [*] Waiting for messages. To exit press CTRL+C')
        try:
            self.channel.start_consuming()
        except Exception as e:
            logger.exception("Caught exception while waiting of RMQ")
            logger.info("Ignoring error and continuing")
            pass

def start_file_logger(filename, name='parsl-rabbitmq', level=logging.DEBUG, format_string=None):
    """Add a stream log handler.

    Args:
        - filename (string): Name of the file to write logs to
        - name (string): Logger name
        - level (logging.LEVEL): Set the logging level.
        - format_string (string): Set the format string

    Returns:
       -  None
    """
    if format_string is None:
        format_string = "%(asctime)s.%(msecs)03d %(name)s:%(lineno)d [%(levelname)s]  %(message)s"

    global logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def cli():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--port", default=8080,
                        help="Port at which the service will listen on")
    parser.add_argument("-a", "--address", default="localhost",
                        help="RabbitMQ address to connect to")
    parser.add_argument("-q", "--qname", default="hello",
                        help="RabbitMQ queue to listen on")
    parser.add_argument("--tasks_qname", default="task_updates",
                        help="RabbitMQ queue to publish task updates on")
    parser.add_argument("-c", "--config", default=None,
                        help="Config file")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")

    # from parsl.configs.htex_local import config
    from cori import config

    dfk = parsl.load(config)
    parsl_run_dir = dfk.run_dir

    args = parser.parse_args()
    
    if args.logfile:
        logfile_path = args.logfile
    else:
        logfile_path = f'{parsl_run_dir}/parsl-RabbitMQ.log'

    os.makedirs(os.path.dirname(logfile_path), exist_ok=True)
    
    start_file_logger(logfile_path, level=logging.DEBUG if args.debug else logging.INFO)
    logger.info("Starting")

    global G_UPDATES_CHANNEL
    G_UPDATES_CHANNEL = UpdatesChannel(args.address, args.tasks_qname)
    tasks_channel = TasksChannel(args.address, args.qname)
    tasks_channel.listen()

    logger.info("Exiting")


if __name__ == "__main__":
    cli()
