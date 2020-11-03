import argparse
import logging
import pika
from jaws_parsl import config


def start_file_logger(filename, name='task_updates', level=logging.DEBUG, format_string=None):
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


def on_message_callback(ch, method, properties, body):
    logger.debug("Received message script.")
    message = body.decode('utf-8')
    logger.debug(f" Body : {message}")
    print(f"[MESSAGE] {message}")


class UpdatesChannel():

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
        except Exception:
            logger.exception("Caught exception while waiting for RMQ")
            logger.info("Ignoring error and continuing")
            pass


if __name__ == '__main__':

    rpc_params = config.conf.get_rpc_params()

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--address", default=rpc_params["host"],
                        help="RabbitMQ address to connect to")
    parser.add_argument("--tasks_qname", default=rpc_params["queue"],
                        help="RabbitMQ queue to publish task updates on")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    args = parser.parse_args()

    if args.logfile:
        logfile_path = args.logfile
    else:
        logfile_path = 'task_updates.log'

    start_file_logger(logfile_path, level=logging.DEBUG if args.debug else logging.INFO)
    logger.info("Starting")

    updates_channel = UpdatesChannel(args.address, args.tasks_qname)
    updates_channel.listen()
