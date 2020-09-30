import argparse
import logging
import os
import pika
from jaws_rpc import rpc_client

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
        except Exception as e:
            logger.exception("Caught exception while waiting of RMQ")
            logger.info("Ignoring error and continuing")
            pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--address", default="localhost",
                        help="RabbitMQ address to connect to")
    parser.add_argument("--tasks_qname", default="task_updates",
                        help="RabbitMQ queue to publish task updates on")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")
    args = parser.parse_args()
    
    if args.logfile:
        logfile_path = args.logfile
    else:
        logfile_path = f'task_updates.log'

    start_file_logger(logfile_path, level=logging.DEBUG if args.debug else logging.INFO)
    logger.info("Starting")
    
    updates_channel = UpdatesChannel(args.address, args.tasks_qname)
    updates_channel.listen()

"""
# --------------------------------------------------------------------------------------------------
def send_update_task_status_msg(task_id: int, status_from, status_to: int, fail_code=None):
    
    # Publish a message for pushing task status change to JAWS Site

    # :param task_id:
    # :param status_from: None or status code 0 ~ 4 or -2 if failed
    # :param status_to: status code 0 ~ 4 or -2 if failed
    # :param fail_code: fail code if failed -1 ~ -7
    # :return: None
  
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    run_id = extract_cromwell_run_id(task_id)
    reversed_task_status = dict(map(reversed, CONFIG.constants.TASK_STATUS.items()))
    reversed_done_flags = dict(map(reversed, CONFIG.constants.DONE_FLAGS.items()))
    data = {
        "cromwell_run_id": run_id,  # this is not the JAWS run_id
        "cromwell_job_id": task_id,
        "status_from": reversed_task_status[status_from] if status_from is not None else "",
        "status_to": reversed_task_status[status_to] if status_to is not None else "",
        "timestamp": now,
        "reason": reversed_done_flags[fail_code] if fail_code else None,
    }

    # send message to Site
    try:
        with rpc_client.RPC_Client({"host": CONFIG.configparser.get("SITE_RPC_CLIENT", "host"),
                                    "vhost": CONFIG.configparser.get("SITE_RPC_CLIENT", "vhost"),
                                    "port": CONFIG.configparser.get("SITE_RPC_CLIENT", "port"),
                                    "user": CONFIG.configparser.get("SITE_RPC_CLIENT", "user"),
                                    "queue": CONFIG.configparser.get("SITE_RPC_CLIENT", "queue"),
                                    "password": CONFIG.configparser.get("SITE_RPC_CLIENT", "password")}
                                   ) as rpc_cl:
            wait_count = 0
            response = rpc_cl.request("update_job_status", data)
            logger.debug(f"Return msg from JAWS Site: {response}")
            while "error" in response and response["error"]["message"] == "Server timeout":
                wait_count += 1
                if wait_count == 60:  # try for 1min
                    logger.error("RPC reply timeout!")
                    break
                logger.debug(f"RPC reply delay. Wait for a result from JAWS Site RPC server: {response}")
                time.sleep(1.0)
                response = rpc_cl.request("update_job_status", data)
    except Exception as error:
        logger.error(f"RPC call failed: {error}")
        raise

    if "result" in response:
        logger.debug(f"Status change message sent successfully: {data}")
        pass
    else:
        logger.error(f"Status update failed: {response['error']['message']}")
"""
 
