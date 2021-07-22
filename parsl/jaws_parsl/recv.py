import argparse
import logging
import os
import re
import time
import datetime
import parsl
import jaws_parsl
from parsl import bash_app, AUTO_LOGNAME
from jaws_rpc.rpc_client import RpcClient
from multiprocessing.connection import Listener
from functools import partial

logger = None
rpc_params = {}
executor = ''
cpus = 0
mem = ''


@bash_app(executors=[executor])
def run_script(script, stdout=AUTO_LOGNAME, stderr=AUTO_LOGNAME):
    cmd = 'bash ' + script
    return cmd


def get_cromwell_run_id(msg_cmd):
    cromwell_id_regex = r"[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}"
    try:
        srch = re.search(cromwell_id_regex, msg_cmd)
        run_id = srch.group()
        logger.debug(f"[Cromwell Run ID: {run_id}")
    except Exception as e:
        logger.exception(f"Failed to get Cromwell run ID for task: {e}")
        raise
    return run_id


def start_file_logger(filename, name='jaws-parsl-recv.log', level=logging.DEBUG, format_string=None):
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


def on_task_callback(task_id, run_id, future):
    logger.debug(f"[Task:{task_id}] Received callback")

    try:
        result = future.result()
    except Exception as e:
        logger.exception(f"[Task:{task_id}] failed with exception : {e}")
        update_site('LAUNCHED', 'FAILED', task_id, run_id)
    else:
        logger.info(f"[Task:{task_id}] completed successfully")
        logger.info(f"Result: {result}")
        update_site('LAUNCHED', 'COMPLETED', task_id, run_id)


def update_site(status_from, status_to, task_id, run_id):
    global rpc_params
    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    data = {
        "cromwell_run_id": run_id,
        "cromwell_job_id": task_id,
        "status_from": status_from,
        "status_to": status_to,
        "timestamp": now,
    }

    # send message to Site
    try:
        with RpcClient({"user": rpc_params["user"],
                        "password": rpc_params["password"],
                        "host": rpc_params["host"],
                        "vhost": rpc_params["vhost"],
                        "port": rpc_params["port"],
                        "queue": rpc_params["queue"]},
                       logger
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


class TasksChannel():

    def __init__(self, host, port, password):
        self.listener = Listener((host, port), authkey=bytes(password, encoding='utf-8'))

    def listen(self):
        logger.info(' [*] Waiting for messages. To exit press CTRL+C')
        running = True
        while running:
            conn = self.listener.accept()
            msg = conn.recv()
            try:
                global cpus, mem, site, executor
                cpus = msg['cpus']
                mem = msg['memory']

                # assume memory is spec'd in GB via memory_gb in WDL
                # check memory request and route on Cori accordingly
                if site == "CORI":
                    if int(mem) <= 128:
                        executor = 'cori_genepool'
                    else:
                        executor = 'cori_exvivo'
                elif site == "JGI":
                    executor = 'lbl'
                future = run_script(msg)
                task_id = future.tid
                logger.debug(f"[Task:{task_id} Launched")
                run_id = get_cromwell_run_id(msg)
                future.add_done_callback(partial(on_task_callback, task_id, run_id))
                update_site('', 'LAUNCHED', task_id, run_id)
            except Exception as e:
                logger.exception(f"Caught exception while waiting for message: {e}")
                running = False
        self.listener.close()


def cli():
    global rpc_params
    rpc_params = jaws_parsl.config.conf.get_rpc_params()

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None,
                        help="Config file")
    parser.add_argument("-l", "--logfile", default=None,
                        help="Path to logfile")
    parser.add_argument("-d", "--debug", action='store_true',
                        help="Enables debug logging")

    site_id = rpc_params = jaws_parsl.config.conf.get_site_id()
    if site_id == "CORI":
        from execs import config_cori as config
    elif site_id == "LBL":
        from execs import config_lbl as config
    else:
        e = "Unknown site_id specified in Parsl backend config"
        raise ValueError(e)

    config.retries = 3
    dfk = parsl.load(config)
    parsl_run_dir = dfk.run_dir
    # parsl_run_id = dfk.run_id

    args = parser.parse_args()

    if args.logfile:
        logfile_path = args.logfile
    else:
        logfile_path = f'{parsl_run_dir}/jaws-parsl-recv.log'

    os.makedirs(os.path.dirname(logfile_path), exist_ok=True)

    start_file_logger(logfile_path, level=logging.DEBUG if args.debug else logging.INFO)
    logger.info("Starting")

    mp_params = jaws_parsl.config.conf.get_mp_params()
    mp_password = mp_params["password"]
    mp_host = mp_params["host"]
    mp_port = mp_params["port"]

    tasks_channel = TasksChannel(mp_host, mp_port, mp_password)
    tasks_channel.listen()

    logger.info("Exiting")


if __name__ == "__main__":
    cli()
