# run like:
# python send.py -c <# of CPUs> -m <memory> -cmd $(pwd)/hello.sh

import click
import pika
import json
import logging
from jaws_parsl import config


@click.command()
@click.option("-c", "--cpus", help="Number of CPUs needed")
@click.option("-m", "--memory", help="Amount of memory needed")
@click.option("-cmd", "--command", help="Command script")
def send(cpus, memory, command):
    job_info = {
        "cpus": cpus,
        "memory": memory,
        "command": command
    }

    rmq_params = config.conf.get_rmq_params()

    rmq_user = rmq_params["user"]
    rmq_password = rmq_params["password"]
    rmq_host = rmq_params["host"]
    rmq_vhost = rmq_params["vhost"]
    rmq_port = rmq_params["port"]
    rmq_queue = rmq_params["queue"]
    rmq_exch = rmq_params["exchange"]

    creds = pika.PlainCredentials(rmq_user, rmq_password)
    params = pika.ConnectionParameters(credentials=creds,
                                       host=rmq_host,
                                       virtual_host=rmq_vhost,
                                       port=rmq_port)
    connection = pika.BlockingConnection(params)
    channel = connection.channel()
    channel.basic_publish(exchange=rmq_exch,
                          routing_key=rmq_queue,
                          body=json.dumps(job_info))
    logging.debug(" [x] Sent job info to RMQ server.")
