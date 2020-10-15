# run like:
# python send.py -cl cori -cmd $(pwd)/hello.sh

import click
import pika
import json
import logging
from jaws_parsl import config

@click.command()
@click.option("-cl", "--cluster", help="Compute cluster name")
@click.option("-cmd", "--command", help="Command script")
def send(cluster, command):
    job_info = {
        "cluster": cluster,
        "command": command
    }

    rmq_params = config.conf.get_rmq_params()

    rmq_user = rmq_params["user"]
    rmq_password = rmq_params["password"]
    rmq_host = rmq_params["host"]
    rmq_vhost = rmq_params["vhost"]
    rmq_port = rmq_params["port"]

    creds = pika.PlainCredentials(rmq_user, rmq_password)
    params = pika.ConnectionParameters(credentials=creds,
                                       host=rmq_host,
                                       virtual_host=rmq_vhost,
                                       port=rmq_port)
    connection = pika.BlockingConnection(params)
    channel = connection.channel()
    channel.queue_declare(queue='hello')
    channel.basic_publish(exchange='',
                          routing_key='hello',
                          body=json.dumps(job_info))
    logging.debug(" [x] Sent job info to RMQ server.")

