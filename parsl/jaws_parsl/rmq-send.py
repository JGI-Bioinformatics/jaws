# run like:
# python rmq-send.py -cl cori -cmd $(pwd)/hello.sh

import click
import pika
import json

@click.command()
@click.option("-cl", "--cluster", help="Compute cluster name")
@click.option("-cmd", "--command", help="Command script")
def main(cluster, command):
    job_info = {
        "cluster": cluster,
        "command": command
    }
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue='hello')
    channel.basic_publish(exchange='',
                          routing_key='hello',
                          body=json.dumps(job_info))
    print(" [x] Sent job info to RMQ server.")

if __name__ == "__main__":
    main()

