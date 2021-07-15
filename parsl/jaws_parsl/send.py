# run like:
# python send.py -c <# of CPUs> -m <memory> -s <site> -cmd $(pwd)/hello.sh

import click
import json
import logging
from jaws_parsl import config
from multiprocessing.connection import Client


@click.command()
@click.option("-c", "--cpus", help="Number of CPUs needed")
@click.option("-m", "--memory", help="Amount of memory needed")
@click.option("-s", "--site", help="Compute site")
@click.option("-cmd", "--command", help="Command script")
def send(cpus, memory, site, command):
    job_info = {
        "cpus": cpus,
        "memory": memory,
        "command": command,
        "site": site
    }

    mp_params = config.conf.get_mp_params()
    mp_password = mp_params["password"]
    mp_host = mp_params["host"]
    mp_port = mp_params["port"]

    conn = Client((mp_host, mp_port), authkey=bytes(mp_password, encoding='utf-8'))
    conn.send(json.dumps(job_info))
    conn.close()

    logging.debug(" [x] Sent job info.")
