# run like:
# python send.py -cl cori -cmd $(pwd)/hello.sh

import click
import pika
import json

#from jaws_rpc import rpc_client

@click.command()
@click.option("-cl", "--cluster", help="Compute cluster name")
@click.option("-cmd", "--command", help="Command script")
def send(cluster, command):
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

#def send_update_task_status_msg(task_id: int, status_from, status_to: int, fail_code=None):
#    """
#    Publish a message for pushing task status change to JAWS Site
#
#    :param task_id:
#    :param status_from: None or status code 0 ~ 4 or -2 if failed
#    :param status_to: status code 0 ~ 4 or -2 if failed
#    :param fail_code: fail code if failed -1 ~ -7
#    :return: None
#    """
#    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
#    run_id = extract_cromwell_run_id(task_id)
#    reversed_task_status = dict(map(reversed, CONFIG.constants.TASK_STATUS.items()))
#    reversed_done_flags = dict(map(reversed, CONFIG.constants.DONE_FLAGS.items()))
#    data = {
#        "cromwell_run_id": run_id,  # this is not the JAWS run_id
#        "cromwell_job_id": task_id,
#        "status_from": reversed_task_status[status_from] if status_from is not None else "",
#        "status_to": reversed_task_status[status_to] if status_to is not None else "",
#        "timestamp": now,
#        "reason": reversed_done_flags[fail_code] if fail_code else None,
#    }
#
#    # send message to Site
#    try:
#        with rpc_client.RPC_Client({"host": CONFIG.configparser.get("SITE_RPC_CLIENT", "host"),
#                                    "vhost": CONFIG.configparser.get("SITE_RPC_CLIENT", "vhost"),
#                                    "port": CONFIG.configparser.get("SITE_RPC_CLIENT", "port"),
#                                    "user": CONFIG.configparser.get("SITE_RPC_CLIENT", "user"),
#                                    "queue": CONFIG.configparser.get("SITE_RPC_CLIENT", "queue"),
#                                    "password": CONFIG.configparser.get("SITE_RPC_CLIENT", "password")}
#                                   ) as rpc_cl:
#            wait_count = 0
#            response = rpc_cl.request("update_job_status", data)
#            logger.debug(f"Return msg from JAWS Site: {response}")
#            while "error" in response and response["error"]["message"] == "Server timeout":
#                wait_count += 1
#                if wait_count == 60:  # try for 1min
#                    logger.error("RPC reply timeout!")
#                    break
#                logger.debug(f"RPC reply delay. Wait for a result from JAWS Site RPC server: {response}")
#                time.sleep(1.0)
#                response = rpc_cl.request("update_job_status", data)
#    except Exception as error:
#        logger.error(f"RPC call failed: {error}")
#        raise
#
#    if "result" in response:
#        logger.debug(f"Status change message sent successfully: {data}")
#        pass
#    else:
#        logger.error(f"Status update failed: {response['error']['message']}")
#

