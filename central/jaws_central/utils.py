from jaws_central import rpc_manager, config
import amqpstorm
import urllib.parse


def status() -> dict:
    """Check all systems to see if they're available.

    :return: A dict of server:status
    :rtype: dict
    """
    result = {"JAWS-Central": "UP"}
    rpcm = rpc_manager.manager
    for site_id in rpcm.get_sites():
        client = rpcm.get_client(site_id)
        response = client.request("server_status")
        if "error" not in response:
            result[site_id + "-Site"] = "UP"
            result[site_id + "-Cromwell"] = "UP"
        elif response["error"]["code"] == 500:
            result[site_id + "-Site"] = "DOWN"
            result[site_id + "-Cromwell"] = "Unknown"
        else:
            result[site_id + "-Site"] = "UP"
            result[site_id + "-Cromwell"] = "DOWN"
        result[f"{site_id}-RMQ"] = _rmq_server_status(config.conf.sites[site_id])
    return result, 200


def _rmq_server_status(params):
    """Check RMQ server status.

    :param params: A dictionary containing host, user, password, vhost, queue
    :type params: dict
    :return: "UP" or "DOWN"
    :rtype: str
    """
    if params.get("amqp_vhost", None):
        uri = "amqp://%s:%s@%s:5672/%s?heartbeat=60" % (
            params["amqp_user"],
            urllib.parse.quote_plus(params["amqp_password"]),
            params["amqp_host"],
            params["amqp_vhost"],
        )
        try:
            connection = amqpstorm.UriConnection(uri)
            connection.check_for_errors()
        except amqpstorm.AMQPConnectionError:
            return "DOWN"
    else:
        try:
            connection = amqpstorm.Connection(
                params["amqp_host"], params["amqp_user"], params["amqp_password"]
            )
            connection.check_for_errors()
        except amqpstorm.AMQPConnectionError:
            return "DOWN"
    return "UP"
