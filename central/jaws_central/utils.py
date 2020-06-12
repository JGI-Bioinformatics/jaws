from jaws_central import rpc_index, config
import amqpstorm


def status() -> dict:
    """Check all systems to see if they're available.

    :return: A dict of server:status
    :rtype: dict
    """
    result = {"JAWS-Central": "UP"}
    rpci = rpc_index.index
    for site_id in rpci.get_sites():
        client = rpci.get_client(site_id)
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
    if params.get("vhost", None):
        try:
            with amqpstorm.Connection(
                params["host"],
                params["user"],
                params["password"],
                int(params.get("port", 5672)),
                virtual_host=params["vhost"],
            ) as connection:
                connection.check_for_errors()
        except Exception:
            return "DOWN"
    else:
        try:
            connection = amqpstorm.Connection(
                params["host"], params["user"], params["password"]
            )
            connection.check_for_errors()
        except amqpstorm.AMQPConnectionError:
            return "DOWN"
    return "UP"
