from .config import conf


def status():
    """Check all systems to see if they're available.

    :return: A dict of server:status
    :rtype: dict
    """
    result = {"JAWS-Central": "UP"}
    for site_id, rpc in conf.rpc_clients.items():
        response = rpc.request("server_status")
        if "error" not in response:
            result[site_id + "-Site"] = "UP"
            result[site_id + "-Cromwell"] = "UP"
        elif response["error"]["code"] == 500:
            result[site_id + "-Site"] = "DOWN"
            result[site_id + "-Cromwell"] = "Unknown"
        else:
            result[site_id + "-Site"] = "UP"
            result[site_id + "-Cromwell"] = "DOWN"
    return result, 200
