from jaws_central import rpc_manager


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
    return result, 200
