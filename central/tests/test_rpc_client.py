import jaws_central.rpc_client
import json


def test_rpc_client_send_request(rpc_dict, mock_connection, mock_thread, mock_message):
    rpc_client = jaws_central.rpc_client.RPC_Client(params=rpc_dict)
    r = rpc_client.send_request("method", params=rpc_dict)
    assert r == "corr_id"


def test_rpc_client_get_response(rpc_dict, mock_connection, mock_thread, mock_message):
    rpc_client = jaws_central.rpc_client.RPC_Client(params=rpc_dict)
    rpc_client.queue = {
        "q1": json.dumps("response"),
        "q2": None,
        "q3": "bad response"
    }
    r1 = rpc_client.get_response("q1")
    r2 = rpc_client.get_response("q2")
    r3 = rpc_client.get_response("q3")
    assert r1 == "response"
    assert r2 is None
    assert r3["error"]
