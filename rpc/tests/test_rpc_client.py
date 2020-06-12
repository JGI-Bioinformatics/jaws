import jaws_rpc.rpc_client


invalid_test_config = {
    "host": "rmq.spaceforce.gov",
    "user": "blightyear",
    "vhost": "pplanet",
}


def test_init_invalid_config():
    try:
        rpc_client = jaws_rpc.rpc_client.RPC_Client(invalid_test_config)
        rpc_client = rpc_client  # so flake8 doesn't complain rpc_client is not used
    except Exception:
        pass
    else:
        assert -1
