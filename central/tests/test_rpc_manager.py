import jaws_central.config
import jaws_central.rpc_manager
from jaws_central.rpc_client import RPC_Client


def test_rpc_manager_get_sites(config_file, mock_connection):

    # call destructor to remove old references
    jaws_central.config.Configuration._destructor()

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)

    expected_sites = ["LBNL", "NERSC"]

    rpc_manager = jaws_central.rpc_manager.JawsRpc(conf=cfg)
    sites = rpc_manager.get_sites()

    assert sites == expected_sites


def test_rpc_manager_get_client(config_file, mock_connection):

    # call destructor to remove old references
    jaws_central.config.Configuration._destructor()

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)

    rpc_manager = jaws_central.rpc_manager.JawsRpc(conf=cfg)
    lbnl_client = rpc_manager.get_client("LBNL")
    nersc_client = rpc_manager.get_client("NERSC")

    assert isinstance(lbnl_client, RPC_Client)
    assert isinstance(nersc_client, RPC_Client)


def test_rpc_manager_valid_site(config_file, mock_connection):

    # call destructor to remove old references
    jaws_central.config.Configuration._destructor()

    config_path = config_file
    cfg = jaws_central.config.Configuration(config_path)

    rpc_manager = jaws_central.rpc_manager.JawsRpc(conf=cfg)

    assert rpc_manager.is_valid_site("LBNL")
    assert rpc_manager.is_valid_site("NERSC")
    assert not rpc_manager.is_valid_site("NASA")
