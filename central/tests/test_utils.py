import collections
from jaws_rpc import rpc_index
from jaws_central import utils


def test_status():
    rpci1 = rpc_index.RPC_Index({})
    rpci2 = rpc_index.rpc_index
    assert(rpci1 == rpci2)
    sites = rpci1.get_sites()
    assert isinstance(sites, collections.abc.KeysView)
    result, code = utils.status()
    assert result["JAWS-Central"] == "UP"
