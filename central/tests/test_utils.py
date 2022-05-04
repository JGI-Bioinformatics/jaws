import collections
from jaws_rpc import rpc_index
from jaws_central import utils
import jaws_central.config
import pytest
import logging


@pytest.fixture()
def configuration(config_file):
    return jaws_central.config.Configuration(config_file)


def test_status():
    logger = logging.getLogger(__package__)
    rpci1 = rpc_index.RpcIndex({}, logger)
    rpci2 = rpc_index.rpc_index
    assert rpci1 == rpci2
    sites = rpci1.get_sites()
    assert isinstance(sites, collections.abc.KeysView)
    result, code = utils.status()
    assert result["JAWS-Central"] == "UP"


def test_info(configuration):
    for key in ("name", "version", "docs_url"):
        value = configuration.get("JAWS", key)
        assert value is not None
