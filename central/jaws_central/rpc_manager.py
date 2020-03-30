import logging
from typing import List
from jaws_central import rpc_client, config


rpc = None


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class JawsRpcError(Exception):
    def __init__(self, message):
        super().__init__(message)


class JawsRpc(metaclass=Singleton):
    """Singleton which contains dictionary of site_id => rpc_client objects"""

    clients = {}

    def __init__(self, conf: config.JawsConfig):
        """Initialize an RPC client object for each configured Site.

        :param config: JawsConfig object
        :type: obj
        :return: JawsRpc object
        :rtype: obj
        """
        logger = logging.getLogger(__package__)
        for site_id in conf.sites.keys():
            logger.info(f"Initializing RPC client for {site_id}")
            self.clients[site_id] = rpc_client.RPC_Client(conf.sites[site_id])
        global rpc
        rpc = self

    def get_sites(self) -> List[str]:
        """Return list of JAWS-Site IDs.

        :return: list of JAWS-Site IDs (str)
        :rtype: list
        """
        return self.clients.keys()

    def get_client(self, site_id: str) -> rpc_client:
        """Get RPC client object of a Site.

        :param site_id: ID of the JAWS-Site
        :type site_id: str
        :return" rpc_client object
        :rtype: obj
        """
        site_id = site_id.upper()
        if site_id not in self.clients:
            raise JawsRpcError(f"Unknown Site, {site_id}")
        return self.clients[site_id]

    def is_valid_site(self, site_id: str) -> bool:
        """Determine if the user-specified site ID is valid.

        :param site_id: Unique idenfitier of the JAWS-Site
        :type site_id: str
        :return: True if site_id is valid, false otherwise
        :rtype: bool
        """
        return True if site_id.upper() in self.clients else False
