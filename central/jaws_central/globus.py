import logging
import os
import globus_sdk
import jaws_central.config


logger = logging.getLogger(__package__)


class GlobusService:
    """
    Connects to Globus with a TransferClient using the client app authentication model.
    This uses a client id and a client secret to connect to a shared endpoint. Client ID
    must be given access to the shared endpoint in order for transfers to successfully complete
    """

    def __init__(self):
        self.globus_config = jaws_central.config.Configuration()

    def _authorize_transfer_client(self):
        """
        Helper method to create a Globus transfer client using
        client id and client secret for credentials.

        :return globus_sdk.TransferClient:
        """
        client = globus_sdk.ConfidentialAppAuthClient(self.globus_config.get("GLOBUS", "client_id"),
                                                      self.globus_config.get("GLOBUS", "client_secret"))
        scopes = "urn:globus:auth:scope:transfer.api.globus.org:all"
        authorizer = globus_sdk.ClientCredentialsAuthorizer(client, scopes)
        return globus_sdk.TransferClient(authorizer=authorizer)

    def _create_transfer_client(self):
        """
        Creates the transfer client using client id and client secret for
        credentials
        :return: globus_sdk.TransferClient
        """
        transfer_client = self._authorize_transfer_client()
        return transfer_client

    def transfer_status(self, task_id):
        """
        Query Globus transfer service for transfer task status.
        Uploads are done only by JAWS's globus credentials (i.e. confidential app).

        :param task_id: Globus task id
        :return: str
        """
        transfer_client = self._create_transfer_client()
        task = transfer_client.get_task(task_id)
        globus_status = task["status"]
        return globus_status

    def virtual_transfer_path(self, full_path, host_path):
        """Return an absolute path used by Globus transfer service that uses the host_path as root.

        :param full_path: The complete absolute path
        :type full_path: str
        :param host_path: Host path is root path of the globus endpoint
        :type host_path: str
        :return: virtual absolute path
        :rtype: str
        """
        if not full_path.startswith(host_path):
            raise ValueError(f"Path, {full_path}, is not accessible via Globus endpoint")
        return os.path.join("/", os.path.relpath(full_path, host_path))

    def submit_transfer(self, label, host_paths, src_endpoint, dest_endpoint, manifest_file):
        """
        Submit a transfer to Globus

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param dest_endpoint: destination endpoint UUID
        :param manifest_file: manifest of all the files to be transferred
        :return:
        """
        logger.debug(f"Globus xfer {label}")

        src_host_path = host_paths["src"]
        dest_host_path = host_paths["dest"]

        transfer_client = self._create_transfer_client()

        tdata = globus_sdk.TransferData(
            transfer_client,
            src_endpoint,
            dest_endpoint,
            label=label,
            sync_level="mtime",
            verify_checksum=True,
            preserve_timestamp=True,
            notify_on_succeeded=False,
            notify_on_failed=True,
            notify_on_inactive=True,
            skip_activation_check=False,
        )

        for line in manifest_file:
            line = line.decode("UTF-8")
            source_path, dest_path, inode_type = line.split("\t")
            logger.debug(f"add transfer: {source_path} -> {dest_path}")
            virtual_src_path = self.virtual_transfer_path(source_path, src_host_path)
            virtual_dest_path = self.virtual_transfer_path(dest_path, dest_host_path)

            if inode_type == "D":
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=True)
            else:
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=False)

        transfer_result = transfer_client.submit_transfer(tdata)
        return transfer_result["task_id"]
