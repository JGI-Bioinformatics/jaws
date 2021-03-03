import logging
import os
import globus_sdk
import jaws_site.config


logger = logging.getLogger(__package__)


class GlobusService:
    """
    Connects to Globus with a TransferClient using the client app authentication model.
    This uses a client id and a client secret to connect to a shared endpoint. Client ID
    must be given access to the shared endpoint in order for transfers to successfully complete
    """

    def __init__(self):
        self.globus_config = jaws_site.config.Configuration()
        self.endpoint_id = self.globus_config.get("GLOBUS", "endpoint_id")
        self.host_path = self.globus_config.get("GLOBUS", "host_path")

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

    def submit_transfer(self, label, dest_endpoint, src_dir, dest_dir):
        """
        Submit a transfer to Globus

        :param label: label for the data transfer
        :param dest_endpoint: destination endpoint UUID
        :param src_dir: source directory
        :param dest_dir: destination directory
        :return:
        """
        logger.debug(f"Globus xfer {label}")

        transfer_client = self._create_transfer_client()
        tdata = globus_sdk.TransferData(
            transfer_client,
            self.endpoint_id,
            dest_endpoint,
            label=label,
            sync_level="mtime",
            verify_checksum=False,
            preserve_timestamp=True,
            notify_on_succeeded=False,
            notify_on_failed=False,
            notify_on_inactive=False,
            skip_activation_check=True,
        )

        virtual_src_path = self.virtual_transfer_path(src_dir, self.host_path)
        tdata.add_item(virtual_src_path, dest_dir, recursive=True)
        transfer_result = transfer_client.submit_transfer(tdata)
        return transfer_result["task_id"]
