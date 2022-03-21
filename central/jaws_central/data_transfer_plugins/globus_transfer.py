import logging
import os
import globus_sdk
import jaws_central.config
from ..datatransfer_protocol import (
    DataTransferError,
    DataTransferAPIError,
    DataTransferNetworkError,
)


logger = logging.getLogger(__package__)


class DataTransfer:
    """
    Connects to Globus with a TransferClient using the client app authentication model.
    This uses a client id and a client secret to connect to a shared endpoint. Client ID
    must be given access to the shared endpoint in order for transfers to successfully complete
    """

    def __init__(self):
        self.globus_config = jaws_central.config.Configuration()
        self.endpoint_id = self.globus_config.get("GLOBUS", "endpoint_id")
        self.host_path = self.globus_config.get("GLOBUS", "host_path")
        self._transfer_client = None

    def transfer_client(self):
        """
        Return an authorized transfer client object.  Initialize if not exists.

        :return: datatransfer_client
        :rtype: globus_sdk.ConfidentialAppAuthClient
        """
        if not defined(self.transfer_client):
            try:
                client = globus_sdk.ConfidentialAppAuthClient(
                    self.globus_config.get("GLOBUS", "client_id"),
                    self.globus_config.get("GLOBUS", "client_secret"),
                )
                scopes = "urn:globus:auth:scope:transfer.api.globus.org:all"
                authorizer = globus_sdk.ClientCredentialsAuthorizer(client, scopes)
                transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
            except globus_sdk.GlobusAPIError as error:
                raise
            else:
                self._transfer_client = transfer_client
        return self._transfer_client

    def transfer_status(self, transfer_task_id: str) -> int:
        """
        Query Globus transfer service for transfer task status.
        Uploads are done only by JAWS's globus credentials (i.e. confidential app).

        :param transfer_task_id: Globus transfer task id
        :ptype: str
        :return: integer value representing the status of the transfer
        :rtype: DataTransfer.Status
        """
        try:
            transfer_client = self.transfer_client()
            task = transfer_client.get_task(transfer_task_id)
        except globus_sdk.GlobusError as error:
            raise DataTransferError(error)

        status = None
        globus_status = task["status"]
        if globus_status == "FAILED":
            status = "failed"
        elif globus_status == "INACTIVE":
            status = "inactive"
        else:
            status = "succeeded"
        return status

    @staticmethod
    def virtual_transfer_path(full_path: str, host_path: str) -> str:
        """Return an absolute path used by Globus transfer service that uses the host_path as root.

        :param full_path: The complete absolute path
        :type full_path: str
        :param host_path: Host path is root path of the globus endpoint
        :type host_path: str
        :return: virtual absolute path
        :rtype: str
        """
        if not full_path.startswith(host_path):
            raise ValueError(
                f"Path, {full_path}, is not accessible via Globus endpoint"
            )
        return os.path.join("/", os.path.relpath(full_path, host_path))

    def submit_transfer(self, manifest: list, **kwargs) -> str:
        """
        Submit a data transfer

        :param kwargs: src_endpoint_id, dest_endpoint_id, label, etc.
        :type kwargs: dict
        :param manifest: list of files/folders to transfer
        :type manifest: list
        :return: transfer task id
        :rtype: str
        """
        # required parameters
        if "src_endpoint_id" not in kwargs:
            raise KeyError("src_endpoint_id is not defined.")
        elif "dest_endpoint_id" not in kwargs:
            raise KeyError("dest_endpoint_id is not defined.")

        # optional parameters
        if "src_host_path" not in kwargs:
            kwargs["src_host_path"] = "/"
        if "dest_host_path" not in kwargs:
            kwargs["dest_host_path"] = "/"
        if "label" not in kwargs:
            kwargs["label"] = f"Xfer {len(manifest)} files/folders"

        try:
            transfer_task_id = self._submit_transfer(manifest, kwargs)
        except globus_sdk.GlobusAPIError as error:
            raise DataTransferAPIError(error)
        except globus_sdk.NetworkError as error:
            raise DataTransferNetworkError(error)
        except globus_sdk.GlobusError as error:
            raise DataTransferError(error)

        return transfer_task_id

    def _submit_transfer(self, manifest: list, **kwargs: dict):
        """
        Submit a transfer to Globus

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param dest_endpoint: destination endpoint UUID
        :param manifest: manifest of all the files to be transferred
        :return:
        """
        logger.debug(f"Globus xfer")

        transfer_client = self.transfer_client()

        tdata = globus_sdk.TransferData(
            transfer_client,
            src_endpoint,
            dest_endpoint,
            label=label,
            sync_level="mtime",
            verify_checksum=False,
            preserve_timestamp=False,
            notify_on_succeeded=False,
            notify_on_failed=False,
            notify_on_inactive=False,
            skip_activation_check=False,
            skip_source_errors=True,
        )

        for line in manifest:
            line = line.decode("UTF-8")
            source_path, dest_path, inode_type = line.split("\t")
            logger.debug(f"add transfer: {source_path} -> {dest_path}")
            virtual_src_path = self.virtual_transfer_path(
                source_path, kwargs["src_host_path"]
            )
            virtual_dest_path = self.virtual_transfer_path(
                dest_path, kwargs["dest_host_path"]
            )

            if inode_type == "D":
                # folders are always transferred recursively
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=True)
            else:
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=False)

        transfer_result = transfer_client.submit_transfer(tdata)
        return transfer_result["task_id"]

    def cancel_transfer(self, transfer_task_id: str) -> str:
        """Cancel a Globus transfer.

        :param transfer_task_id: Globus transfer task id
        :type transfer_task_id: str
        """
        try:
            transfer_client = self.transfer_client()
            transfer_response = transfer_client.cancel_task(transfer_task_id)
        except globus_sdk.GlobusAPIError as error:
            logger.error(
                f"Error cancelling Globus transfer, {transfer_task_id}: {error}"
            )
            return f"{error}"
        else:
            return transfer_response
