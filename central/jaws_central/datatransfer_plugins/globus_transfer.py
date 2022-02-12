import logging
import os
import globus_sdk
import jaws_central.config
from ..datatransfer_protocol import (
    SiteTransfer,
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

    def _authorize_transfer_client(self) -> globus_sdk.TransferClient:
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

    def _create_transfer_client(self) -> globus_sdk.TransferClient:
        """
        Creates the transfer client using client id and client secret for
        credentials
        :return: globus_sdk.TransferClient
        """
        transfer_client = self._authorize_transfer_client()
        return transfer_client

    def transfer_status(self, task_id: str) -> int:
        """
        Query Globus transfer service for transfer task status.
        Uploads are done only by JAWS's globus credentials (i.e. confidential app).

        :param task_id: Globus task id
        :return: integer value representing the different statuses. The statuses are defined in
         DataTransfer_protocol.Status class.
        :rtype: integer
        """
        try:
            transfer_client = self._create_transfer_client()
            task = transfer_client.get_task(task_id)
        except globus_sdk.GlobusError as error:
            raise DataTransferError(error)

        status = None
        globus_status = task["status"]
        if globus_status == "FAILED":
            status = SiteTransfer.status.failed
        elif globus_status == "INACTIVE":
            status = SiteTransfer.status.inactive
        else:
            status = SiteTransfer.status.succeeded
        return status

    def virtual_transfer_path(self, full_path: str, host_path: str) -> str:
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

    def submit_download(self, metadata: dict, src_dir: str, dest_dir: str) -> str:
        """
        Submit a data transfer

        :param metadata: dict containing entries for data transfer, i.e. label.
        :type metadata: dict
        :param src_dir: source directory
        :type src_dir: str
        :param dest_dir: destination directory
        :return:
        """
        if 'label' not in metadata:
            raise KeyError("Missing label key in input metadata.")
        label = metadata['label']

        logger.debug(f"Globus xfer {label}")

        if 'dest_endpoint' not in metadata:
            raise KeyError("Missing dest_endpoint key in input metadata.")
        if 'label' not in metadata:
            raise KeyError("Missing label key in input metadata.")

        try:
            transfer_client = self._create_transfer_client()
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.endpoint_id,
                metadata['dest_endpoint'],
                label=metadata['label'],
                sync_level="mtime",
                verify_checksum=False,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=False,
                notify_on_inactive=False,
                skip_activation_check=True,
                skip_source_errors=True,
            )
        except globus_sdk.GlobusError as error:
            raise DataTransferError(error)

        virtual_src_path = self.virtual_transfer_path(src_dir, self.host_path)
        tdata.add_item(virtual_src_path, dest_dir, recursive=True)
        transfer_result = transfer_client.submit_transfer(tdata)
        return transfer_result["task_id"]

    def submit_upload(self, metadata: dict, manifest_files: list) -> str:
        """
        Submit a data transfer

        :param metadata: dict containing entries for data transfer, i.e. label.
        :type metadata: dict
        :param manifest_files: list of files to transfer
        :type manifest_files: list
        :return: upload task id
        :rtype: str
        """
        logger.debug("GLOBUS: submit_upload ...")
        if 'host_paths' not in metadata:
            raise KeyError("host_paths is not defined in input metadata.")
        elif 'input_endpoint' not in metadata:
            raise KeyError("input_endpoint is not defined in input metadata.")
        elif 'compute_endpoint' not in metadata:
            raise KeyError("compute_endpoint is not defined in input metadata.")
        elif 'run_id' not in metadata:
            raise KeyError("run_id is not defined in input metadata.")

        try:
            upload_task_id = self._submit_upload(
                f"Upload run {metadata['run_id']}",
                metadata["host_paths"],
                metadata["input_endpoint"],
                metadata["compute_endpoint"],
                manifest_files,
            )
        except globus_sdk.GlobusAPIError as error:
            raise DataTransferAPIError(error)
        except globus_sdk.NetworkError as error:
            raise DataTransferNetworkError(error)
        except globus_sdk.GlobusError as error:
            raise DataTransferError(error)

        logger.debug(f"{upload_task_id=}")
        return upload_task_id

    def cancel_transfer(self, transfer_task_id: str) -> str:
        """Cancel a Globus transfer.

        :param transfer_task_id: Globus transfer task id
        :type transfer_task_id: str
        """
        try:
            transfer_client = self._authorize_transfer_client()
            transfer_response = transfer_client.cancel_task(transfer_task_id)
        except globus_sdk.GlobusAPIError as error:
            logger.error(f"Error cancelling Globus transfer, {transfer_task_id}: {error}")
            return f"{error}"
        else:
            return transfer_response

    def _submit_upload(self, label: str, host_paths: str, src_endpoint: str, dest_endpoint: str,
                       manifest_files: list):
        """
        Submit a transfer to Globus

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param dest_endpoint: destination endpoint UUID
        :param manifest_files: manifest of all the files to be transferred
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

        for line in manifest_files:
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
