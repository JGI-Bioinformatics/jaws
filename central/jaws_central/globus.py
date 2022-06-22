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

    def _create_transfer_client(self):
        """
        Creates the transfer client using client id and client secret for
        credentials
        :return: globus_sdk.TransferClient
        """
        client = globus_sdk.ConfidentialAppAuthClient(
            self.globus_config.get("GLOBUS", "client_id"),
            self.globus_config.get("GLOBUS", "client_secret"),
        )
        scopes = "urn:globus:auth:scope:transfer.api.globus.org:all"
        authorizer = globus_sdk.ClientCredentialsAuthorizer(client, scopes)
        return globus_sdk.TransferClient(authorizer=authorizer)

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
        if "fatal_error" in task and "description" in task["fatal_error"]:
            reason = task["fatal_error"]["description"]
        return globus_status, reason

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
            raise ValueError(
                f"Path, {full_path}, is not accessible via Globus endpoint"
            )
        return os.path.join("/", os.path.relpath(full_path, host_path))

    def submit_transfer(
        self,
        label,
        src_endpoint,
        src_host_path,
        src_base_dir,
        dest_endpoint,
        dest_host_path,
        dest_base_dir,
        manifest,
    ):
        """
        Submit a transfer to Globus

        :param label: label for the data transfer
        :param src_endpoint: globus source endpoint UUID
        :param src_host_path: The endpoint basedir is used to modify the paths and create virtual transfer paths.
        :param dest_endpoint: destination endpoint UUID
        :param dest_host_path: The endpoint basedir is used to modify the paths and create virtual transfer paths.
        :param manifest: table of all the files to be transferred
        :ptype manifest: list
        :return:
        """
        logger.debug(f"Globus xfer {label}")
        transfer_client = self._create_transfer_client()
        # sync level 1 : Copy files if the size of the destination does not match the size of the source.
        tdata = globus_sdk.TransferData(
            transfer_client,
            src_endpoint,
            dest_endpoint,
            label=label,
            sync_level=1,
            verify_checksum=True,
            preserve_timestamp=True,
            notify_on_succeeded=False,
            notify_on_failed=True,
            notify_on_inactive=True,
            recursive_symlinks="copy",
        )

        # the manifest is empty for complete download, add path and do recursive
        if len(manifest) == 0:
            virtual_src_path = self.virtual_transfer_path(src_base_dir, src_host_path)
            virtual_dest_path = self.virtual_transfer_path(dest_base_dir, dest_host_path)
            tdata.add_item(virtual_src_path, virtual_dest_path, recursive=True)
        else:
            for relpath in manifest:
                source_path = f"{src_base_dir}/{relpath}"
                dest_path = f"{dest_base_dir}/{relpath}"
                virtual_src_path = self.virtual_transfer_path(source_path, src_host_path)
                virtual_dest_path = self.virtual_transfer_path(dest_path, dest_host_path)
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=False)

        transfer_result = transfer_client.submit_transfer(tdata)
        return transfer_result["task_id"]
