import os
import globus_sdk
import logging
from jaws_site import config


logger = logging.getLogger(__package__)


class Globus:
    def __init__(self):
        conf = config.conf
        self.globus_root_dir = conf.get("GLOBUS", "root_dir")
        self.globus_default_dir = conf.get("GLOBUS", "default_dir")
        self.uploads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "uploads_subdirectory")
        )
        self.downloads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "downloads_subdirectory")
        )
        self.globus_endpoint = conf.get("GLOBUS", "endpoint_id")

    def _authorize_transfer_client(self, token):
        client_id = config.conf.get("GLOBUS", "client_id")
        client = globus_sdk.NativeAppAuthClient(client_id)
        authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
        return globus_sdk.TransferClient(authorizer=authorizer)

    def transfer_status(self, transfer_rt, task_id):
        """
        Query Globus transfer service for transfer task status.
        """
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
            task = transfer_client.get_task(task_id)
            globus_status = task["status"]
        except Exception:
            logger.exception("Failed to check Globus upload status", exc_info=True)
            raise
        return globus_status

    def transfer_folder(self, transfer_rt, src_dir, dest_endpoint, dest_dir, label):
        """
        Send folder (recursively) via Globus
        """
        logger.debug(f"Xfer {src_dir}")
        if not src_dir.startswith(self.globus_root_dir):
            logger.error(f"Dir is not accessible via Globus: {src_dir}")
            return None
        rel_src_dir = os.path.relpath(src_dir, self.globus_default_dir)
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to get Globus transfer client for {src_dir}", exc_info=True
            )
            return None
        try:
            tdata = globus_sdk.TransferData(
                transfer_client,
                self.globus_endpoint,
                dest_endpoint,
                label=label,
                sync_level="exists",
                verify_checksum=False,
                preserve_timestamp=True,
                notify_on_succeeded=False,
                notify_on_failed=False,
                notify_on_inactive=False,
                skip_activation_check=True,
            )
            if self.globus_root_dir == "/":
                tdata.add_item(src_dir, dest_dir, recursive=True)
            else:
                tdata.add_item(rel_src_dir, dest_dir, recursive=True)
        except Exception:
            logger.warning(
                f"Failed to prepare download manifest for {src_dir}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Xfer failed for {src_dir}", exc_info=True,
            )
            return None
        return transfer_result["task_id"]
