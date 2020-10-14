"""
JAWS Analysis Service API

This service stores persistent Run information in a db and interacts with Cromwell.
"""

import logging
import os
import requests
import globus_sdk
from jaws_run import config
from jaws_run.cromwell import Cromwell
from jaws_rpc.rpc_client import RpcClient, RpcError


# config and logging must be initialized before importing this module
cromwell = Cromwell(config.conf.get("CROMWELL", "url"))
logger = logging.getLogger(__package__)


class GlobusError(Exception):
    pass


class Transfer:
    """
    Class representing a Globus transfer task.
    It is used by Run to start upload/download transfers and check if they succeed.
    """

    user_service_client = None  # rpc client for jaws-user service
    transfer_client = None  # globus transfer client obj

    def __init__(self, **kwargs):
        if "transfer_id" in kwargs.items:
            # an existing Globus transfer task
            self.transfer_id = kwargs.get("transfer_id")
            self.__load__()
        else:
            # a new Globus transfer task
            self.submit(**kwargs)

    def __load__(self):
        pass  # TODO

    def submit(self, **kwargs):
        """
        Submit a transfer to Globus, get transfer task id.
        """
        transfer_rt = self.get_transfer_refresh_token
        transfer_task_id = self.__transfer_folder(
            f"Run {self.run_id}",
            transfer_rt,
            self.data.cromwell_workflow_dir,
            self.data.output_endpoint,
            self.data.output_dir,
        )
        if transfer_task_id:
            self.data.download_task_id = transfer_task_id
            self.__update_status__(
                "downloading", f"download_task_id={self.data.download_task_id}"
            )
            self.session.commit()
        self.transfer_id = None  # TODO

    def user_service(self):
        if not self.user_service_client:
            self.user_service_client = RpcClient.client(
                config.conf.user_service_params()
            )
        return self.user_service_client

    def transfer_client(self):
        """
        Get a Globus transfer client, authorized to transfer files on the user's behalf.
        """
        if not self.transfer_client:
            # get Globus token from JAWS User service
            user_id = self.get_user_id()
            user_service = self.get_user_service()
            response = user_service.call("get_transfer_token", {"user_id": user_id})
            token = response["result"]

            # get Globus transfer client
            client = globus_sdk.NativeAppAuthClient(self.client_id)
            authorizer = globus_sdk.RefreshTokenAuthorizer(token, client)
            self.transfer_client = globus_sdk.TransferClient(authorizer=authorizer)
        return self.transfer_client

    def status(self):
        """
        Query Globus transfer service for transfer task status.
        """
        transfer_client = self.get_transfer_client()
        try:
            task = transfer_client.get_task(self.transfer_id)
            self.status = task["status"]
        except Exception:
            logger.exception("Failed to check Globus upload status", exc_info=True)
            raise
        return self.status

    def __transfer_folder(self, label, transfer_rt, src_dir, dest_endpoint, dest_dir):
        """
        Recursively transfer folder via Globus
        :param label: Label to attach to transfer (e.g. "Run 99")
        :type label: str
        :param transfer_rt: User's Globus transfer refresh token
        :type transfer_rt: str
        :param src_dir: Folder to transfer
        :type src_dir: str
        :param dest_endpoint: Globus endpoint for destination
        :type dest_endpoint: str
        :param dest_dir: Destination path
        :type dest_dir: str
        :return: Globus transfer task id
        :rtype: str
        """
        logger.debug(f"Globus xfer {label}")
        if not src_dir.startswith(self.globus_root_dir):
            logger.error(f"Dir is not accessible via Globus: {src_dir}")
            return None
        rel_src_dir = os.path.relpath(src_dir, self.globus_default_dir)
        try:
            transfer_client = self._authorize_transfer_client(transfer_rt)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to get Globus transfer client to xfer {label}", exc_info=True
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
                f"Failed to prepare download manifest for {label}", exc_info=True
            )
        try:
            transfer_result = transfer_client.submit_transfer(tdata)
        except globus_sdk.GlobusAPIError:
            logger.warning(
                f"Failed to download results with Globus for {label}", exc_info=True,
            )
            return None
        return transfer_result["task_id"]
