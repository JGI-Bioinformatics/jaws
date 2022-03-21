import logging
import os
import shutil
import uuid
from ..datatransfer_protocol import DataTransferAPIError

logger = logging.getLogger(__package__)


class DataTransfer:
    """
    Represents the default local transfer method for running in a local instance. It uses the built-in
    shutil python module for creating data transfers. This should not be used for production-sized data
    transfers.
    """

    def _copy(self, src, dest, inode_type):
        logger.info(f"Local copy: {src} -> {dest}")
        if inode_type == "D":
            # folders are copied recursively
            if not os.path.exists(dest):
                shutil.copytree(src, dest)
        else:
            try:
                shutil.copy2(src, dest)
            except IOError:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy2(src, dest)

    def submit_transfer(self, manifest, **kwargs):
        """
        Copies files and/or folders listed in transfer manifest table.
        """
        try:
            for line in manifest:
                line = line.decode("UTF-8")
                source_path, dest_path, inode_type = line.strip("\n").split("\t")
                logger.info(f"COPYING THE FILE {source_path} to {dest_path}")
                self._copy(source_path, dest_path, inode_type)
        except Exception as e:
            logging.error(e, exc_info=True)
            raise DataTransferAPIError("Problem reading file")
        else:
            # an ID is expected by the API, so return a random string
            transfer_task_id = uuid.uuid4()
            return str(transfer_task_id)

    def cancel_transfer(self, task_id):
        """
        Notify the user that cancel cannot happen.

        Since we are using shutil in the main thread, we cannot really cancel
        a file transfer. This returns a notification that the transfer cannot be cancelled.
        """
        _ = task_id  # task_id is not created for local transfer
        return "Cannot cancel local python transfer {task_id}"

    def transfer_status(self, task_id):
        """
        Return the status code for a successful transfer.
        Assumes that the transfer is successful.
        """
        _ = task_id  # task_id is not created for local transfer
        # status will always be successful since local copy is a blocking, synchronous operation
        # (i.e. it either completed or failed &submit_transfer method already)
        return "succeeded"
