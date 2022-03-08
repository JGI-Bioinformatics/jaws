import logging
import os
import shutil
import uuid
from ..datatransfer_protocol import DataTransferAPIError, SiteTransfer

logger = logging.getLogger(__package__)


class DataTransfer:
    """
    Represents the default local transfer method for running in a local instance. It uses the built-in
    shutil python module for creating data transfers. This should not be used for production-sized data
    transfers.
    """
    def _upload(self, src, dest, inode_type):
        logger.info(f"running transfer: {src} -> {dest}")
        if inode_type == "D":
            if not os.path.exists(dest):
                shutil.copytree(src, dest)
        else:
            try:
                shutil.copy2(src, dest)
            except IOError:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy2(src, dest)

    def submit_upload(self, metadata, manifest_files):
        """
        Submits transfer uploads from TSV file.
        """
        _ = metadata  # metadata is not needed for the local transfer
        try:
            for line in manifest_files:
                line = line.decode("UTF-8")
                source_path, dest_path, inode_type = line.strip("\n").split("\t")
                logger.info(f"COPYING THE FILE {source_path} to {dest_path}")
                self._upload(source_path, dest_path, inode_type)
            upload_task_id = uuid.uuid4()
            return str(upload_task_id)
        except Exception as e:
            logging.error(e, exc_info=True)
            raise DataTransferAPIError("Problem reading file")

    def submit_download(self, metadata, source_dir, dest_dir):
        """
        Submits a transfer download.

        Transfers a source directory to a destination directory.
        """
        _ = metadata  # metadata is not needed for local transfer
        if os.path.isdir(source_dir):
            try:
                shutil.copytree(source_dir, dest_dir)
            except IOError as io_error:
                logging.error(io_error, exc_info=True)

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
        return SiteTransfer.status.succeeded
