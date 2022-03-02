import logging
import os
import shutil
import uuid
from ..datatransfer_protocol import DataTransferAPIError, SiteTransfer

logger = logging.getLogger(__package__)


class DataTransfer:
    """
    Represents the default local transfer method for running in a local instance. This assumes that rsync is installed
    on the host.
    """
    def _upload(self, src, dest, inode_type):
        logger.info(f"running transfer: {src} -> {dest}")
        if not inode_type == "D":
            try:
                shutil.copy2(src, dest)
            except IOError as io_error:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy2(src, dest)
        else:
            if not os.path.exists(dest):
                    shutil.copytree(src, dest)

    def submit_upload(self, metadata, manifest_files):
        try:
            for line in manifest_files:
                line = line.decode("UTF-8")
                source_path, dest_path, inode_type = line.strip("\n").split("\t")
                logger.info(f"COPYING THE FILE {source_path} to {dest_path}")
                self._upload(source_path, dest_path, inode_type)
            upload_task_id = uuid.uuid4()
            return str(upload_task_id)
        except Exception as e:
            raise DataTransferAPIError("Problem reading file")

    def submit_download(self, metadata, source_dir, dest_dir):
        if os.path.isdir(source_dir):
            try:
                shutil.copytree(source_dir, dest_dir)
            except IOError as io_error:
                logging.error(io_error)

    def cancel_transfer(self, task_id):
        return "Cannot cancel local python transfer {task_id}"

    def transfer_status(self, task_id):
        return SiteTransfer.status.succeeded
