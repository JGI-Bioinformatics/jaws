import logging
import jaws_central.config
from jaws_central.globus import GlobusService


logger = logging.getLogger(__package__)


class XferError(Exception):
    pass


class XferQueue:
    """
    Central file transfer (xfer) service.  Extends GlobusService with a transfer queue.
    """

    def __init__(self):
        self.globus = GlobusService()

    def submit_transfer(
        self, label, host_paths, src_endpoint, dest_endpoint, manifest_file
    ):
        """
        Queue a transfer to Globus but do not submit.

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param dest_endpoint: destination endpoint UUID
        :param manifest_file: manifest of all the files to be transferred
        :return: xfer_id (primary key)
        :rtype: int
        """
        logger.debug(f"Globus xfer {label}")
        # TODO
        return xfer.id

    def _insert_xfer(
        self, src_endpoint, dest_endpoint, manifest, user_id, label, size_gb
    ):
        """
        Insert a xfer into the RDb and return the primary key.

        :param src_endpoint: source Globus endpoint ID
        :type src_endpoint: str
        :param dest_endpoint: destination Globus endpoint ID
        :type dest_endpoint: str
        :param manifest: transfer manifest (list of src, dest paths)
        :type manifest: list
        :param user_id: JAWS user ID (None if not Run data)
        :type user_id: str
        :param label: transfer type (e.g. upload, download, refdata)
        :type label: str
        :param size_gb: total size of all items to transfer, in gigabytes
        :type size_gb: int
        :return: xfer_id (primary key)
        :rtype: int
        """
        # TODO
        return xfer.id

    def transfer_status(self, xfer_id):
        """
        Get the status a transfer.  It returns the status from the object's db, rather than querying Globus, as
        it is assumed this object is periodically updated via the update method.

        :param xfer_id: XferQueue's PK of the xfer task
        :type xfer_id: int
        """
        # TODO
        return status

    def virtual_transfer_path(self, full_path, host_path):
        return self.globus.virtual_transfer_path(full_path, host_path)

    def cancel_transfer(self, xfer_id):
        """
        Cancel a transfer.

        :param xfer_id: XferQueue's PK of the xfer task
        :type xfer_id: int
        :return: None
        """
        # TODO

    def update_status(self):
        """
        Query Globus and update the status of transfer tasks.  This should be called periodically (e.g. every 10s).
        """
        active_transfers = self._active_transfers()

        globus_active_transfers = self.globus.task_list(
            100, status="ACTIVE", type="TRANSFER,DELETE"
        )
        for transfer_task in globus_active_transfers:
            task_id = transfer_task["task_id"]
            globus_status = transfer_task["status"]
            if task_id in active_transfers:
                xfer = active_transfers[task_id]
                if globus_status != xfer.status:
                    xfer.status = globus_status
                    xfer.updated = utcnow()

    def _active_transfers(self):
        """
        Query db for all active transfers.

        :return: Active transfers, where { globus_transfer_task_id : xfer model }
        :rtype: dict
        """
        active_transfers = []  # TODO
        return active_transfers

    def submit_xfers(self, num_xfers: int):
        """
        Submit the indicated number of transfer tasks to Globuds.
        :param num_xfers: Number of transfer tasks to submit
        :type num_xfers: int
        """

    def _select_highest_priority_tasks(self, num_xfers: int):
        """
        Return the indicated number of transfer tasks with the highest priority.
        The queue is FIFO, so the earliest tasks have highest priority.
        :param num_xfers: Number of transfer tasks to return
        :type num_xfers: int
        """
        # TODO
        return xfers
