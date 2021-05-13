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

    def __init__(self, session):
        """
        Constructor
        :param session: database handle
        :type session: sqlalchemy.session
        """
        self.session = session
        self.globus = GlobusService()
        self.max_globus_queue_size = config.conf.get("GLOBUS", "max_queue_size", 10)

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
        xfer = Xfer(
            run_id=run_id,
            status_from=status_from,
            status_to=status_to,
            timestamp=timestamp,
            reason=reason,
        )
        return self._insert_xfer(xfer)

    def _insert_xfer(self, xfer):
        """
        Insert a xfer into the RDb and return the primary key.
        :param xfer: xfer ORM model
        :type xfer: sqlalchemy.model
        """
        try:
            self.session.add(xfer)
            self.session.commit()
        except Exception as error:
            self.session.rollback()
            return None
        return xfer.id

    def _select_xfer(self, xfer_id):
        """
        Select xfer record from db.
        :param xfer_id: primary key for record
        :type xfer_id: int
        :return: ORM model for xfer
        :rtype: sqlalchemy.model
        """
        try:
            xfer = self.session.query(Xfer).get(xfer_id)
        except SQLAlchemyError as error:
            logger.error(error)
            return None
        return xfer

    def transfer_status(self, xfer_id):
        """
        Get the status a transfer.  It returns the status from the object's db, rather than querying Globus, as
        it is assumed this object is periodically updated via the update method.

        :param xfer_id: XferQueue's PK of the xfer task
        :type xfer_id: int
        """
        xfer = self._select_xfer(xfer_id)
        if xfer:
            return xfer.status
        else:
            return None

    def virtual_transfer_path(self, full_path, host_path):
        return self.globus.virtual_transfer_path(full_path, host_path)

    def cancel_transfer(self, xfer_id):
        """
        Cancel a transfer.

        :param xfer_id: XferQueue's PK of the xfer task
        :type xfer_id: int
        :return: None
        """
        xfer = self._select_xfer(xfer_id)
        if xfer.status in active_states:
            self._cancel_globus_transfer(xfer)

    def update(self):
        """
        Update status of active transfers, determine queue size, and submit the appropriate number of new transfers to Globus.
        This is called by the jaws_central.daemon periodically (e.g. every 10 sec).
        """
        self.update_status()
        delta = self.max_globus_queue_size - self.num_active_transfers()
        if delta > 0:
            self._submit_xfers_to_globus(delta)

    def update_status(self):
        """
        Query Globus and update the status of transfer tasks.  This should be called periodically (e.g. every 10s).
        """
        active_transfers = self._active_transfers()

        globus_active_transfers = self.globus.task_list(
            100, status="ACTIVE", type="TRANSFER"
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
        return self.session.query(Xfer).filter(Xfer.status == "transferring").all()

    def _submit_xfers_to_globus(self, num_xfers: int):
        """
        Submit the indicated number of transfer tasks to Globuds.
        :param num_xfers: Number of transfer tasks to submit
        :type num_xfers: int
        """
        xfers = self._select_highest_priority_xfers(num_xfers)
        for xfer in xfers:
            self._submit_xfer_to_globus(xfer)

    def _submit_xfer_to_globus(self, xfer):
        """
        Submit one xfer to Globus and update the db record with the transfer task id.
        :param xfer: xfer ORM object
        :type xfer: sqlalachemy.model
        """
        xfer.transfer_task_id = self.globus.submit_transfer()
        self.session.commit()

    def _select_highest_priority_xfers(self, num_xfers: int):
        """
        Return the indicated number of transfer tasks with the highest priority.
        Intra-site uploads are always processed first (because they are no-op).
        Otherwise, the queue is FIFO, where the earliest tasks have highest priority.
        :param num_xfers: Number of transfer tasks to return
        :type num_xfers: int
        """
        xfers = (
            self.session.query(Xfer)
            .filter(Xfer.src_endpoint_id == Xfer.dest_endpoint_id)
            .filter(Xfer.status == "created")
            .filter(Xfer.type == "upload")
            .order_by(Xfer.submitted)
            .limit(num_xfers)
            .all()
        )
        n = num_xfers - len(xfers)
        if n > 0:
            query = (
                self.session.query(Xfer)
                .filter(Xfer.src_endpoint_id != Xfer.dest_endpoint_id)
                .filter(Xfer.status == "created")
                .order_by(Xfer.submitted)
                .limit(n)
                .all()
            )
            xfers.extend(query)
        return xfers
