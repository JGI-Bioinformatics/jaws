import logging
import jaws_central.config
from datetime.datetime import utcnow
from sqlalchemy.exc import SQLAlchemyError
from jaws_central.globus import GlobusService


logger = logging.getLogger(__package__)


class XferQueueError(Exception):
    pass


class XferNotFoundError(XferQueueError):
    pass


class XferQueue:
    """
    Central file transfer (xfer) service.  Extends GlobusService with a transfer queue.
    """

    def __init__(self, session, globus_service=None):
        """
        Constructor
        :param session: database handle
        :type session: sqlalchemy.session
        :param globus_service: obj for interacting with Globus REST API
        :type globus_service: jaws_central.globus.GlobusService
        """
        self.session = session
        if globus_service is None:
            globus_service = GlobusService()
        self.globus = globus_service
        self.max_globus_queue_size = config.conf.get("GLOBUS", "max_queue_size", 10)

    def _insert_xfer(self, xfer):
        """
        Insert an xfer into the RDb and return the primary key.
        :param xfer: xfer ORM model
        :type xfer: sqlalchemy.model
        """
        try:
            self.session.add(xfer)
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.error(error)
            raise
        else:
            self.session.commit()
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
            raise
        if xfer:
            return xfer
        else:
            raise XferNotFoundError(f"Xfer {xfer_id} not found")

    def _update_xfer(self, xfer):
        """
        Update an xfer record in the RDb.
        :param xfer: ORM model for a transfer task
        :type xfer: sqlalchemy.model
        """
        try:
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.error(error)
            raise

    def submit_transfer(
        self, user_id, label, host_paths, src_endpoint, dest_endpoint, manifest
    ):
        """
        Queue a transfer to Globus but do not submit.

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param dest_endpoint: destination endpoint UUID
        :param manifest: list of items to transfer where item is tuple of src, dest, type
        :type manifest: list
        :return: xfer_id (primary key)
        :rtype: int
        """
        logger.debug(f"Submit xfer {user_id}:{label}")
        xfer = Xfer(
            user_id=user_id,
            label=label,
            src_endpoint_id=src_endpoint_id,
            dest_endpoint_id=dest_endpoint_id,
            manifest=manifest,
            created=utcnow(),
        )
        return self._insert_xfer(xfer)

    def transfer_status(self, xfer_id):
        """
        Get the status a transfer.  It returns the status from the object's db, rather than querying Globus, as
        it is assumed this object is periodically updated via the update method.

        :param xfer_id: primary key for transfer task
        :type xfer_id: int
        :return: status
        :rtype: str
        """
        xfer = self._select_xfer(xfer_id)
        return xfer.status

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
            self.globus.cancel_transfer(xfer.transfer_task_id)
        xfer.status = "cancelled"
        xfer.updated = utcnow()
        self._update_xfer()

    def update(self):
        """
        Update status of active transfers, determine queue size, and submit the appropriate number of new transfers to Globus.
        This is called by the daemon periodically (e.g. every 10 sec).
        """
        self.update_status()
        delta = self.max_globus_queue_size - self.num_active_transfers()
        if delta > 0:
            self._submit_xfers_to_globus(delta)

    def update_status(self):
        """
        Query Globus and update the status of transfer tasks.  This should be called periodically (e.g. every 10s).
        """
        globus_transfers = self._active_transfers()
        for transfer_task in globus_active_transfers:
            task_id = transfer_task["task_id"]
            globus_status = transfer_task["status"]
            if task_id in active_transfers:
                xfer = active_transfers[task_id]
                if globus_status != xfer.status:
                    xfer.status = globus_status
                    xfer.updated = utcnow()

    def _globus_transfers(self):
        """
        Query db for all unfinished transfers which have been submitted to Globus.

        :return: Globus transfer tasks, where { globus_transfer_task_id : xfer model }
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

    def _select_highest_priority_xfers(self, limit: int):
        """
        Return the indicated number of transfer tasks with the highest priority,
        where the earliest xfers have the highest priority (i.e. FIFO).
        :param limit: Maximum number of transfer tasks to return
        :type limit: int
        """
        return (
            self.session.query(Xfer)
            .filter(Xfer.status == "created")
            .order_by(Xfer.submitted)
            .limit(limit)
            .all()
        )
