import logging
import globus_sdk
from globus_sdk import GlobusAPIError
from datetime.datetime import utcnow
from sqlalchemy.exc import SQLAlchemyError
import jaws_central.config
from jaws_central.globus import GlobusService


logger = logging.getLogger(__package__)


class XferQueueError(Exception):
    pass


class XferNotFoundError(XferQueueError):
    pass


class XferQueue:
    """
    Central file transfer (xfer) service.
    """

    def __init__(self, session, config=None)
        """
        Constructor
        :param session: database handle
        :type session: sqlalchemy.session
        :param config: configuration object
        :type config: jaws_central.config
        """
        self.session = session
        self.config = config if config else jaws_central.config.conf
        self.max_globus_queue_size = self.config.get("GLOBUS", "max_queue_size", 10)
        self.globus_transfer_client = self._globus_transfer_client()

    def _globus_transfer_client(self):
        """
        Helper method to create a Globus transfer client using
        client id and client secret for credentials.

        :return globus_sdk.TransferClient:
        """
        scopes = "urn:globus:auth:scope:transfer.api.globus.org:all"
        try:
            client = globus_sdk.ConfidentialAppAuthClient(self.config.get("GLOBUS", "client_id"),
                                                          self.config.get("GLOBUS", "client_secret"))
            authorizer = globus_sdk.ClientCredentialsAuthorizer(client, scopes)
        except GlobusAPIError as error:
            logger.error(f"Unable to get Globus transfer client: {error}")
            raise
        else:
            return globus_sdk.TransferClient(authorizer=authorizer)

    def _globus_submit_transfer(self, label, src_endpoint, src_host_path, dest_endpoint, dest_host_path, manifest) -> str:
        tdata = globus_sdk.TransferData(
            self.globus_transfer_client,
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
        for src_path, dest_path, inode_type in manifest:
            logger.debug(f"add transfer: {src_path} -> {dest_path}")
            virtual_src_path = self.virtual_transfer_path(src_path, src_host_path)
            virtual_dest_path = self.virtual_transfer_path(dest_path, dest_host_path)
            if inode_type == "D":
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=True)
            else:
                tdata.add_item(virtual_src_path, virtual_dest_path, recursive=False)
        try:
            transfer_result = self.transfer_client.submit_transfer(tdata)
        except GlobusAPIError as error:
            logger.error(f"Error submitting transfer: {error}")
            raise
        else:
            return transfer_result["task_id"]

    def _globus_transfer_status(self, globus_task_id):
        """
        Query Globus transfer service for transfer task status.
        Uploads are done only by JAWS's globus credentials (i.e. confidential app).

        :param task_id: Globus task id
        :return: str
        """
        try:
            task = self.globus_transfer_client.get_task(globus_task_id)
        except GlobusAPIError as error:
            logger.error(f"Error getting transfer status: {error}")
        else:
            return task["status"]

    def _globus_cancel_transfer(self, globus_task_id: str) -> None:
        """Cancel a transfer.

        :param globus_task_id: Globus' transfer task id
        :type globus_task_id: str
        """
        try:
            self.globus_transfer_client.cancel_task(globus_task_id)
        except GlobusAPIError as error:
            logger.error(f"Error cancelling Globus transfer, {globus_task_id}: {error}")
            raise

    def _db_insert_xfer(self, xfer):
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

    def _db_select_xfer(self, xfer_id: int):
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

    def _db_commit(self, xfer):
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
            raise ValueError(f"Path, {full_path}, is not accessible via Globus endpoint")
        return os.path.join("/", os.path.relpath(full_path, host_path))

    def submit_transfer(
        self, user_id, label, src_endpoint, src_host_path, dest_endpoint, dest_host_path, manifest
    ):
        """
        Queue a transfer to Globus but do not submit.

        :param label: label for the data transfer
        :param host_path: source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
        :param src_endpoint: globus source endpoint UUID
        :param src_host_path: host path (basepath) of the endpoint, used to construct virtual paths
        :param dest_endpoint: destination endpoint UUID
        :param dest_host_path: host path (basepath) of the endpoint, used to construct virtual paths
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
            manifest=manifest_virtual_paths,
            created=utcnow(),
        )
        return self._db_insert_xfer(xfer)

    def transfer_status(self, xfer_id):
        """
        Get the status a transfer.  It returns the status from the object's db, rather than querying Globus, as
        it is assumed this object is periodically updated via the update method.

        :param xfer_id: primary key for transfer task
        :type xfer_id: int
        :return: status
        :rtype: str
        """
        xfer = self._db_select_xfer(xfer_id)
        return xfer.status

    def cancel_transfer(self, xfer_id):
        """
        Cancel a transfer.

        :param xfer_id: XferQueue's PK of the xfer task
        :type xfer_id: int
        :return: None
        """
        xfer = self._db_select_xfer(xfer_id)
        if xfer.status in active_states:
            self._globus_cancel_transfer(xfer.globus_transfer_id)
        xfer.status = "cancelled"
        xfer.updated = utcnow()
        self._db_commit()

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
