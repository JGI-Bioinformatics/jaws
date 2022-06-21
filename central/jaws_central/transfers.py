import logging
import json
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from jaws_rpc import rpc_index
from jaws_central import config
from jaws_central.globus import GlobusService
from jaws_central import models

logger = logging.getLogger(__package__)


class TransferError(Exception):
    # base class for all errors
    pass


class TransferInputError(TransferError):
    # base class for all user input data errors (e.g. ValueError);
    # these require user to fix their Run
    pass


class TransferSystemError(TransferError):
    # base class for all system related errors (i.e. not user data);
    # these require investigation by jaws admins
    pass


class TransferNotFoundError(TransferInputError):
    pass


class TransferDbError(TransferSystemError):
    pass


class TransferRpcError(TransferSystemError):
    pass


class TransferGlobusError(TransferSystemError):
    pass


class TransferSiteError(TransferSystemError):
    pass


class Transfer:
    def __init__(self, session, data):
        self.session = session
        self.data = data
        self.src_site_config = config.conf.get_site(self.data.src_site_id)
        self.dest_site_config = config.conf.get_site(self.data.dest_site_id)

    @classmethod
    def from_params(cls, session, params):
        """Create new transfer and save in RDb"""
        try:
            data = models.Transfer(
                status="created",
                src_site_id=params["src_site_id"],
                src_base_dir=params["src_base_dir"],
                dest_site_id=params["dest_site_id"],
                dest_base_dir=params["dest_base_dir"],
                manifest_json=json.dumps(params["manifest"]),
            )
        except SQLAlchemyError as error:
            raise (f"Error creating model for new Transfer: {params}: {error}")
        try:
            session.add(data)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            raise (error)
        else:
            return cls(session, data)

    @classmethod
    def from_id(cls, session, id):
        """Select transfer record from Rdb"""
        try:
            data = session.query(models.Transfer).get(id)
        except IntegrityError as error:
            logger.error(f"Transfer {id} not found", error)
            raise TransferNotFoundError(f"Transfer {id} not found")
        except SQLAlchemyError as error:
            logger.error(f"Unable to select Transfer {id}", error)
            raise TransferDbError("Error selecting Transfer {id}: {error}")
        else:
            return cls(session, data)

    def status_rpc(self, site_id):
        """
        Ask jaws-site for current status.
        :param site_id: The jaws-site that is doing the transfer.
        :ptype site_id: str
        :return: new status
        :rtype: str
        """
        rpc_client = rpc_index.rpc_index.get_client(site_id)
        try:
            response = rpc_client.request(
                "transfer_status", {"transfer_id": self.data.id}
            )
        except Exception as error:
            msg = f"RPC error(1) for check status, transfer {self.data.id}: {error}"
            logger.error(msg)
            raise TransferRpcError(msg)
        else:
            if "error" in response:
                msg = f"RPC error(2) for check status, transfer {self.data.id}: {response['error']['message']}"
                raise TransferRpcError(msg)
            else:
                return response["result"]["status"], response["result"]["reason"]

    def status_globus(self):
        """
        Ask Globus service for current status
        :return: new status
        :rtype: str
        """
        globus_client = GlobusService()
        try:
            new_status, reason = globus_client.transfer_status(
                self.data.globus_transfer_id
            ).lower()
        except Exception as error:
            msg = f"Globus error checking status of transfer {self.data.id}: {error}"
            raise TransferGlobusError(msg)
        else:
            return new_status, reason

    def status(self) -> str:
        if self.data.status in ["submission failed", "failed", "succeeded", "cancelled"]:
            # terminal states don't change, no need to query
            return self.data.status, self.data.reason
        site_id = self.responsible_site_id()
        original_status = self.data.status
        new_status = None
        try:
            if site_id:
                new_status, reason = self.status_rpc(site_id)
            else:
                new_status, reason = self.status_globus()
        except TransferError as error:
            logger.error(f"Unable to retrieve transfer status: {error}")
        else:
            if new_status != original_status:
                logger.debug(f"Task {self.data.id} status = {new_status}")
                self.update_status(new_status, reason)
        return self.data.status, self.data.reason

    def manifest(self) -> list:
        """
        The manifest JSON string stored in the RDb must be deserialized.
        """
        files = json.loads(self.data.manifest_json)
        logger.debug(f"MANIFEST = {files}")
        return files

    def cancel(self) -> None:
        if self.data.status == "created":
            self.update_status("cancelled")

    #        elif self.data.status == "queued":
    #            if self.data.globus_transfer_id is not None:
    #                self.cancel_globus_transfer()
    #            else:
    #                self.cancel_rpc_transfer()
    #            self.update_status("cancelled")

    def update_status(self, new_status) -> None:
        try:
            self.data.status = new_status
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.error(f"Unable to update Transfer {self.data.id}: {error}")
            raise TransferDbError(error)
        else:
            logger.info(f"Transfer {self.data.id}: now {new_status}")

    def responsible_site_id(self):
        """
        Which jaws-site is responsible for this transfer, if any.
        When None, that means this (jaws-central) service is responsible for this transfer.
        """
        site_id = None  # if None, then this service (jaws-central) is responsible for the transfer
        if self.data.src_site_id == self.data.dest_site_id:
            raise TransferError(f"No transfer expected for intra-site transfers; Transfer {self.data.id}")
        elif (
            self.src_site_config["globus_endpoint"]
            and self.dest_site_config["globus_endpoint"]
        ):
            pass  # do not send; central will submit to globus
        elif self.data.src_base_dir.startswith(
            "s3://"
        ) and self.data.dest_base_dir.startswith("s3://"):
            # we don't currently support retrieving user input data from S3 buckets
            raise TransferError("S3:S3 transfers are not supported")
        elif self.data.src_base_dir.startswith("s3://"):
            # download from S3->NFS must be done by destination site (with NFS access)
            site_id = self.data.dest_site_id
        elif self.data.dest_base_dir.startswith("s3://"):
            # upload from NFS->S3 must be done by source site (with NFS access)
            site_id = self.data.src_site_id
        else:
            raise TransferError(f"No transfer method known for Transfer {self.data.id}")
        return site_id

    def submit_transfer(self) -> None:
        """
        If the file transfer is between jaws-sites with Globus endpoints, then jaws-central can
        submit to Globus service directly.  For AWS-S3 copy operations, the jaws-site
        must do operation because central doesn't have access to the jaws-sites' file systems.
        REST requests to Globus are sent via the SDK.  Requests to the jaws-site are sent via RPC.
        """
        manifest = self.manifest()
        src_base_dir = self.data.src_base_dir
        dest_base_dir = self.data.dest_base_dir
        responsible_site_id = self.responsible_site_id()

        if responsible_site_id:
            # jaws-sites handle AWS-S3 copy operations
            params = {
                "transfer_id": self.data.id,
                "src_site_id": self.data.src_site_id,
                "src_base_dir": self.data.src_base_dir,
                "dest_site_id": self.data.dest_site_id,
                "dest_base_dir": dest_base_dir,
                "manifest": manifest,
            }
            try:
                status = self._rpc(responsible_site_id, "submit_transfer", params)
            except TransferRpcError as error:
                logger.error(f"Submit site transfer {self.data.id} RPC error: {error}")
                # do not update status; keep retrying
            except TransferSiteError as error:
                logger.error(f"Submit site transfer {self.data.id} rejected: {error}")
                self.update_status("submission failed")
            except Exception as error:
                logger.error(f"Site transfer {self.data.id} unexpected error: {error}")
                self.update_status("submission failed")
            else:
                logger.debug(f"Submitted RPC transfer {self.data.id}: {status}")
                self.update_status("queued")
        else:
            # jaws-central with submit to Globus
            label = f"Transfer {self.data.id}"
            src_endpoint = self.src_site_config["globus_endpoint"]
            src_host_path = self.src_site_config["globus_host_path"]
            dest_endpoint = self.dest_site_config["globus_endpoint"]
            dest_host_path = self.dest_site_config["globus_host_path"]
            try:
                globus_client = GlobusService()
                globus_transfer_id = globus_client.submit_transfer(
                    label,
                    src_endpoint,
                    src_host_path,
                    src_base_dir,
                    dest_endpoint,
                    dest_host_path,
                    dest_base_dir,
                    manifest,
                )
            except Exception as error:
                logger.error(f"Globus transfer {self.data.id} failed: {error}")
                self.update_status("submission failed")
                raise TransferRpcError(error)
            else:
                self.data.globus_transfer_id = globus_transfer_id
                self.update_status("queued", f"globus_transfer_id={globus_transfer_id}")

    def _rpc(self, site_id, function, params):
        rpc_client = rpc_index.rpc_index.get_client(site_id)
        try:
            response = rpc_client.request(function, params)
        except Exception as error:
            reason = f"RPC {function} failed: {error}"
            logger.exception(reason)
            raise TransferRpcError(reason)
        if "error" in response:
            reason = (
                f"Site {site_id} RPC {function} error: {response['error']['message']}"
            )
            logger.error(reason)
            raise TransferSiteError(reason)
        return response["result"]
