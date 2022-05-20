"""
Transfer is a collection of files/folders to transfer (e.g. via Globus, FTP, etc.)
Items are stored in a relational database.
"""

import logging
import os
import pathlib
from datetime import datetime
from sqlalchemy.exc import SQLAlchemyError
import json
import boto3
from jaws_site import config, models


logger = logging.getLogger(__package__)


def mkdir(folder):
    pathlib.Path(folder).mkdir(parents=True, exist_ok=True)


class TransferError(Exception):
    # base class for all errors in this package
    pass


class TransferDbError(TransferError):
    pass


class TransferNotFoundError(TransferError):
    pass


class TransferValueError(TransferError):
    pass


class Transfer:
    """Class representing a transfer (set of files to transfer) associated with one Run."""

    def __init__(self, session, data):
        """
        Initialize transfer object.
        :param session: database handle
        :ptype session: SqlAlchemy.session
        :param data: ORM record for this object
        :ptype: sqlalchemy.ext.declarative.declarative_base
        """
        self.session = session
        self.data = data

    @classmethod
    def from_params(cls, session, params):
        """Create new transfer from parameter values and save in RDb."""
        manifest_json = "[]"
        if "manifest" in params:
            assert type(params["manifest"]) == list
            manifest_json = json.dumps(params["manifest"])
        elif "manifest_json" in params:
            assert (params["manifest_json"]) == str
            manifest_json = params["manifest_json"]
        try:
            data = models.Transfer(
                id=params["transfer_id"],
                status="queued",
                src_base_dir=params["src_base_dir"],
                dest_base_dir=params["dest_base_dir"],
                manifest_json=manifest_json,
            )
        except SQLAlchemyError as error:
            raise TransferValueError(
                f"Error creating model for new Transfer: {params}: {error}"
            )
        try:
            session.add(data)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            raise TransferDbError(error)
        else:
            return cls(session, data)

    @classmethod
    def from_id(cls, session, transfer_id: int):
        """Select existing Transfers record from RDb by primary key"""
        try:
            data = (
                session.query(models.Transfer)
                .filter_by(id=int(transfer_id))
                .one_or_none()
            )
        except SQLAlchemyError as error:
            raise TransferDbError(f"Error selecting Transfer {transfer_id}: {error}")
        except Exception as error:
            raise TransferError(f"Error selecting Transfer {transfer_id}: {error}")
        else:
            if data is None:
                raise TransferNotFoundError(f"Transfer {transfer_id} not found")
            else:
                return cls(session, data)

    def status(self) -> str:
        """Return the current state of the transfer."""
        return self.data.status

    def manifest(self) -> list:
        return json.loads(self.data.manifest_json)

    def cancel(self) -> None:
        """Cancel a transfer by changing the status in the db to prevent it from being picked up
        by the transfer daemon."""
        # changing the status to cancel will prevent the transfer_daemon from picking up the task but
        # will not cancel a transfer task that has already begun transferring
        if self.data.status == "queued":
            self.update_status("cancelled")
            return True
        else:
            return False

    def update_status(self, new_status) -> None:
        """
        Update Transfers' status.
        """
        logger.info(f"Transfers {self.data.id}: now {new_status}")
        timestamp = datetime.utcnow()
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Transfer {self.data.id}: {error}")
            raise (error)

    def transfer_files(self) -> None:
        """
        Do the transfer and return when done -- the operation is blocking.
        """
        logger.debug(f"Begin transfer {self.data.id}")
        self.update_status("transferring")
        try:
            if self.data.src_base_dir.startswith("s3://"):
                self.s3_download_folder()
            elif self.data.dest_base_dir.startswith("s3://"):
                self.s3_upload()
            else:
                logger.error(
                    f"Transfer {self.data.id} failed because neither src/dest start with s3://"
                )
                self.update_status("failed")
        except IOError as error:
            logger.error(f"Transfer {self.data.id} failed: {error}")
            self.update_status("failed")
        else:
            self.update_status("succeeded")

    def aws_s3_resource(self):
        aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
        aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
        aws_region_name = config.conf.get("AWS", "aws_region_name")
        aws_session = boto3.Session(
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name=aws_region_name,
        )
        s3_resource = aws_session.resource("s3")
        return s3_resource

    def aws_s3_client(self):
        aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
        aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
        client = boto3.client(
            "s3",
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
        )
        return client

    @staticmethod
    def s3_parse_path(full_path):
        full_path = full_path.replace("s3://", "", 1)  # discard S3 identifier
        folders = full_path.split("/")
        s3_bucket = folders.pop(0)
        path = "/".join(folders)
        print(f"S3 BUCKET={s3_bucket}; PATH={path}")
        return s3_bucket, path

    def s3_upload_old(self):
        manifest = self.manifest()
        num_files = len(manifest)
        logger.debug(f"Transfer {self.data.id} begin s3 upload of {num_files} files")
        aws_s3_resource = self.aws_s3_resource()
        s3_bucket, dest_base_dir = self.s3_parse_path(self.data.dest_base_dir)
        bucket_obj = aws_s3_resource.Bucket(s3_bucket)
        for rel_path in manifest:
            src_path = os.path.normpath(os.path.join(self.data.src_base_dir, rel_path))
            dest_path = os.path.normpath(os.path.join(dest_base_dir, rel_path))
            logger.debug(f"S3 upload to {s3_bucket}: {src_path} -> {dest_path}")
            try:
                with open(src_path, "rb") as fh:
                    bucket_obj.upload_fileobj(fh, dest_path)
            except Exception as error:
                raise IOError(error)

    def s3_file_size(self, bucket, file_key, aws_s3_client=None):
        """
        If a file_key exists then return it's size, otherwise None.
        :param bucket: name of the S3 bucket
        :ptype bucket: str
        :param file_key: identifier for the file (similar to path)
        :ptype file_key: str
        :param aws_s3_client: S3 client object
        :ptype aws_s3_client: boto3.client
        :return: size
        :rtype: int
        """
        if aws_s3_client is None:
            aws_s3_client = self.aws_s3_client()
        try:
            result = aws_s3_client.list_objects_v2(Bucket=bucket, Prefix=file_key)
        except boto3.ClientError as error:
            logger.error(f"Error getting S3 obj stats for {file_key}: {error}")
            raise
        size = None
        if "Contents" in result:
            contents = result["Contents"]
            if len(contents) > 1:
                raise ValueError(
                    f"Unexpectedly >1 obj for S3 key {file_key}: {contents}"
                )
            file_obj = contents[0]
            size = file_obj["Size"]
            # mtime = file_obj["LastModified"]
        return size

    def s3_upload(self):
        manifest = self.manifest()
        num_files = len(manifest)
        logger.debug(f"Transfer {self.data.id} begin s3 upload of {num_files} files")
        aws_client = self.aws_s3_client()
        aws_s3_resource = self.aws_s3_resource()
        s3_bucket, dest_base_dir = self.s3_parse_path(self.data.dest_base_dir)
        bucket_obj = aws_s3_resource.Bucket(s3_bucket)
        for rel_path in manifest:
            src_path = os.path.normpath(os.path.join(self.data.src_base_dir, rel_path))
            src_file_size = os.path.getsize(src_path)
            dest_path = os.path.normpath(os.path.join(dest_base_dir, rel_path))
            dest_file_size = self.s3_file_size(s3_bucket, dest_path, aws_client)
            if src_file_size == dest_file_size:
                logger.debug(f"S3 upload: Skipping cached file {dest_path}")
            else:
                logger.debug(f"S3 upload to {s3_bucket}: {src_path} -> {dest_path}")
                try:
                    with open(src_path, "rb") as fh:
                        bucket_obj.upload_fileobj(fh, dest_path)
                except Exception as error:
                    raise IOError(error)

    def s3_download(self):
        manifest = self.manifest()
        num_files = len(manifest)
        logger.debug(f"Transfer {self.data.id} begin s3 download of {num_files} files")
        aws_s3_resource = self.aws_s3_resource()
        s3_bucket, src_base_dir = self.s3_parse_path(self.data.src_base_dir)
        bucket_obj = aws_s3_resource.Bucket(s3_bucket)
        for rel_path in manifest:
            src_path = os.path.normpath(os.path.join(src_base_dir, rel_path))
            dest_path = os.path.normpath(
                os.path.join(self.data.dest_base_dir, rel_path)
            )
            logger.debug(
                f"S3 download from {s3_bucket}: {rel_path} as {src_path} -> {dest_path}"
            )
            dest_folder = os.path.dirname(dest_path)
            try:
                mkdir(dest_folder)
            except IOError as error:
                logger.error(f"Unable to make download dir, {dest_folder}: {error}")
            try:
                with open(dest_path, "wb") as fh:
                    bucket_obj.download_fileobj(src_path, fh)
            except Exception as error:
                raise IOError(error)

    def s3_download_folder(self):
        aws_s3_client = self.aws_s3_client()
        s3_bucket, src_base_dir = self.s3_parse_path(self.data.src_base_dir)
        paginator = aws_s3_client.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=s3_bucket, Prefix=src_base_dir):
            for obj in page["Contents"]:
                rel_path = obj["Key"]
                size = obj["Size"]
                dest_rel_path = rel_path.removeprefix(src_base_dir)
                dest_path = os.path.normpath(
                    os.path.join(self.data.dest_base_dir, f"./{dest_rel_path}")
                )
                logger.debug(f"S3 download {rel_path} -> {dest_path}")
                if rel_path.endswith("/") and size == 0:
                    try:
                        mkdir(dest_path)
                    except Exception as error:
                        msg = f"Unable to make download dir, {dest_path}: {error}"
                        logger.error(msg)
                        raise IOError(msg)
                else:
                    try:
                        basedir = os.path.dirname(dest_path)
                        mkdir(basedir)
                    except Exception as error:
                        msg = f"Unable to make download dir, {basedir}: {error}"
                        logger.error(msg)
                        raise IOError(msg)
                    try:
                        aws_s3_client.download_file(s3_bucket, rel_path, dest_path)
                    except Exception as error:
                        msg = f"S3 download error, {rel_path}: {error}"
                        logger.error(msg)
                        raise IOError(msg)


def check_queue(session) -> None:
    """
    Check the transfer queue and start the oldest transfer task, if any.  This only does one task because
    transfers typically take many minutes and the queue may change (e.g. a transfer is cancelled).
    """
    try:
        rows = (
            session.query(models.Transfer)
            .filter(models.Transfer.status == "queued")
            .order_by(models.Transfer.id)
            .limit(1)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
    if len(rows):
        transfer = Transfer(session, rows[0])
        transfer.transfer_files()
