"""
Transfer is a collection of files/folders to transfer (e.g. via Globus, FTP, etc.)
Items are stored in a relational database.
"""
import concurrent.futures
import json
import logging
import os
from datetime import datetime

import boto3
from parallel_sync import rsync
from sqlalchemy.exc import SQLAlchemyError

from jaws_site import config, models
import time

logger = logging.getLogger(__package__)


FILES_PER_THREAD = 10000
MAX_ERROR_STRING_LEN = 1024


def mkdir(path, mode=None):
    if mode is None:
        mode = int(config.conf.get("SITE", "folder_permissions"), base=8)
    if os.path.isdir(path):
        return
    else:
        (head, tail) = os.path.split(path)
        mkdir(head, mode)
        if not os.path.exists(path):
            os.mkdir(path)
            os.chmod(path, mode)


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
        if (
            not isinstance(params["transfer_id"], int)
            or not isinstance(params["src_base_dir"], str)
            or not isinstance(params["dest_base_dir"], str)
        ):
            raise TransferValueError
        manifest_json = "[]"
        if "manifest" in params:
            assert isinstance(params["manifest"], list)
            manifest_json = json.dumps(params["manifest"])
        elif "manifest_json" in params:
            assert isinstance(params["manifest_json"], str)
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

    def reason(self) -> str:
        """Return the failure message of the transfer."""
        return self.data.reason

    def manifest(self) -> list:
        """
        Return list of files to transfer -- may be empty list if complete folder is to be transferred.
        """
        manifest = []
        if self.data.manifest_json is not None:
            manifest = json.loads(self.data.manifest_json)
        return manifest

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

    def update_status(self, new_status: str, reason: str = None) -> None:
        """
        Update Transfers' status.
        """
        logger.info(f"Transfers {self.data.id}: now {new_status}")
        if reason is not None:
            reason = reason[:MAX_ERROR_STRING_LEN]
        timestamp = datetime.utcnow()
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            if reason is not None:
                self.data.reason = reason
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
        result = None
        reason = None
        try:
            if self.data.src_base_dir.startswith("s3://"):
                self.s3_download()
            elif self.data.dest_base_dir.startswith("s3://"):
                self.s3_upload()
            else:
                self.local_copy()
        except Exception as error:
            result = "failed"
            reason = str(error)
        else:
            result = "succeeded"

        # session may be stale, so close it to get a new connection
        self.session.remove()
        self.update_status(result, reason)

    def aws_s3_resource(self):
        aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
        aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
        aws_region_name = config.conf.get("AWS", "aws_region_name")
        aws_session = boto3.Session(
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name=aws_region_name,
        )
        return aws_session.resource("s3")

    def aws_s3_client(self):
        aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
        aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
        return boto3.client(
            "s3",
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
        )

    @staticmethod
    def s3_parse_path(full_path):
        full_path = full_path.replace("s3://", "", 1)  # discard S3 identifier
        folders = full_path.split("/")
        s3_bucket = folders.pop(0)
        path = "/".join(folders)
        return s3_bucket, path

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
        aws_s3_client = self.aws_s3_client()
        result = aws_s3_client.list_objects_v2(Bucket=bucket, Prefix=file_key)
        size = None
        if "Contents" in result:
            contents = result["Contents"]
            if len(contents) == 1 and "Size" in contents[0]:
                size = contents[0]["Size"]
                # mtime = file_obj["LastModified"]
            elif len(contents) > 1:
                for rec in contents:
                    if rec["Key"] == file_key and "Size" in rec:
                        size = rec["Size"]
                        # mtime = file_obj["LastModified"]
                        break
        return size

    def s3_upload(self):
        """
        Upload files from NFS->S3.  Skip files which exist and have the same size.
        """
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
            if dest_file_size is not None and src_file_size == dest_file_size:
                logger.debug(f"S3 upload: Skipping existing file {dest_path}")
            else:
                logger.debug(f"S3 upload to {s3_bucket}: {src_path} -> {dest_path}")
                with open(src_path, "rb") as fh:
                    bucket_obj.upload_fileobj(fh, dest_path)

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
            mkdir(dest_folder)
            with open(dest_path, "wb") as fh:
                bucket_obj.download_fileobj(src_path, fh)

    def local_copy(self) -> None:
        """
        Copy files and folders (recursively).
        """
        logger.debug(f"Transfer {self.data.id}: Begin local copy")
        manifest = self.manifest()
        src = f"{self.data.src_base_dir}/"
        if not os.path.isdir(src):
            raise FileNotFoundError(f"Source directory not found: {src}")
        dest = f"{self.data.dest_base_dir}/"
        try:
            mkdir(dest)
        except IOError as error:
            logger.error(f"Transfer {self.data.id} failed: {error}")
            raise IOError(f"Transfer {self.data.id} failed: {error}")
        rel_paths = abs_to_rel_paths(src, get_abs_files(src, manifest))

        num_files = len(rel_paths)
        parallelism = calculate_parallelism(num_files)
        logger.debug(
            f"Transfer {self.data.id}: Copy {num_files} files using {parallelism} threads"
        )
        parallel_copy_files_only(rel_paths, src, dest, parallelism=parallelism)

        logger.debug(f"Transfer {self.data.id}: Chmod files")
        file_mode = int(config.conf.get("SITE", "file_permissions"), base=8)
        folder_mode = int(config.conf.get("SITE", "folder_permissions"), base=8)
        parallel_chmod(dest, file_mode, folder_mode, parallelism)

        # temporarily add for testing
        logger.debug("Sleeping to test db connection")
        time.sleep(3700)


def check_transfer_queue(session) -> None:
    """
    Check the transfer queue and start the oldest transfer task, if any.  This only does one task because
    transfers typically take many minutes and the queue may change (e.g. a transfer is cancelled).
    """
    rows = []
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


def get_abs_files(root, rel_paths) -> list:
    """
    Create list of all source files by recursively expanding any folders.
    :param root: Folder under which all the paths in rel_paths exist
    :ptype root: str
    :param rel_paths: list of paths (files/folders) under root dir
    :ptype rel_paths:
    :return: Any folders in the input list shall be replaced with all the files under the folder.
    :rtype: list
    """
    abs_files = set()
    for item in rel_paths:
        full_path = os.path.normpath(os.path.join(root, item))
        if os.path.isfile(full_path):
            abs_files.add(full_path)
        elif os.path.isdir(full_path):
            abs_files.update(list_all_files_under_dir(full_path))
        # else symlink -- we don't use symlinks so they are not supported
    return list(abs_files)


def list_all_files_under_dir(root) -> list:
    """
    Walk down directory tree and return a list of all files therein.
    :param root: starting folder
    :ptype root: str
    :return: absolute paths of every file contained under root dir
    :rtype: list
    """
    files = []
    for path, dirnames, filenames in os.walk(root):
        for file in filenames:
            files.append(os.path.join(path, file))
    return files


def abs_to_rel_paths(root: str, paths: list) -> list:
    """
    Convert a list of paths from absolute to relative, given root dir.
    :param root: All paths shall be relative to this root folder.
    :ptype root: str
    :param paths: List of pathnames; all must be under root.
    :ptype paths: list
    :return: List of relative paths
    :rtype: list
    """
    rel_paths = []
    for abs_path in paths:
        rel_paths.append(os.path.relpath(abs_path, start=root))
    return rel_paths


def calculate_parallelism(num_files):
    """
    Calculate the amount of parallelism needed by the parallel_sync process.
    Each 10k files will require 1 thread with a max number of threads limited at 7.
    """
    if num_files < 0:
        raise ValueError("num_files cannot be negative")

    max_threads = int(config.conf.get("SITE", "max_transfer_threads"))

    if max_threads < 0:
        raise ValueError("max_threads must be greater than zero")

    upper_limit_files = max_threads * FILES_PER_THREAD
    min_threads = 1

    if num_files >= upper_limit_files:
        return max_threads
    parallelism = num_files // FILES_PER_THREAD
    return max(parallelism, min_threads)


def parallel_copy_files_only(manifest: list, src: str, dest: str, **kwargs):
    """
    Given list of files, copy them in parallel using parallel_sync.  Copies regular files only, skips others.
    :param manifest: list of file relative paths
    :ptype manifest: list
    :param src: source root directory
    :ptype src: str
    :param dest: destination root directory
    :ptype dest: str
    """
    parallelism = kwargs.get("parallelism", 1000)
    paths = []
    for rel_path in manifest:
        s = os.path.join(src, rel_path)
        if os.path.exists(s) and os.path.isfile(s):
            d = os.path.join(dest, rel_path)
            paths.append((s, d))
    rsync.local_copy(paths, parallelism=parallelism, extract=False, validate=False)


def parallel_chmod(path, file_mode, folder_mode, parallelism=3, **kwargs):
    """
    Recursively copy folder and set permissions.
    """
    if not os.path.isdir(path):
        if kwargs.get("ok_not_exists", False) is True:
            return
        else:
            raise IOError(f"Cannot chmod; path does not exist: {path}")
    if kwargs.get("chmod_parent", False) is True:
        try:
            parent = os.path.dirname(path)
            os.chmod(parent, folder_mode)
        except Exception as error:
            logger.warning(f"Error changing permissions of {parent}: {error}")
    try:
        os.chmod(path, folder_mode)
    except Exception as error:
        logger.warning(f"Error changing permissions of {path}: {error}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=parallelism) as executor:
        root_dir = os.path.abspath(path)
        for src_dir, dirs, files in os.walk(root_dir):
            for subdir in dirs:
                subdir_path = os.path.join(src_dir, subdir)
                try:
                    executor.submit(os.chmod, subdir_path, folder_mode)
                except Exception as e:
                    logger.warning(f"Error changing permissions of {subdir_path}: {e}")
            for file in files:
                file_path = os.path.join(src_dir, file)
                try:
                    executor.submit(os.chmod, file_path, file_mode)
                except Exception as e:
                    logger.warning(f"Error changing permissions of {file_path}: {e}")


def reset_queue(session):
    rows = []
    try:
        rows = (
            session.query(models.Transfer)
            .filter(models.Transfer.status == "transferring")
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
    for row in rows:
        transfer = Transfer(session, row)
        transfer.update_status("queued")


class FixPermsError(Exception):
    # base class for all errors in this package
    pass


class FixPermsDbError(FixPermsError):
    pass


class FixPermsNotFoundError(FixPermsError):
    pass


class FixPermsValueError(FixPermsError):
    pass


class FixPerms:
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
        try:
            data = models.Fix_Perms(
                base_dir=params["base_dir"],
            )
        except SQLAlchemyError as error:
            raise FixPermsValueError(
                f"Error creating model for new FixPerms: {params}: {error}"
            )
        try:
            session.add(data)
            session.commit()
        except SQLAlchemyError as error:
            session.rollback()
            raise FixPermsDbError(error)
        else:
            return cls(session, data)

    @classmethod
    def from_id(cls, session, fix_perms_id: int):
        """Select existing FixPerms record from RDb by primary key"""
        try:
            data = (
                session.query(models.Fix_Perms)
                .filter_by(id=int(fix_perms_id))
                .one_or_none()
            )
        except SQLAlchemyError as error:
            raise FixPermsDbError(f"Error selecting FixPerms {fix_perms_id}: {error}")
        except Exception as error:
            raise FixPermsError(f"Error selecting FixPerms {fix_perms_id}: {error}")
        else:
            if data is None:
                raise FixPermsNotFoundError(f"FixPerms {fix_perms_id} not found")
            else:
                return cls(session, data)

    def fix_perms(self, parallelism=3):
        """
        Recursively change the permissions of folders and files.
        """
        file_mode = int(config.conf.get("SITE", "file_permissions"), base=8)
        folder_mode = int(config.conf.get("SITE", "folder_permissions"), base=8)
        try:
            parallel_chmod(self.data.base_dir, file_mode, folder_mode, parallelism)
        except Exception as error:
            logger.error(f"Fix perms {self.data.id} failed: {error}")
            self.update_status("failed", str(error))
        else:
            self.udpate("succeeded")

    def update_status(self, new_status: str, reason: str = None) -> None:
        """
        Update status.
        """
        logger.info(f"Fix Perms {self.data.id}: now {new_status}")
        if reason is not None:
            reason = reason[:MAX_ERROR_STRING_LEN]
        timestamp = datetime.utcnow()
        try:
            self.data.status = new_status
            self.data.updated = timestamp
            if reason is not None:
                self.data.reason = reason
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Fix Perms {self.data.id}: {error}")


def check_fix_perms_queue(session) -> None:
    """
    Do any chmod tasks for Globus transfers.
    """
    rows = []
    try:
        rows = (
            session.query(models.Fix_Perms)
            .filter(models.Fix_Perms.status == "queued")
            .order_by(models.Fix_Perms.id)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
    for row in rows:
        fix_perms = FixPerms(session, row)
        fix_perms.fix_perms()
