"""
Transfer is a collection of files/folders to transfer (e.g. via Globus, FTP, etc.)
Items are stored in a relational database.
"""
import concurrent.futures
import logging
import os
from datetime import datetime
from sqlalchemy.exc import SQLAlchemyError
import json
import boto3
from jaws_site import config, models
from jaws_site.globus import GlobusService
import botocore
from parallel_sync import rsync


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


class TransferGlobusError(TransferError):
    pass


class TransferNotFoundError(TransferError):
    pass


class TransferKeyError(TransferError):
    pass


class TransferValueError(TransferError):
    pass


class Transfer:
    """
    Class representing a transfer (set of files to transfer) associated with one Run.
    States are:
    - queued
    - active
    - failed
    - succeeded
    - done
    - cancelled
    """

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
        self.operations = {
            "queued": self.begin_transfer,
            "active": self.check_globus_status,
            "succeeded": self.change_permissions,
            "failed": self.change_permissions,
        }

    @classmethod
    def from_params(cls, session, params):
        """Create new transfer from parameter values and save in RDb."""
        # manifest of files to transfer is saved as json string that must be decoded before use
        manifest_json = "[]"
        if "manifest" in params:
            assert isinstance(params["manifest"], list)
            manifest_json = json.dumps(params["manifest"])
        elif "manifest_json" in params:
            assert isinstance(params["manifest_json"], str)
            manifest_json = params["manifest_json"]

        # determine if the transfer type is valid
        transfer_id = params.get("transfer_id", 0)
        transfer_type = params.get("transfer_type", None)
        if transfer_type is None:
            err = f"Transfer {transfer_id}: Missing transfer type"
            logger.error(err)
            raise TransferKeyError(err)
        if transfer_type not in ["local", "globus", "s3"]:
            err = f"Transfer {transfer_id}: Unsupported type, {transfer_type}"
            logger.error(err)
            raise TransferValueError(err)

        # verify required Globus params are defined, if applicable
        if transfer_type == "globus":
            required_params = [
                "src_globus_endpoint",
                "src_globus_host_path",
                "dest_globus_endpoint",
                "dest_globus_host_path",
            ]
            for key in required_params:
                value = params.get(key, None)
                if value is None:
                    raise TransferKeyError(f"Missing required parameter: {key}")

        # verify datatypes
        if not isinstance(params["transfer_id"], int):
            raise TransferValueError(f"Transfer {transfer_id}: Invalid id")
        if not isinstance(params["src_base_dir"], str):
            raise TransferValueError(f"Transfer {transfer_id}: Invalid src_base_dir")
        if not isinstance(params["dest_base_dir"], str):
            raise TransferValueError(f"Transfer {transfer_id}: Invalid dest_base_dir")

        # insert row
        try:
            data = models.Transfer(
                id=transfer_id,
                transfer_type=transfer_type,
                status="queued",
                src_base_dir=params["src_base_dir"],
                dest_base_dir=params["dest_base_dir"],
                src_globus_endpoint=params.get("src_globus_endpoint", None),
                dest_globus_endpoint=params.get("dest_globus_endpoint", None),
                src_globus_host_path=params.get("src_globus_host_path", None),
                dest_globus_host_path=params.get("dest_globus_host_path", None),
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
        """Return key fields about the status of the transfer."""
        info = {
            "status": self.data.status,
            "reason": self.data.reason,
            "result": self.data.result,
        }
        return info

    def check_globus_status(self):
        """
        If transfer is a Globus transfer, query service for current status.  Otherwise, do nothing.
        :return: new status
        :rtype: str
        """
        if self.data.transfer_type != "globus":
            return
        globus_client = GlobusService()
        try:
            new_status, reason = globus_client.transfer_status(
                self.data.globus_transfer_id
            )
            new_status = new_status.lower()
        except Exception as error:
            logger.warn(
                f"Transfer {self.data.id}: Unable to check Globus status: {error}"
            )
            return
        else:
            if new_status != self.data.status:
                logger.debug(f"Transfer {self.data.id} status = {new_status}")
                # this will also update the "result" field, if applicable
                self.update_status(new_status)

    def check_status(self) -> None:
        """Check the status of the transfer and promote to next state if ready."""
        logger.debug(f"Transfer {self.data.id} is {self.data.status}")
        if self.data.status in self.operations:
            self.operations[self.data.status]()

    def manifest(self) -> list:
        """
        Return list of files to transfer -- may be empty list if complete folder is to be transferred.
        """
        return json.loads(self.data.manifest_json)

    def cancel(self) -> None:
        """
        Attempt to cancel a transfer.
        :return: Whether or not the transfer was cancelled.
        :rtype: bool
        """
        result = False
        if self.data.status == "queued":
            self.update_status("cancelled")
            result = True
        elif self.data.status == "active" and self.data.transfer_type == "globus":
            globus_client = GlobusService()
            result = globus_client.cancel_task(self.data.globus_transfer_id)
        if result is True:
            logger.debug(f"Transfer {self.data.id}: cancelled")
        return result

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
            if new_status in ["succeeded", "failed"]:
                self.data.result = new_status
            self.session.commit()
        except SQLAlchemyError as error:
            self.session.rollback()
            logger.exception(f"Unable to update Transfer {self.data.id}: {error}")
            raise (error)

    def begin_transfer(self) -> None:
        """
        Call appropriate transfer method.
        """
        if self.data.transfer_type == "globus":
            self.submit_globus_transfer()
        elif self.data.transfer_type == "s3":
            self.s3_copy()
        elif self.data.transfer_type == "local":
            self.local_copy()
        else:
            self.update_status("failed", "Invalid transfer type")

    def submit_globus_transfer(self) -> None:
        logger.debug(f"Submit Globus transfer {self.data.id}")
        label = f"Transfer {self.data.id}"
        manifest = self.manifest()
        try:
            globus_client = GlobusService()
            globus_transfer_id = globus_client.submit_transfer(
                label,
                self.data.src_globus_endpoint,
                self.data.src_globus_host_path,
                self.data.src_base_dir,
                self.data.dest_globus_endpoint,
                self.data.dest_globus_host_path,
                self.data.dest_base_dir,
                manifest,
            )
        except Exception as error:
            logger.error(f"Globus transfer {self.data.id} failed: {error}")
            self.update_status("failed", "Globus submission failed")
            raise TransferGlobusError(error)
        else:
            self.data.globus_transfer_id = globus_transfer_id
            self.update_status("active", f"globus_transfer_id={globus_transfer_id}")

    def s3_copy(self) -> None:
        """
        Do the S3 transfer and return when done -- the operation is blocking.
        """
        logger.debug(f"Begin S3 transfer {self.data.id}")
        self.update_status("active")
        if self.data.src_base_dir.startswith("s3://"):
            self.s3_download()
        else:
            self.s3_upload()

    def aws_s3_resource(self):
        aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
        aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
        aws_region_name = config.conf.get("AWS", "aws_region_name")
        try:
            aws_session = boto3.Session(
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key,
                region_name=aws_region_name,
            )
        except Exception as error:
            raise (f"Error getting aws session: {error}")
        try:
            s3_resource = aws_session.resource("s3")
        except Exception as error:
            raise (f"Error getting s3 sources: {error}")
        else:
            return s3_resource

    def aws_s3_client(self):
        try:
            aws_access_key_id = config.conf.get("AWS", "aws_access_key_id")
            aws_secret_access_key = config.conf.get("AWS", "aws_secret_access_key")
            client = boto3.client(
                "s3",
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key,
            )
        except Exception as error:
            raise (f"Failed to get AWS client: {error}")
        else:
            return client

    @staticmethod
    def s3_parse_path(full_path):
        full_path = full_path.replace("s3://", "", 1)  # discard S3 identifier
        folders = full_path.split("/")
        s3_bucket = folders.pop(0)
        path = "/".join(folders)
        # print(f"S3 BUCKET={s3_bucket}; PATH={path}")
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
        if aws_s3_client is None:
            aws_s3_client = self.aws_s3_client()
        try:
            result = aws_s3_client.list_objects_v2(Bucket=bucket, Prefix=file_key)
        except botocore.exceptions.ClientError as error:
            logger.error(f"Error getting S3 obj stats for {file_key}: {error}")
            raise
        except botocore.exceptions.ParamValidationError as error:
            raise ValueError(
                "The parameters you provided are incorrect: {}".format(error)
            )
        except Exception as error:
            logger.error(f"Error getting S3 obj stats for {file_key}: {error}")
            raise
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
        try:
            aws_client = self.aws_s3_client()
        except Exception as error:
            # do not fail, retry later
            logger.error(f"{error}")
            return
        try:
            aws_s3_resource = self.aws_s3_resource()
        except Exception as error:
            # do not fail, retry later
            logger.error(f"{error}")
            return
        try:
            s3_bucket, dest_base_dir = self.s3_parse_path(self.data.dest_base_dir)
        except Exception as error:
            msg = f"Error parsing s3 uri, {self.data.dest_base_dir}: {error}"
            logger.error(msg)
            self.update_status("failed", msg)
            return
        try:
            bucket_obj = aws_s3_resource.Bucket(s3_bucket)
        except Exception as error:
            msg = f"Error accessing bucket, {s3_bucket}: {error}"
            logger.error(msg)
            self.update_status("failed", msg)
            return
        for rel_path in manifest:
            src_path = os.path.normpath(os.path.join(self.data.src_base_dir, rel_path))
            src_file_size = os.path.getsize(src_path)
            dest_path = os.path.normpath(os.path.join(dest_base_dir, rel_path))
            try:
                dest_file_size = self.s3_file_size(s3_bucket, dest_path, aws_client)
            except Exception as error:
                self.update_status("failed", f"{error}")
                return
            if dest_file_size is not None and src_file_size == dest_file_size:
                logger.debug(f"S3 upload: Skipping cached file {dest_path}")
            else:
                logger.debug(f"S3 upload to {s3_bucket}: {src_path} -> {dest_path}")
                try:
                    with open(src_path, "rb") as fh:
                        bucket_obj.upload_fileobj(fh, dest_path)
                except Exception as error:
                    msg = f"Failed to upload to S3, file {src_path} -> {dest_path}: {error}"
                    logger.error(msg)
                    self.update_status("failed", msg)
                    return
        self.update_status("succeeded")

    def s3_download(self):
        return (
            self.s3_download_files()
            if len(self.manifest())
            else self.s3_download_folder()
        )

    def s3_download_files(self):
        manifest = self.manifest()
        num_files = len(manifest)
        logger.debug(f"Transfer {self.data.id} begin s3 download of {num_files} files")
        aws_s3_resource = self.aws_s3_resource()
        s3_bucket, src_base_dir = self.s3_parse_path(self.data.src_base_dir)
        try:
            bucket_obj = aws_s3_resource.Bucket(s3_bucket)
        except Exception as error:
            msg = f"Error accessing bucket, {s3_bucket}: {error}"
            logger.error(msg)
            self.update_status("failed", msg)
            return
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
                msg = f"Unable to make download dir, {dest_folder}: {error}"
                logger.error(msg)
                self.update_status("failed", msg)
            try:
                with open(dest_path, "wb") as fh:
                    bucket_obj.download_fileobj(src_path, fh)
            except Exception as error:
                msg = f"Failed to download s3 file, {src_path}: {error}"
                logger.error(msg)
                self.update_status("failed", msg)
                return
        self.update_status("succeeded")

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
                        self.update_status("failed", msg)
                        return
                else:
                    try:
                        basedir = os.path.dirname(dest_path)
                        mkdir(basedir)
                    except Exception as error:
                        msg = f"Unable to make download dir, {basedir}: {error}"
                        logger.error(msg)
                        self.update_status("failed", msg)
                        return
                    try:
                        aws_s3_client.download_file(s3_bucket, rel_path, dest_path)
                    except Exception as error:
                        msg = f"S3 download error, {rel_path}: {error}"
                        logger.error(msg)
                        self.update_status("failed", msg)
                        return
        self.update_status("succeeded")

    def local_copy(self) -> None:
        """
        Copy files and folders (recursively).
        """
        manifest = self.manifest()
        src = f"{self.data.src_base_dir}/"
        if not os.path.isdir(src):
            err = f"Source directory not found: {src}"
            logger.error(f"Transfer {self.data.id}: {err}")
            self.update_status("failed", err)
            return
        rel_paths = abs_to_rel_paths(src, get_abs_files(src, manifest))
        num_files = len(rel_paths)
        self.data.num_files = num_files
        parallelism = calculate_parallelism(num_files)
        logger.debug(
            f"Transfer {self.data.id}: Local copy {num_files} files using {parallelism} threads"
        )
        self.update_status("active")

        dest = f"{self.data.dest_base_dir}/"
        try:
            mkdir(dest)
        except IOError as error:
            err = f"Mkdir failed: {error}"
            logger.error(f"Transfer {self.data.id} failed: {err}")
            self.update_status("failed", err)
            return

        try:
            parallel_copy_files(rel_paths, src, dest, parallelism=parallelism)
        except Exception as error:
            err = f"Local copy failed: {error}"
            logger.error(f"Transfer {self.data.id}: {err}")
            self.update_status("failed", err)
        else:
            self.update_status("succeeded", f"{num_files} copied")

    def change_permissions(self) -> None:
        """
        After copying/downloading files, the permissions shall be changed.
        """
        logger.debug(f"Transfer {self.data.id}: Chmod files")
        file_mode = int(config.conf.get("SITE", "file_permissions"), base=8)
        folder_mode = int(config.conf.get("SITE", "folder_permissions"), base=8)
        dest = self.data.dest_base_dir
        num_files = self.data.num_files
        if num_files is None:
            manifest = self.manifest()
            num_files = len(get_abs_files(dest, manifest))
        parallelism = calculate_parallelism(num_files)
        parallel_chmod(dest, file_mode, folder_mode, parallelism=parallelism)
        self.update_status("done")


def check_active_transfers(session) -> None:
    """
    Get active transfers from db, check, and update their status.
    Transfers are queried in a particular order of states.
    """
    # Check all Globus transfers first since they are quick, require no resources.
    try:
        rows = (
            session.query(models.Transfer)
            .filter(models.Transfer.transfer_type == "globus")
            .filter(models.Transfer.status.in_(["queued", "active"]))
            .order_by(models.Transfer.id)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
        return
    else:
        for row in rows:
            transfer = Transfer(session, row)
            transfer.check_status()

    # Check for completed transfers requiring chmod
    try:
        rows = (
            session.query(models.Transfer)
            .filter(models.Transfer.status == "succeeded")  # also failed?
            .order_by(models.Transfer.id)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
        return
    else:
        for row in rows:
            transfer = Transfer(session, row)
            transfer.check_status()

    # Check for local copy tasks and do only ONE because copying typically takes
    # many minutes and the queue may change (e.g. a transfer is cancelled)
    try:
        rows = (
            session.query(models.Transfer)
            .filter(models.Transfer.transfer_type == "local")
            .filter(models.Transfer.status == "queued")
            .order_by(models.Transfer.id)
            .limit(1)
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
    else:
        if len(rows):
            transfer = Transfer(session, rows[0])
            transfer.check_status()


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


def parallel_copy_files(manifest: list, src: str, dest: str, **kwargs):
    """
    Given list of files, copy them in parallel using rsync.  Copies regular files only, skips others.
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


def parallel_chmod(path, file_mode, folder_mode, parallelism=1):
    """
    Recursively copy folder and set permissions.
    """
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
            .filter(models.Transfer.transfer_type == "local")
            .filter(models.Transfer.status == "active")
            .all()
        )
    except SQLAlchemyError as error:
        logger.warning(
            f"Failed to select transfer task from db: {error}", exc_info=True
        )
    for row in rows:
        transfer = Transfer(session, row)
        transfer.update_status("queued")
