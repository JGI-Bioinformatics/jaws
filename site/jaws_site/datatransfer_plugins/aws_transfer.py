import logging
import os
from typing import List, Tuple
from jaws_site import config
from glob import glob
import threading

# For creating random uuid
import uuid

# For aws transfers we'll use boto3
import boto3
from botocore.exceptions import ClientError

from ..datatransfer_protocol import SiteTransfer

logger = logging.getLogger(__package__)


class DataTransfer:
    def __init__(self) -> None:
        """DataTransferS3 is based on DataTransferProtocol"""
        # Init the DataTransfer Factory to get self._config
        self._config = config.Configuration()
        logger.debug("Creating DataTransferS3")

        self.aws_access_key_id = self._config.get("AWS", "aws_access_key_id")
        self.aws_secret_access_key = self._config.get("AWS", "aws_secret_access_key")

        # Setup aws clinet based on site configuration
        self._aws_client = boto3.client(
            "s3",
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
        )
        self.s3_bucket = self._config.get("AWS", "s3_bucket")
        self._transfer_threads = dict()

    def _add_transfer(self) -> str:
        """
        Return an id associated with transfer

        Create a random uuid for aws transfer
            -> Could also send to a database
        """
        return str(uuid.uuid4())

    def submit_upload(self, metadata: dict, manifest_files: list) -> List[str]:
        """
        Submit a transfer to S3

        :param metadata: dict containing entries for data transfer, i.e. label.
        :type metadata: dict
        :param manifest_files: list of files to transfer
        :type manifest_files: list
        :return list: Return list of transfer ID's
        """

        label = metadata.get('label', '')
        transfer_ids = []

        # Get source and destinations from manifest
        source_paths, dest_paths = self._get_manifest_paths(manifest_files)

        # For each file call _submit_transfer
        for src, dest in zip(source_paths, dest_paths):
            transfer_id = self._add_transfer()
            logger.debug(f"S3 Transfer starting for {transfer_id} {label}")
            transfer_ids.append(transfer_id)
            # If we're submitting to aws we need to upload data
            self._submit_upload(transfer_id, src, dest, label=label)

        # The client code expects a transfer id string to be returned, not a list. Hence, we're returning 0 in
        # this case instead of a list of transfer ids.

        # return transfer_ids
        return 0

    def submit_download(self, metadata: dict, src_path: str, dest_path: str) -> List[str]:
        """
        Submit a transfer to S3

        :param metadata: dict containing entries for data transfer, i.e. label.
        :type metadata: dict
        :param src_path: path containing the source files to transfer
        :param dest_path: path containing the destnation files to tranfser
        :return list: Return list of transfer ID's
        """

        label = metadata.get('label', '')
        transfer_ids = []

        # Get source and destinations from manifest
        source_paths, dest_paths = self._get_file_paths(src_path, dest_path)

        # For each file call _submit_transfer
        for src, dest in zip(source_paths, dest_paths):
            transfer_id = self._add_transfer()
            logger.debug(f"S3 Transfer starting for {transfer_id} {label}")
            transfer_ids.append(transfer_id)
            # If we're submitting to aws we need to upload data
            self._submit_download(transfer_id, src, dest, label=label)

        # client is expecting one transfer id, not a list so here, we don't return anything
        # return transfer_ids

    def cancel_transfer(self, task_id: str) -> None:
        """TODO: Need to implement. As of current 2/2022, boto3 doesn't support cancelling s3 transfers."""
        pass

    def _get_manifest_paths(self, manifest_files: list) -> Tuple[List[str], List[str]]:
        """
        Given a manifest of a list of:
            source_path, dest_path, inode_type
        return these modifed for transfering to S3

        :param manifest_files: manifest of all the files to be transferred
        :return: list of source and destination file paths
        :rtype: list
        """
        source_paths = []
        dest_paths = []
        for line in manifest_files:
            line = line.decode("UTF-8")
            source_path, dest_path, inode_type = line.split("\t")

            if inode_type == "D" or os.path.isdir(source_path):
                # If it's a directory get all the files and transfer them
                files = glob(f"{source_path}/*")
                for fil in files:
                    _fil = os.path.basename(fil)
                    source_paths.append(f"{source_path}/{_fil}")
                    dest_paths.append(f"{dest_path}")
            else:
                source_paths.append(f"{source_path}")
                dest_paths.append(f"{dest_path}")

        return source_paths, dest_paths

    def _get_file_paths(self, src_path: str, dest_path: str) -> Tuple[List[str], List[str]]:
        """
        Given a source file path and a destination file path, construct a list of source files and destination
        files recursively.

        :param src_path: source file path
        :type src_path: str
        :param dst_path: source file path
        :type dst_path: str
        :return: list of source and destination files
        :rtype: tuple of two lists
        """
        src_files = []
        dest_files = []
        for root, dirs, files in os.walk(src_path):
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_path, os.path.basename(src_file))
                src_files.append(src_file)
                dest_files.append(dest_file)

        return src_files, dest_files

    def _submit_upload(self, transfer_id: str, source_path: str, dest_path: str, label: str = None) -> bool:
        """PRIVATE, used to upload a file to aws"""
        # Get the filename from the src_path
        filename = os.path.basename(source_path)

        # Object name inside of S3 bucket
        # Can include a folder name inside bucket as dest_dir
        object_name = f"{dest_path}/{filename}"

        logger.debug(f"add transfer: {source_path} -> {object_name}")
        extra_args = {"Metadata": {"label": label}}
        try:
            # Create a thread to transfer data
            #   Could run into issue with lots of threads transfering data
            #   Look into making a queue to start transfers
            self._transfer_threads[transfer_id] = threading.Thread(
                target=self._aws_client.upload_file,
                args=(source_path, self.s3_bucket, object_name),
                kwargs={"ExtraArgs": extra_args},
            )

            # Start thread running and return
            self._transfer_threads[transfer_id].start()
            return True
        except RuntimeError as e:
            logger.error("Error uploading to aws", e)
            return False

    def _submit_download(self, transfer_id: str, source_path: str, dest_path: str, label=None) -> bool:
        """PRIVATE, used to download a file from aws"""
        # Get the filename from the src_path
        filename = os.path.basename(source_path)

        # Object name inside of S3 bucket
        # Can include a folder name inside bucket as dest_dir
        object_name = f"{dest_path}/{filename}"

        logger.debug(f"add transfer: {source_path} -> {object_name}")
        extra_args = {"Metadata": {"label": label}}
        try:
            # Create a thread to transfer data
            #   Could run into issue with lots of threads transfering data
            #   Look into making a queue to start transfers
            self._transfer_threads[transfer_id] = threading.Thread(
                target=self._aws_client.download_file,
                args=(source_path, self.s3_bucket, dest_path),
                kwargs={"ExtraArgs": extra_args},
            )

            # Start thread running and return
            self._transfer_threads[transfer_id].start()
            return True
        except RuntimeError as e:
            logger.error("Error uploading to aws", e)
            return False

    def transfer_status(self, transfer_id: str) -> str:
        """
        Return the status of the transfer

        :param transfer_id: ID of the tansfer we want to check on.
        """
        # If we already checked the status it is stored as the transfer ID
        # Just return the string value
        if isinstance(self._transfer_threads[transfer_id], str):
            return self._transfer_threads[transfer_id]

        try:
            # While transfering thread "is_alive" -> True
            # When done transfering thread "is_alive" -> False
            alive = self._transfer_threads[transfer_id].is_alive()
        except KeyError:
            # If the uuid is not in the thread table anymore
            logging.error(f"UUID: {transfer_id} not in table")
            self._transfer_threads[transfer_id] = "failed"
            return SiteTransfer.status.failed
        except Exception as e:
            self._transfer_threads[transfer_id] = "failed"
            logging.error(e)
            return SiteTransfer.status.failed

        if not alive:
            try:
                self._transfer_threads[transfer_id].join()
                self._transfer_threads[transfer_id] = "upload complete"
                return SiteTransfer.status.succeeded
            except ClientError as e:
                logging.error(e)
                self._transfer_threads[transfer_id] = "failed"
                return SiteTransfer.status.failed
        else:
            return SiteTransfer.status.transferring
