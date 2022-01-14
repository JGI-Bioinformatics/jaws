import logging
import os
from typing import Dict
from jaws_central import config
from glob import glob
import threading

# For creating random uuid
import uuid

# For aws transfers we'll use boto3
import boto3
from botocore.exceptions import ClientError

from ..datatransfer_protocol import DataTransferError

logger = logging.getLogger(__package__)


class DataTransfer:
    def __init__(self) -> None:
        """DataTransferS3 is based on DataTransferProtocol"""
        # Init the DataTransfer Factory to get self._config
        self._config = config.Configuration()
        logger.debug("Creating DataTransferS3")

        # self.host_path = self._config.get("AWS", "host_path")
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

    def submit_transfer(self, label, src_site_id, dest_site_id, manifest_file) -> Dict:
        """
        Submit a transfer to S3

        :param label: label for the data transfer
        :param src_site_id: Source of the data currently. (cori, jgi, tahoma, aws)
        :param dest_site_id: Where the data will be transfered to. (cori, jgi, tahoma, aws)
        :param manifest_file: manifest of all the files to be transferred
        :return list: Return list of transfer ID's
        """

        transfer_ids = []

        # Get source and destinations from manifest
        source_paths, dest_paths = self._get_transfer_paths(manifest_file)

        # For each file call _submit_transfer
        for src, dest in zip(source_paths, dest_paths):
            transfer_id = self._add_transfer()
            logger.debug(f"S3 Transfer starting for {transfer_id} {label}")
            transfer_ids.append(transfer_id)
            # If we're submitting to aws we need to upload data
            if dest_site_id == "aws":
                self._submit_upload(transfer_id, src, dest, label=label)
            # If we're gettting from aws we need to download
            elif src_site_id == "aws":
                self._submit_download(transfer_id, src, dest, label=label)
            else:
                raise DataTransferError("Not an AWS transfer")

        return transfer_ids

    def _get_transfer_paths(self, manifest_file):
        """
        Given a manifest of a list of:
            source_path, dest_path, inode_type
        return these modifed for transfering to S3

        :param manifest_file: manifest of all the files to be transferred
        """
        source_paths = []
        dest_paths = []
        for line in manifest_file:
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

    def _submit_upload(self, transfer_id, source_path, dest_path, label=None):
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

    def _submit_download(self, transfer_id, source_path, dest_path, label=None):
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

    def transfer_status(self, transfer_id):
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
            return "failed"
        except Exception as e:
            self._transfer_threads[transfer_id] = "failed"
            logging.error(e)
            return "failed"

        if not alive:
            try:
                self._transfer_threads[transfer_id].join()
                self._transfer_threads[transfer_id] = "upload complete"
                return "upload complete"
            except ClientError as e:
                logging.error(e)
                self._transfer_threads[transfer_id] = "failed"
                return "failed"
        else:
            return "uploading"
