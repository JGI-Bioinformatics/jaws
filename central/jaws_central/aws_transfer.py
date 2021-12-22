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

logger = logging.getLogger(__package__)


class DataTransferFactory:
    def __init__(self) -> None:
        self._config = config.Configuration()

    def submit_transfer(self, label, src_site_id, dest_site_id, manifest_file) -> Dict:
        """
        label : human readable label (e.g. "Upload Run 2552") -- optional
        src_site_id (e.g. "cori", "jgi", "tahoma")
        dest_site_id (e.g. "aws")
        manifest_file (list of paths)
        kwargs
        "Save the transfer in the queue (database) and return transfer ID (pk)"
        -> returns dictionary of { input : URI }
        """
        pass

    def _add_transfer(self) -> str:
        """Insert transfer into table with "queued" initial state"""
        pass

    def _submit_transfer(self, transfer_id):
        """Submit the transfer and wait until done"""
        pass

    def complete_transfer(self):
        """Called via callback (from REST server or via RPC from REST server),
        update row to change state to "completed" or "failed"""
        pass

    def transfer_status(self, transfer_id):
        """Query db and return current state"""
        pass

    def cancel(self):
        """Cancel transfer (optional in first version)"""
        pass


class DataTransferS3(DataTransferFactory):
    def __init__(self) -> None:
        # Init the DataTransfer  Factory to get self._config
        super().__init__(self)
        logger.debug("Creating DataTransferS3")

        self.host_path = self._config.get("AWS", "host_path")
        self.aws_access_key_id = self._config.get("AWS", "aws_access_key_id")
        self.aws_secret_access_key = self._config.get("AWS", "aws_secret_access_key")

        # Setup aws clinet based on site configuration
        _aws_client = boto3.client(
            "s3",
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
        )
        bucket_name = self._config.get("AWS", "s3_bucket")
        self.s3_bucket = _aws_client.Bucket(bucket_name)
        self._transfer_threads = dict()

    def _add_transfer(self) -> str:
        """
        Return an id associated with transfer

        Create a random uuid for aws transfer
            -> Could also send to a database
        """
        return str(uuid.uuid4())

    def submit_transfer(self, label, manifest_file) -> Dict:
        """
        Submit a transfer to S3

        :param label: label for the data transfer
        :param host_paths: a dictionary that includes source host path and compute host path. This is used to modify
        the paths and create virtual transfer paths.
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
            self._submit_transfer(transfer_id, src, dest, label=label)

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
            if inode_type == "D":
                # If it's a directory get all the files and transfer them
                files = glob(f"{source_path}/*")
                for fil in files:
                    source_paths.append(f"{source_path}/{fil}")
                    dest_paths.append(f"{dest_path}/{fil}")
            else:
                source_paths.append(f"{source_path}/{fil}")
                dest_paths.append(f"{dest_path}/{fil}")

        return source_paths, dest_paths

    def _submit_transfer(self, transfer_id, source_path, dest_path, label=None):
        # Get the filename from the src_path
        filename = os.path.basename(source_path)

        # Object name inside of S3 bucket
        # Can include a folder name inside bucket as dest_dir
        object_name = f"{dest_path}/{filename}"

        logger.debug(f"add transfer: {source_path} -> {object_name}")
        extra_args = {"Metadata": {"label": label}}

        # Create a thread to transfer data
        #   Could run into issue with lots of threads transfering data
        #   Look into making a queue to start transfers
        self._transfer_threads[transfer_id] = threading.Thread(
            target=self._aws_client.upload_file,
            args=(source_path, self.s3_bucket, object_name),
            kwargs={"ExtraArgs": extra_args},
        )

        # Start thread running and return
        self._transfer_thread[transfer_id].start()
        return True

    def transfer_status(self, transfer_id):
        """
        Return the status of the transfer
        """
        try:
            # While transfering thread "is_alive" -> True
            # When done transfering thread "is_alive" -> False
            alive = self._transfer_thread[transfer_id].is_alive()
        except KeyError:
            # If the uuid is out of the
            logging.error(f"UUID: {transfer_id} not in table")
            return "failed"
        except Exception as e:
            logging.error(e)
            return "failed"

        if not alive:
            try:
                self._transfer_thread[transfer_id].join()
                self._transfer_thread.pop(transfer_id, None)
                return "upload complete"
            except ClientError as e:
                logging.error(e)
                return "failed"
        else:
            return "uploading"
