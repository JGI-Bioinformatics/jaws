import logging
import os
from typing import List, Tuple
from jaws_central import config
import uuid
import boto3
from botocore.exceptions import ClientError
from ..datatransfer_protocol import SiteTransfer

logger = logging.getLogger(__package__)


class DataTransfer:
    def __init__(self, **kwargs) -> None:
        """DataTransferS3 is based on DataTransferProtocol"""

        if "aws_access_key_id" not in kwargs:
            raise (ValueError("aws_access_key_id required"))
        if "aws_secret_access_key" not in kwargs:
            raise (ValueError("aws_secret_access_key required"))
        if "s3_bucket" not in kwargs:
            raise (ValueError("s3_bucket required"))
        self.aws_access_key_id = kwargs["aws_access_key_id"]
        self.aws_secret_access_key = kwargs["aws_secret_access_key"]
        self.s3_bucket = kwargs["s3_bucket"]

        # Setup aws clinet based on site configuration
        self._aws_client = boto3.client(
            "s3",
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
        )

    def submit_transfer(self, manifest: list, **kwargs: dict) -> str:
        """
        Submit a transfer between local file system and S3.
        :param manifest: table of files/folders to transfer
        :type manifest: list
        :param kwargs: src/dest endpoints, base paths, and optional label
        :type kwargs: dict
        :return: transfer ID
        :rtype: str
        """
        label = kwargs.get("label", f"Transfer {len(manifest)} files/folders")
        logger.debug(f"S3 Transfer starting for {label}")
        transfer_id = str(uuid.uuid4())

        # Get source and destinations from manifest
        source_paths, dest_paths = self._get_manifest_paths(manifest)

        # For each file call _submit_transfer
        for src, dest in zip(source_paths, dest_paths):
            if kwargs["direction"] == "upload":
                try:
                    self._aws_client.upload_file(src, self.s3_bucket, dest)
                except Exception as error:
                    logger.error(f"Error uploading {src}", error)
            else:
                try:
                    self._aws_client.download_file(src, self.s3_bucket, dest)
                except Exception as error:
                    logger.error(f"Error uploading {src}", error)

        return transfer_id

    def cancel_transfer(self, task_id: str) -> None:
        """TODO: Need to implement. As of current 2/2022, boto3 doesn't support cancelling s3 transfers."""
        pass

    def _get_manifest_paths(self, manifest: list) -> Tuple[List[str], List[str]]:
        """
        Given a manifest of a list of:
            source_path, dest_path, inode_type
        return these modifed for transfering to S3

        :param manifest: manifest of all the files to be transferred
        :return: list of source and destination file paths
        :rtype: list
        """
        source_paths = []
        dest_paths = []
        for line in manifest:
            line = line.decode("UTF-8")
            source_path, dest_path, inode_type = line.split("\t")

            if inode_type == "D" or os.path.isdir(source_path):
                source_paths, dest_paths = self._get_file_paths(source_path, dest_path)
            else:
                source_paths.append(f"{source_path}")
                dest_paths.append(f"{dest_path}")

        return source_paths, dest_paths

    def _get_file_paths(
        self, src_path: str, dest_path: str
    ) -> Tuple[List[str], List[str]]:
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

    def transfer_status(self, transfer_id: str) -> str:
        """
        Transfer status is not applicable to this class.

        :param transfer_id: ID of the tansfer we want to check on.
        """
        return "succeeded"
