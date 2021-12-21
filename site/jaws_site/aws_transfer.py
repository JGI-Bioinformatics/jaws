import logging
import os
from typing import Dict
import jaws_site.config

# For aws transfers we'll use boto3
import boto3
from botocore.exceptions import ClientError

logger = logging.getLogger(__package__)


class DataTransferFactory:
    def __init__(self) -> None:
        self._config = jaws_site.config.Configuration()
        pass

    def submit_transfer(self) -> Dict:
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

    def _add_transfer(self):
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

        self.host_path = self._config.get("AWS", "host_path")
        self.s3_bucket = self._config.get("AWS", "s3_bucket")

        self.aws_access_key_id = self._config.get("AWS", "aws_access_key_id")
        self.aws_secret_access_key = self._config.get("AWS", "aws_secret_access_key")

        # Setup aws clinet based on site configuration
        self._aws_client = boto3.client(
            "s3",
            aws_access_key_id=self.aws_access_key_id,
            aws_secret_access_key=self.aws_secret_access_key,
        )

    def submit_transfer(self, label, dest_endpoint, src_dir, dest_dir) -> Dict:
        self._add_transfer()
        self._submit_transfer(transfer_id="")

    def _submit_transfer(self, transfer_id):
        self._aws_client.upload_file(
            "FILE_NAME",
            self.s3_bucket,
            "OBJECT_NAME",
            ExtraArgs={"Metadata": {"mykey": "myvalue"}},
        )
