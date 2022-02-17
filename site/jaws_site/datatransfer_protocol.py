from typing import Callable, Protocol
from dataclasses import dataclass
import importlib


class DataTransferError(Exception):
    pass


class DataTransferAPIError(DataTransferError):
    pass


class DataTransferNetworkError(DataTransferError):
    pass


class DataTransferProtocol(Protocol):
    """Interface for the data transfer object using the Protocol structural subtyping class."""
    def submit_upload(self, metadata: dict, manifest_files: list) -> str:
        ...

    def submit_download(self, metadata: dict, src_dir: str, dst_dir: str) -> str:
        ...

    def transfer_status(self, task_id: str) -> str:
        ...

    def cancel_transfer(self, task_id: str):
        ...


@dataclass
class Status:
    """Data transfer statuses."""
    succeeded = 1
    failed = 2
    transferring = 3
    inactive = 4


@dataclass
class SiteTransfer:
    """Defines the data transfer object types for each site. The values are used as arguments to the
    DataTransferFactory class.
    """
    type = {
        'CORI': 'globus_transfer',
        'JGI': 'globus_transfer',
        'TAHOMA': 'globus_transfer',
        'AWS': 'aws_s3_transfer',
    }
    status = Status()


class DataTransferFactory:
    """Factory class with plugin architecture to return a data transfer object. Based on the input object type
    which is the name of the module sans .py within the datatransfer_plugins directory, the factory dynamically
    imports the module and returns an instantiated object of the DataTransfer class from the imported module.

    Example of input object types:
      "globus_transfer" -> datatransfer_plugins/globus_transfer.py
      "aws_transfer" -> datatransfer_plugins/aws_transfer.py
    """
    def __new__(cls, obj_type: str) -> Callable:
        classname = "DataTransfer"
        modulename = f"jaws_site.datatransfer_plugins.{obj_type}"

        try:
            imported_module = importlib.import_module(modulename)
        except ImportError:
            raise ImportError(f"No module {modulename} found")

        try:
            obj = getattr(imported_module, classname)
        except AttributeError as err:
            msg = f"{imported_module} cannot be found or has no class named {classname}: {err}"
            raise DataTransferError(msg)

        return obj()


# Notes on implementing
"""
from datatransfer_protocol import DataTransferProtocol, DataTransferError

# Setup methods that use a DataTransferProtocol

# File upload
def submit_upload(obj: DataTransferProtocol, metadata: dict, manifest_files: list ) -> List[str]:
    return obj.submit_transfer(metadata, manifest_files)

# File download
def submit_download(obj: DataTransferProtocol, metadata: dict, src_dir: str, dst_dir:str) -> List[str]:
    return obj.submit_download(metadata, src_dir, dst_dir)


# Get status
def transfer_status(obj: DataTransferProtocol, transfer_id: str) -> str:
    return obj.transfer_status(transfer_id)

def main():
    try:
        # call based on filename/modulename in datatransfer_plugins
        obj = DataTransferFactory("aws_transfer")
    except DataTransferError as err:
        print(err)
        raise SystemExit()

    metadata = {'label': 'this is a label'}

    # Upload files to remote site
    tansfer_id = submit_upload(obj, metadata, manifest_files)

    # Download files from remote site
    tansfer_id = submit_download(obj, metadata, src_dir, dst_dir)

    # Check status
    print(transfer_status(obj, tansfer_id))
"""
