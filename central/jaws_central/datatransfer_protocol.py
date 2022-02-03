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
    def submit_upload(self, metadata: str, manifest_file: list):
        ...

    def submit_download(self, metadata: dict, src_dir: str, dst_dir: str):
        ...

    def transfer_status(self, task_id: str):
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
        'AWS': 'aws_transfer',
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
        modulename = f"datatransfer_plugins.{obj_type}"

        try:
            imported_module = importlib.import_module(modulename)
        except ImportError:
            raise DataTransferError(f"No module {modulename} found")

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

# Transfer arguments as *args, **kwargs
def transfer_file(obj: DataTransferProtocol, *args, **kwargs) -> List[str]:
    return obj.submit_transfer(*args, **kwargs)

# Or directly specifying them
def transfer_status(obj: DataTransferProtocol, transfer_id: str) -> str:
    return obj.transfer_status(transfer_id)

def main():
    try:
        # call based on filename/modulename in datatransfer_plugins
        obj = DataTransferFactory("aws_transfer")
    except DataTransferError as err:
        print(err)
        raise SystemExit()

    tansfer_ids = transfer_file(obj, label, src_site_id, dest_site_id, manifest_file)
    print(transfer_status(obj, tansfer_ids[0]))
"""
