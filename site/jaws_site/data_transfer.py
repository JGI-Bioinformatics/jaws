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

    def submit_transfer(self, manifest_files: list, **kwargs: dict) -> str:
        ...

    def transfer_status(self, transfer_task_id: str) -> str:
        ...

    def cancel_transfer(self, transfer_task_id: str):
        ...


class DataTransferFactory:
    """Factory class with plugin architecture to return a data transfer object. Based on the input object type
    which is the name of the module sans .py within the data_transfer_plugins directory, the factory dynamically
    imports the module and returns an instantiated object of the DataTransfer class from the imported module.

    Example of input object types:
      "globus_transfer" -> data_transfer_plugins/globus_transfer.py
      "aws_s3_transfer" -> data_transfer_plugins/aws_s3_transfer.py
    """

    def __new__(cls, obj_type: str) -> Callable:
        classname = "DataTransfer"
        modulename = f"jaws_central.data_transfer_plugins.{obj_type}"

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
