from typing import Protocol
import importlib


class DataTransferError(Exception):
    pass


class DataTransferProtocol(Protocol):
    def submit_transfer(self):
        ...

    def transfer_status(self):
        ...

    def cancel(self):
        ...


class DataTransferFactory:
    def __new__(cls, obj_type: str):
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
