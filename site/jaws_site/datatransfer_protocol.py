from typing import Callable, Protocol
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
