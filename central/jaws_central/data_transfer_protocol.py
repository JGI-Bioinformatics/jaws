from typing import Protocol, Dict


class DataTransferException(Exception):
    pass


class DataTransferProtocol(Protocol):
    def submit_transfer(self, label, src_site_id, dest_site_id, manifest_file) -> Dict:
        """
        Save the transfer in the queue (database) and return transfer ID (pk)
        """

    def _add_transfer(self) -> str:
        """Insert transfer into table with "queued" initial state"""

    def _submit_transfer(self, transfer_id):
        """Submit the transfer and wait until done"""

    def complete_transfer(self):
        """
        Called via callback (from REST server or via RPC from REST server),
        update row to change state to 'completed' or 'failed'
        """

    def transfer_status(self, transfer_id):
        """Query db and return current state"""

    def cancel(self):
        """Cancel transfer (optional in first version)"""


class DataTransferFactory:
    def start_transfer(self, obj: DataTransferProtocol) -> None:
        return obj.submit_transfer()
