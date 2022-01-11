from typing import Protocol, Dict


class DataTransferException(Exception):
    pass


class DataTransferFactory(Protocol):
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
