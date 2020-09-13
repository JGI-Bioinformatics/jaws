"""
Run class
"""

import os
import logging
from datetime import datetime
from jaws_site.database import Session
from jaws_site import models
from jaws_site import config
from jaws_site.cromwell import Cromwell

# from jaws_rpc import rpc_client
from jaws_site.Globus import Globus


logger = logging.getLogger(__package__)
globus = Globus
conf = config.conf
cromwell = Cromwell(conf.get("CROMWELL", "url"))


class Run:
    def __init__(self, run_id=None, session=None):
        self.run_id = run_id
        self.session = session
        if run_id:
            if self.session is None:
                self.session = Session()
            self.run = self.session.query(models.Run).get(self.run_id)

    def check_status(self):
        if self.run.status == "uploading":
            self.check_if_upload_complete()
        elif self.run.status == "upload complete":
            self.submit_run()
        elif self.run.status == "submitted":
            self.check_run_cromwell_status()
        elif self.run.status == "queued":
            self.check_run_cromwell_status()
        elif self.run.status == "running":
            self.check_run_cromwell_status()
        elif self.run.status == "succeeded":
            self.transfer_results()
        elif self.run.status == "failed":
            self.transfer_results()
        elif self.run.status == "downloading":
            self.check_if_download_complete()

    def check_if_upload_complete(self):
        """
        Query Globus to see if transfer is complete
        """
        logger.debug(f"Run {self.run.id}: Check upload status")
        try:
            globus_status = globus.transfer_status(
                self.run.transfer_refresh_token, self.run.upload_task_id
            )
        except Exception as error:
            logger.exception(
                f"Failed to check upload {self.run.upload_task_id}: {error}"
            )
            return
        if globus_status == "FAILED":
            self.update_run_status("upload failed")
        elif globus_status == "INACTIVE":
            self.update_run_status(
                "upload inactive", "Your endpoint authorization has expired"
            )
        elif globus_status == "SUCCEEDED":
            self.update_run_status("upload complete")
            self.submit_run()

    def submit_run(self):
        """
        Submit a run to Cromwell.
        """
        logger.debug(f"Run {self.run.id}: Submit to Cromwell")

        # Validate input
        uploads_dir = os.path.join(
            conf.get("GLOBUS", "root_dir"), conf.get("SITE", "uploads_subdirectory")
        )
        file_path = os.path.join(uploads_dir, self.run.user_id, self.run.submission_id)
        wdl_file = file_path + ".wdl"
        json_file = file_path + ".json"
        zip_file = file_path + ".zip"  # might not exist
        if not os.path.exists(json_file):
            logger.warning(f"Missing inputs JSON for run {self.run.id}: {json_file}")
            self.update_run_status("missing input", "Missing inputs JSON file")
            return
        if not os.path.exists(wdl_file):
            logger.warning(f"Missing WDL for run {self.run.id}: {wdl_file}")
            self.update_run_status("missing input", "Missing WDL file")
            return
        if not os.path.exists(zip_file):
            zip_file = None

        # submit to Cromwell
        cromwell_run_id = cromwell.submit(wdl_file, json_file, zip_file)
        if cromwell_run_id:
            self.run.cromwell_run_id = cromwell_run_id
            self.session.commit()
            self.update_run_status(
                "submitted", f"cromwell_run_id={cromwell_run_id}"
            )
        else:
            self.update_run_status("submission failed")

    def check_run_cromwell_status(self):
        """
        Check Cromwell for the status of one Run.
        """
        logger.debug(f"Run {self.run.id}: Check Cromwell status")
        try:
            cromwell_status = cromwell.get_status(self.run.cromwell_run_id)
        except Exception as error:
            logger.error(
                f"Unable to check Cromwell status of Run {self.run.id}: {error}"
            )
            return  # try again next time
        logger.debug(f"Run {self.run.id}: Cromwell status is {cromwell_status}")
        if cromwell_status == "Running":
            if self.run.status == "submitted":
                try:
                    metadata = cromwell.get_metadata(self.run.cromwell_run_id)
                except Exception:
                    return
                if metadata is None:
                    return
                self.run.cromwell_workflow_dir = metadata.get("workflowRoot")
                self.session.commit()
                self.update_run_status("queued")
        elif cromwell_status == "Failed":
            self.update_run_status("failed")
            self.transfer_results()
        elif cromwell_status == "Succeeded":
            if self.run.status == "queued":
                self.update_run_status("running")
            if self.run.status == "running":
                self.update_run_status("succeeded")
            self.transfer_results()
        elif cromwell_status == "Aborted":
            if self.run.status == "queued" or self.run.status == "running":
                self.update_run_status("cancelled")

    def transfer_results(self):
        """
        Send run output via Globus
        """
        logger.debug(f"Run {self.run.id}: Download output")
        if self.run.cromwell_workflow_dir is None:
            try:
                metadata = cromwell.get_metadata(self.run.cromwell_run_id)
            except Exception:
                return
            if metadata is None:
                return
            self.run.cromwell_workflow_dir = metadata.get("workflowRoot")
            self.session.commit()
        try:
            task_id = globus.transfer_folder(
                self.run.transfer_refresh_token,
                self.run.cromwell_workflow_dir,
                self.run.output_endpoint,
                self.run.output_dir,
                f"Download run {self.run.id}",
            )
        except Exception as error:
            logger.warn(f"Run {self.run.id} download error: {error}")
            return
        self.run.download_task_id = task_id
        self.session.commit()
        self.update_run_status(
            "downloading", f"download_task_id={self.run.download_task_id}"
        )

    def check_if_download_complete(self):
        """
        If download is complete, change state.
        """
        logger.debug(f"Run {self.run.id}: Check download status")
        try:
            globus_status = globus.transfer_status(
                self.run.transfer_refresh_token, self.run.download_task_id
            )
        except Exception as error:
            logger.exception(
                f"Failed to check download {self.run.download_task_id}: {error}"
            )
            return
        if globus_status == "SUCCEEDED":
            self.update_run_status("download complete")
        elif globus_status == "FAILED":
            self.update_run_status("download failed")

    def update_run_status(self, new_status, reason=None):
        """
        Update Run's current status in 'runs' table and insert entry into 'run_logs' table.
        """
        logger.info(f"Run {self.run.id}: now {new_status}")
        status_from = self.run.status
        timestamp = datetime.utcnow()

        # update current status in "runs" table
        try:
            self.run.status = new_status
            self.run.updated = timestamp
            self.session.commit()
        except Exception as error:
            logger.exception(f"Unable to update Run {self.run.id}: {error}")

        # populate "reason" field
        if new_status == "submitted":
            reason = f"cromwell_run_id={self.run.cromwell_run_id}"
        elif new_status == "downloading":
            reason = f"download_task_id={self.run.download_task_id}"

        # save state transition in "run_logs" table
        try:
            log_entry = models.Run_Log(
                run_id=self.run.id,
                status_from=status_from,
                status_to=new_status,
                timestamp=timestamp,
                reason=reason,
            )
            self.session.add(log_entry)
            self.session.commit()
        except Exception as error:
            logger.exception(
                f"Failed to create run_log object for Run {self.run.id} : {new_status}, {reason}: {error}"
            )
        # notifying Central of state change is handled by daemon.send_run_status_logs
