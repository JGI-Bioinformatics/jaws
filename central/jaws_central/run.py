import logging
from datetime import datetime, timedelta
from flask import abort, request
from sqlalchemy.exc import SQLAlchemyError, IntegrityError
import globus_sdk
import jaws_central.globus
from jaws_central import config
from jaws_central import jaws_constants
from jaws_central.run_log import RunLog
from jaws_rpc import rpc_index
from jaws_central.models_fsa import db, Run, User


logger = logging.getLogger(__package__)

globus = jaws_central.globus.GlobusService()

run_active_states = [
    "created",
    "uploading",
    "upload inactive",
    "upload complete",
    "submitted",
    "queued",
    "running",
    "succeeded",
    "ready",
    "downloading",
]


run_pre_cromwell_states = [
    "created",
    "uploading",
    "upload inactivte",
    "upload complete",
    "upload stalled",
    "upload failed",
]


class RunError(Exception):
    """Baseclass for all Run exceptions"""

    pass


class RunNotFoundError(RunError):
    """The specified Run does not exist."""

    pass


class RunAccessDenied(RunError):
    """The requesting user does not have access to the Run."""

    pass


class Run:
    def __init__(self, session, user, run_id=None, params={}):
        """Initialize Run object. Either loads existing or inserts new.

        :param session: db handle
        :type session: sqlalchemy.session
        :param user: Current user object
        :type user: jaws_central.user.User
        :param run_id: Run's primary key; required if exists
        :type run_id: int
        :param params: parameters for new Run; required if not exists
        :type params: dict
        """
        self.session = session
        self.user = user
        self.id = run_id
        if run.id:
            self._select_run()
        else:
            self._submit_run()

    def _select_run(self):
        """Retrieve Run's object model from RDb"""
        try:
            run = db.session.query(Run).get(self.id)
        except IntegrityError as error:
            raise RunNotFoundError(f"Run {self.id} not found")
        if not run:
            raise RunNotFoundError(f"Run {self.id} not found")
        if run.user_id != self.user.id and not self.user.is_admin:
            raise RunAccessDenied(f"User {self.user.id} does not own Run {self.id}")
        self.model = run

    def _rpc_call(self, method, params={}):
        """This is not a Flask endpoint, but a helper used by several endpoints.
        It checks a user's permission to access a run, perform the specified RPC function,
        and returns result if OK, aborts if error.

        :param method: the method to execute remotely
        :type method: string
        :param params: parameters for the remote method, depends on method
        :type params: dict
        :return: response in JSON-RPC2 format
        :rtype: dict or list
        """
        a_site_rpc_client = rpc_index.rpc_index.get_client(self.model.site_id)
        params["user_id"] = self.user.id
        params["run_id"] = self.id
        params["cromwell_run_id"] = self.model.cromwell_run_id
        try:
            response = a_site_rpc_client.request(method, params)
        except Exception as error:
            logger.error(f"RPC {method} {params} failed: {error}")
            raise
        return response

    def _insert_run(self, params):
        """Insert new Run record into RDb"""
        site_id = params.get("site_id", None).upper()
        submission_id = params.get("submission_id")
        input_site_id = params.get("input_site_id", None).upper()
        input_endpoint = params.get("input_endpoint", None)
        output_endpoint = params.get("output_endpoint")
        output_dir = params.get("output_dir")
        wdl_file = params.get("wdl_file")
        json_file = params.get("json_file")
        tag = params.get("tag")

        run = Run(
            user_id=self.user.id,
            site_id=site_id,
            submission_id=submission_id,
            input_site_id=input_site_id,
            input_endpoint=input_endpoint,
            output_endpoint=output_endpoint,
            output_dir=output_dir,
            wdl_file=wdl_file,
            json_file=json_file,
            tag=tag,
            status="created",
        )
        try:
            db.session.add(run)
        except IntegrityError as error:
            db.session.rollback()
            logger.exception(f"Error inserting Run: {error}")
            raise
        try:
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            err_msg = f"Unable to insert new run in db: {error}"
            logger.exception(err_msg)
            raise
        self.model = run
        self.id = run.id

    def _submit_run(self, params):
        """
        Record the run submission in the database, with status as "uploading".

        :param params: New run parameters
        :type params: dict
        """
        site_id = params.get("site_id", None).upper()
        submission_id = params.get("submission_id")
        input_site_id = params.get("input_site_id", None).upper()
        input_endpoint = params.get("input_endpoint", None)
        output_endpoint = params.get("output_endpoint")
        output_dir = params.get("output_dir")
        wdl_file = params.get("wdl_file")
        json_file = params.get("json_file")
        tag = params.get("tag")

        compute_endpoint = config.conf.get_site(site_id, "globus_endpoint")

        if compute_endpoint is None:
            raise ValueError(
                f'Unknown Site ID, "{site_id}"; try the "list-sites" command'
            )
        logger.info(f"User {user}: New run submission {submission_id} to {site_id}")

        self._insert_run(params)
        logger.debug(f"User {user}: New run {run.id}")

        # Output directory is a subdirectory that includes the user id, site id and run id.
        # These are all placed in a common location with setgid sticky bits so that all
        # submitting users have access.
        output_dir += f"/{user}/{site_id}/{run.id}"
        src_host_path = config.conf.get_site(input_site_id, "globus_host_path")

        # We modify the output dir path since we know the endpoint of the returning source site. From here
        # a compute site can simply query the output directory and send.
        virtual_output_path = globus.virtual_transfer_path(output_dir, src_host_path)

        # Due to how the current database schema is setup, we have to update the output
        # directory from the model object itself immediately after insert.
        # TODO: Think of a better way to do this
        try:
            run.output_dir = virtual_output_path
        except Exception as error:
            db.session.rollback()
            err_msg = f"Unable to update output_dir in db: {error}"
            logger.exception(err_msg)
            abort(500, err_msg)
        try:
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            err_msg = f"Unable to update output_dir in db: {error}"
            logger.exception(err_msg)
            abort(500, err_msg)
        logger.debug(f"Updating output dir for run_id={self.id}")

        # SUBMIT GLOBUS TRANSFER
        manifest_file = request.files["manifest"]
        host_paths = {}
        host_paths["src"] = src_host_path
        host_paths["dest"] = config.conf.get_site(site_id, "globus_host_path")

        try:
            upload_task_id = globus.submit_transfer(
                f"Upload run {run.id}",
                host_paths,
                input_endpoint,
                compute_endpoint,
                manifest_file,
            )
        except globus_sdk.GlobusAPIError as error:
            run.status = "upload failed"
            db.session.commit()
            if error.code == "NoCredException":
                logger.warning(
                    f"{user} submission {run.id} failed due to Globus {error.code}"
                )
                abort(401, error.message)
            else:
                logger.exception(
                    f"{user} submission {run.id} failed for GlobusAPIError: {error}",
                    exc_info=True,
                )
                abort(error.code, error.message)
        except globus_sdk.NetworkError as error:
            logger.exception(
                f"{user} submission {run.id} failed due to NetworkError: {error}",
                exc_info=True,
            )
            abort(500, f"Network Error: {error}")
        except globus_sdk.GlobusError as error:
            logger.exception(
                f"{user} submission {run.id} failed for unknown error: {error}",
                exc_info=True,
            )
            abort(500, f"Unexpected error: {error}")

        logger.debug(f"User {user}: Run {run.id} upload {upload_task_id}")

        # UPDATE RUN WITH UPLOAD TASK ID AND ADD LOG ENTRY
        run.upload_task_id = upload_task_id
        old_status = run.status
        new_status = run.status = "uploading"
        log = Run_Log(
            run_id=self.id,
            status_to=new_status,
            status_from=old_status,
            timestamp=run.updated,
            reason=f"upload_task_id={upload_task_id}",
        )
        try:
            db.session.add(log)
        except Exception as error:
            db.session.rollback()
            logger.exception(f"Error inserting run log for Run {run.id}: {error}")
        try:
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            err_msg = f"Error inserting run log for Run {run.id}: {error}"
            logger.exception(err_msg)

        # SEND TO SITE
        params = {
            "run_id": self.id,
            "user_id": user,
            "email": current_user.email,
            "submission_id": submission_id,
            "upload_task_id": upload_task_id,
            "output_endpoint": output_endpoint,
            "output_dir": output_dir,
        }
        a_site_rpc_client = rpc_index.rpc_index.get_client(run.site_id)
        logger.debug(f"User {user}: submit run: {params}")
        try:
            result = a_site_rpc_client.request("submit", params)
        except Exception as error:
            reason = f"RPC submit failed: {error}"
            logger.exception(reason)
            _submission_failed(user, run, reason)
            abort(500, reason)
        if "error" in result:
            reason = f"Error sending new run to {site_id}: {result['error']['message']}"
            logger.error(reason)
            _submission_failed(user, run, reason)
            abort(result["error"]["code"], result["error"]["message"])

        # DONE
        result = {
            "run_id": self.id,
            "status": self.model.status,
            "site_id": self.model.site_id,
            "output_dir": self.model.output_dir,
            "tag": self.model.tag,
        }
        logger.info(f"User {self.user.id}: New run: {result}")
        return result, 201

    def _submission_failed(self, reason):
        """Cancel upload and update run status"""
        globus.cancel_transfer(self.model.upload_task_id)
        self._update_run_status("submission failed", reason)

    def _update_run_status(self, new_status, reason=None):
        """Update run table and insert run_logs entry."""
        status_from = self.model.status
        self.model.status = new_status
        try:
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.exception(f"Error updating run status in db: {error}")
        log = Run_Log(
            run_id=self.id,
            status_from=status_from,
            status_to=new_status,
            timestamp=self.model.updated,
            reason=reason,
        )
        try:
            db.session.add(log)
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.error(f"Error insert run log entry: {error}")

    @property
    def pre_cromwell(self):
        """Returns True if Run has not been submitted to Cromwell yet"""
        return True if self.status in run_pre_cromwell_states else False

    def run_status(self):
        """
        Retrieve the current status of a run.

        :return: The status of the run
        :rtype: dict
        """
        logger.info(f"Get status of Run {run.id}")
        result = {
            "id": self.model.id,
            "submission_id": self.model.submission_id,
            "cromwell_run_id": self.model.cromwell_run_id,
            "result": self.model.result,
            "status": self.model.status,
            "status_detail": jaws_constants.run_status_msg.get(self.model.status, ""),
            "site_id": self.model.site_id,
            "submitted": self.model.submitted.strftime("%Y-%m-%d %H:%M:%S"),
            "updated": self.model.updated.strftime("%Y-%m-%d %H:%M:%S"),
            "input_site_id": self.model.input_site_id,
            "input_endpoint": self.model.input_endpoint,
            "upload_task_id": self.model.upload_task_id,
            "output_endpoint": self.model.output_endpoint,
            "output_dir": self.model.output_dir,
            "download_task_id": self.model.download_task_id,
            "user_id": self.model.user_id,
            "tag": self.model.tag,
            "wdl_file": self.model.wdl_file,
            "json_file": self.model.json_file,
        }
        return result

    def task_status(self):
        """
        Retrieve the current status of each Task in a Run.

        :return: The status of each task in a run.
        :rtype: dict
        """
        logger.info(f"User {self.user.id}: Get task-status of Run {self.id}")
        if self.pre_cromwell:
            return None
        else:
            return self._rpc_call("get_task_status")

    def run_log(self):
        """
        Retrieve complete log of a Run's state transitions.

        :return: Table of log entries
        :rtype: list
        """
        logger.info(f"User {self.user.id}: Get log of Run {self.id}")
        log = RunLog(self.id)
        return log.get_log()

    def task_log(self):
        """
        Retrieve log of all task state transitions.

        :return: The complete log of all task state transitions.
        :rtype: dict
        """
        logger.info(f"User {self.user.id}: Get task-log for Run {self.id}")
        if self.pre_cromwell:
            return None
        else:
            return self._rpc_call("get_task_log")

    def run_metadata(self):
        """
        Retrieve the metadata of a run.

        :return: Cromwell metadata for the run, if found; abort otherwise
        :rtype: dict
        """
        logger.info(f"User {self.ser}: Get metadata for Run {self.id}")
        if self.pre_cromwell:
            return None
        else:
            return self._rpc_call("run_metadata")

    def get_errors(self):
        """
        Retrieve error messages and stderr for failed tasks.

        :return: Cromwell error messages and stderr for all failed tasks.
        :rtype: str
        """
        logger.info(f"User {self.user.id}: Get errors for Run {self.id}")
        if self.pre_cromwell:
            return None
        else:
            return self._rpc_call("get_errors")

    def cancel_run(self):
        """
        Cancel a run.  It doesn't cancel Globus transfers, just Cromwell runs.

        :return: OK message upon success; abort otherwise
        :rtype: dict
        """
        logger.info(f"User {self.user.id}: Cancel Run {self.id}")
        if self.model.status in ["cancelled", "download complete"]:
            return
        _cancel_run(run)
        if status.startswith("upload"):
            self._cancel_transfer(self.model.upload_task_id)
        elif status.startswith("download"):
            self._cancel_transfer(self.model.download_task_id)
        self._rpc_call("cancel_run")

    def _cancel_run(self, reason="Cancelled by user"):
        """Update database record."""
        status_from = self.model.status
        self.model.status = "cancelled"
        self.model.result = "cancelled"
        try:
            db.session.commit()
        except Exception as error:
            db.session.rollback()
            logger.exception(f"Error while updating run to 'cancelled': {error}")
        log = RunLog(self.id)
        log.add_log(status_from, self.model.status, self.model.updated, reason)


class SearchRuns:
    """Object to query runs db"""

    def __init__(self, session, user):
        """Initialize Runs query object.

        :param session: db handle
        :type session: sqlalchemy.Session
        :param user: user object
        :type user: jaws_central.User
        """
        self.session = session
        self.user = user

    def queue(self):
        """Return the current user's unfinished runs.

        :return: details about any current runs
        :rtype: list
        """
        logger.info(f"User {self.user.id}: Get queue")
        try:
            queue = (
                db.session.query(Run)
                .filter_by(user_id=self.user.id)
                .filter(Run.status.in_(run_active_states))
                .all()
            )
        except SQLAlchemyError as error:
            logger.error(error)
            raise

        result = []
        for run in queue:
            result.append(
                {
                    "id": run.id,
                    "submission_id": run.submission_id,
                    "cromwell_run_id": run.cromwell_run_id,
                    "result": run.result,
                    "status": run.status,
                    "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
                    "site_id": run.site_id,
                    "submitted": run.submitted.strftime("%Y-%m-%d %H:%M:%S"),
                    "updated": run.updated.strftime("%Y-%m-%d %H:%M:%S"),
                    "input_site_id": run.input_site_id,
                    "input_endpoint": run.input_endpoint,
                    "upload_task_id": run.upload_task_id,
                    "output_endpoint": run.output_endpoint,
                    "output_dir": run.output_dir,
                    "download_task_id": run.download_task_id,
                    "user_id": run.user_id,
                    "tag": run.tag,
                    "wdl_file": run.wdl_file,
                    "json_file": run.json_file,
                }
            )
        return result

    def history(self, delta_days=10):
        """Return the current user's recent runs, regardless of status.

        :param delta_days: number of days in which to search
        :type delta_days: int
        :return: details about any recent runs
        :rtype: list
        """
        start_date = datetime.today() - timedelta(int(delta_days))
        logger.info(f"User {self.user.id}: Get history, last {delta_days} days")
        try:
            history = (
                db.session.query(Run)
                .filter_by(user_id=self.user.id)
                .filter(Run.submitted >= start_date)
                .all()
            )
        except SQLAlchemyError as error:
            logger.exception(f"Failed to select run history: {error}")
            raise

        result = []
        for run in history:
            result.append(
                {
                    "id": run.id,
                    "submission_id": run.submission_id,
                    "cromwell_run_id": run.cromwell_run_id,
                    "result": run.result,
                    "status": run.status,
                    "status_detail": jaws_constants.run_status_msg.get(run.status, ""),
                    "site_id": run.site_id,
                    "submitted": run.submitted.strftime("%Y-%m-%d %H:%M:%S"),
                    "updated": run.updated.strftime("%Y-%m-%d %H:%M:%S"),
                    "input_site_id": run.input_site_id,
                    "input_endpoint": run.input_endpoint,
                    "upload_task_id": run.upload_task_id,
                    "output_endpoint": run.output_endpoint,
                    "output_dir": run.output_dir,
                    "download_task_id": run.download_task_id,
                    "user_id": run.user_id,
                    "tag": run.tag,
                    "wdl_file": run.wdl_file,
                    "json_file": run.json_file,
                }
            )
        return result
