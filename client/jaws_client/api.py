"""
JAWS Client API
"""

import os
import requests
import json
import uuid
import shutil
from jaws_client import log as logging
from jaws_client.config import Configuration

JAWS_LOG_ENV = "JAWS_CLIENT_LOG"
JAWS_USER_LOG = os.path.expanduser("~/jaws.log")
JAWS_CONFIG_ENV = "JAWS_CLIENT_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-client.conf")
JAWS_USER_CONFIG_ENV = "JAWS_USER_CONFIG"
JAWS_USER_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws.conf")


logger = None
config = None


class JawsClientError(Exception):
    pass


class JawsUserError(JawsClientError):
    pass


class JawsServiceError(JawsClientError):
    pass


class Client:
    def __init__(
        self,
        jaws_config_file: str = None,
        user_config_file: str = None,
        log_file: str = None,
        log_level: str = "INFO",
    ):
        if log_file is None:
            log_file = (
                os.environ[JAWS_LOG_ENV]
                if JAWS_LOG_ENV in os.environ
                else JAWS_USER_LOG
            )
        global logger
        logger = logging.setup_logger(__package__, log_file, log_level)

        if jaws_config_file is None:
            jaws_config_file = (
                os.environ[JAWS_CONFIG_ENV]
                if JAWS_CONFIG_ENV in os.environ
                else JAWS_CONFIG_DEFAULT_FILE
            )
        os.environ[JAWS_CONFIG_ENV] = jaws_config_file
        if user_config_file is None:
            user_config_file = (
                os.environ[JAWS_USER_CONFIG_ENV]
                if JAWS_USER_CONFIG_ENV in os.environ
                else JAWS_USER_CONFIG_DEFAULT_FILE
            )
        os.environ[JAWS_USER_CONFIG_ENV] = user_config_file
        global config
        config = Configuration(jaws_config_file, user_config_file)

    def _request(self, rest_op, urn, data={}, files={}) -> dict:
        """Perform specified REST operation.  A JSON response is expected."""
        if config is None:
            raise JawsClientError(
                "The config obj must be initialized before using this function"
            )
        access_token = config.get("USER", "token")
        if not access_token:
            raise JawsServiceError(
                "User access token required; contact an admin to get yours."
            )
        header = {"Authorization": f"Bearer {access_token}"}
        url = f'{config.get("JAWS", "url")}/{urn}'
        response = None
        try:
            response = requests.request(rest_op, url, headers=header, data=data, files=files)
        except requests.exceptions.Timeout as err:
            raise JawsServiceError("Unable to communicate with JAWS server (timeout)", err)
        except requests.exceptions.TooManyRedirects as err:
            raise JawsServiceError(
                "Unable to communicate with JAWS server (too many redirects; bad url?)",
                err,
            )
        except requests.exceptions.HTTPError as err:
            raise JawsServiceError(
                "Unable to communicate with JAWS server (http error)", err
            )
        except requests.exceptions.RequestException as err:
            raise JawsServiceError("Unable to communicate with JAWS server", err)
        if response.status_code < 200 or response.status_code > 299:
            try:
                result = response.json()
            except Exception:
                raise JawsServiceError(response.text)
            if "error" in result:
                raise JawsServiceError(result["error"])
            else:
                raise JawsServiceError(result)
        return response.json()

    def health(self) -> dict:
        """Current system status."""
        url = f'{config.get("JAWS", "url")}/status'
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            raise JawsServiceError("JAWS Central is DOWN")
        if r.status_code != 200:
            raise JawsServiceError(r.text)
        return r.json()

    def info(self) -> dict:
        """JAWS version and info."""
        url = f'{config.get("JAWS", "url")}/info'
        try:
            r = requests.get(url)
        except requests.exceptions.RequestException:
            raise JawsServiceError("JAWS Central is DOWN")
        if r.status_code != 200:
            raise JawsServiceError(r.text)
        return r.json()

    def queue(self, site: str = "ALL") -> dict:
        """List of user's current runs"""
        data = {
            "delta_days": 0,
            "site_id": site.upper(),
            "active_only": True,
            "result": "any",
        }
        return self._request("POST", "search", data)

    def history(self, days: int = 1, site: str = "ALL", result: str = "any") -> dict:
        """Print a list of the user's past runs."""
        if days < 1:
            raise JawsServiceError("User error: --days must be a positive integer")
        if site:
            site = site.upper()
        data = {
            "delta_days": days,
            "site_id": site.upper(),
            "active_only": False,
            "result": result.lower(),
        }
        return self._request("POST", "search", data)

    def status(self, run_id: int, verbose: bool = False) -> dict:
        """Print the current status of a run."""

        if verbose is True:
            urn = f'run/{run_id}/complete'
        else:
            urn = f'run/{run_id}'
        return self._request("GET", urn)

    def task_status(self, run_id: int, fmt: str = "text") -> dict:
        """Show the current status of each Task."""

        urn = f'run/{run_id}/task_status'
        result = self._request("GET", urn)
        header = [
            "#CROMWELL_RUN_ID",
            "TASK_NAME",
            "ATTEMPT",
            "CROMWELL_JOB_ID",
            "STATUS_FROM",
            "STATUS_TO",
            "TIMESTAMP",
            "REASON",
        ]
        if fmt == "json":
            return result
        else:
            result_txt = "\t".join(header)
        for row in result:
            row[2] = str(row[2])
            row[3] = str(row[3])
            result_txt += "\t".join(row)
            return result_txt

    def metadata(self, run_id: int) -> dict:
        """Print detailed metadata for a run, generated by cromwell."""

        urn = f'run/{run_id}/metadata'
        return self._request("GET", urn)

    def log(self, run_id: int, fmt: str = "text") -> dict:
        """View the log of Run state transitions for the workflow as a whole."""

        urn = f'run/{run_id}/run_log'
        result = self._request("GET", urn)
        header = ["#STATUS_FROM", "STATUS_TO", "TIMESTAMP", "REASON"]
        if fmt == "json":
            return result
        elif fmt == "tab":
            result_txt = "\t".join(header)
            for log_entry in result:
                result_txt += "\t".join(log_entry)
            return result_txt
        else:
            """Get the max length of element in every col and add padding (2)"""
            result.insert(0, header)
            col_widths = []
            for idx in range(len(header)):
                col_widths.append(max(len(log_entry[idx]) for log_entry in result) + 2)
            result_txt = ""
            for log_entry in result:
                result_txt += "".join(
                    cell.ljust(col_widths[col_idx])
                    for col_idx, cell in enumerate(log_entry)
                )
            return result_txt

    def task_log(self, run_id: int, fmt: str) -> dict:
        """Get log of each Task's state transitions."""

        urn = f'run/{run_id}/task_log'
        result = self._request("GET", urn)
        header = [
            "#CROMWELL_RUN_ID",
            "TASK_NAME",
            "ATTEMPT",
            "CROMWELL_JOB_ID",
            "STATUS_FROM",
            "STATUS_TO",
            "TIMESTAMP",
            "REASON",
        ]
        if fmt == "json":
            return result
        elif fmt == "tab":
            result_txt = "\t".join(header)
            for row in result:
                row[2] = str(row[2])
                row[3] = str(row[3])
                result_txt += "\t".join(row)
            return result_txt
        else:
            """Get the max length of element in every col and add padding (2)"""
            result.insert(0, header)
            col_widths = []
            for idx in range(len(header)):
                col_widths.append(max(len(log_entry[idx]) for log_entry in result) + 2)
            result_txt = ""
            for log_entry in result:
                result_txt += "".join(
                    cell.ljust(col_widths[col_idx])
                    for col_idx, cell in enumerate(log_entry)
                )
            return result_txt

    def errors(self, run_id: int) -> dict:
        """View error messages and stderr for failed Tasks."""

        urn = f'run/{run_id}/errors'
        return self._request("GET", urn)

    def cancel(self, run_id: int) -> dict:
        """Cancel a run; returns whether aborting was successful or not."""

        urn = f'run/{run_id}/cancel'
        return self._request("PUT", urn)

    def cancel_all(self) -> dict:
        """Cancel all active runs."""
        urn = f'run/cancel-all'
        return self._request("PUT", urn)

    def list_sites(self) -> dict:
        """List available compute Sites"""
        return self._request("GET", "site")

    def submit(
        self,
        wdl_file: str,
        json_file: str,
        site: str,
        tag: str = None,
        no_cache: bool = False,
    ) -> dict:
        """Submit a run for execution at a JAWS-Site.
        Available sites can be found by running 'jaws run list-sites'.
        """
        from jaws_client import workflow
        from jaws_client.workflow import WdlError

        wdl_file = os.path.abspath(wdl_file)
        json_file = os.path.abspath(json_file)

        # the users' jaws id may not match the linux uid where the client is installed
        result = self._request("GET", "user")
        uid = result["uid"]

        staging_subdir = config.get("JAWS", "staging_dir")
        staging_user_subdir = os.path.join(staging_subdir, uid)
        globus_host_path = config.get("GLOBUS", "host_path")
        output_directory = config.get("JAWS", "data_repo_basedir")
        input_site_id = config.get("JAWS", "site_id")
        local_staging_endpoint = workflow.join_path(
            globus_host_path, staging_user_subdir
        )

        # GET SITE INFO
        compute_site_id = site.upper()
        urn = f'site/{compute_site_id}'
        result = self._request("GET", urn)
        compute_basedir = result["globus_host_path"]
        compute_uploads_subdir = result["uploads_dir"]
        compute_max_ram_gb = int(result["max_ram_gb"])

        # VALIDATE WORKFLOW WDLs
        submission_id = str(uuid.uuid4())
        try:
            wdl = workflow.WdlFile(wdl_file, submission_id)
            wdl.validate()
            max_ram_gb = wdl.max_ram_gb
        except WdlError as error:
            raise JawsServiceError(error)
        if max_ram_gb > compute_max_ram_gb:
            raise JawsServiceError(
                f"The workflow requires {max_ram_gb}GB but {compute_site_id} has only {compute_max_ram_gb}GB available"
            )

        # any and all subworkflow WDL files must be supplied to Cromwell in a single ZIP archive
        try:
            staged_wdl, zip_file = wdl.compress_wdls(local_staging_endpoint)
        except IOError as error:
            raise JawsServiceError(f"Unable to copy WDLs to inputs dir: {error}")

        # VALIDATE INPUTS JSON
        try:
            inputs_json = workflow.WorkflowInputs(json_file, submission_id)
        except json.JSONDecodeError as error:
            raise JawsServiceError(
                f"Your file, {json_file}, is not a valid JSON file: {error}"
            )

        staged_json = workflow.join_path(
            local_staging_endpoint, f"{submission_id}.json"
        )
        site_subdir = workflow.join_path(local_staging_endpoint, input_site_id)
        jaws_site_staging_dir = workflow.join_path(
            compute_basedir, compute_uploads_subdir
        )
        jaws_site_staging_site_subdir = workflow.join_path(
            jaws_site_staging_dir, input_site_id
        )

        # copy infiles in inputs json to site's inputs dir so they may be read by jaws user and
        # transferred to the compute site via globus
        moved_files = inputs_json.move_input_files(site_subdir)

        # the paths in the inputs json file are changed to their paths at the compute site
        modified_json = inputs_json.prepend_paths_to_json(jaws_site_staging_site_subdir)
        modified_json.write_to(staged_json)

        # the original inputs json file is kept with the run submission for record-keeping only
        orig_json = workflow.join_path(
            local_staging_endpoint, f"{submission_id}.orig.json"
        )
        try:
            shutil.copy(json_file, orig_json)
        except IOError as error:
            raise JawsServiceError(f"Error copying JSON to {orig_json}: {error}")
        try:
            os.chmod(orig_json, 0o0664)
        except PermissionError as error:
            raise JawsServiceError(f"Unable to chmod {orig_json}: {error}")

        # turning off call-caching requires a Cromwell options json file be created
        options_json_file = None
        if no_cache is True:
            options_json_file = workflow.join_path(
                local_staging_endpoint, f"{submission_id}.options.json"
            )
            with open(options_json_file, "w") as fh:
                fh.write('{"read_from_cache": false, "write_to_cache": false}')

        # write the file transfer manifest; jaws-central shall submit the transfer to globus
        manifest_file = workflow.Manifest(
            local_staging_endpoint, compute_uploads_subdir
        )
        manifest_file.add(staged_wdl, zip_file, staged_json, orig_json, *moved_files)
        if options_json_file:
            manifest_file.add(options_json_file)
        staged_manifest = workflow.join_path(
            staging_user_subdir, f"{submission_id}.tsv"
        )
        manifest_file.write_to(staged_manifest)

        # SUBMIT RUN TO CENTRAL
        local_endpoint_id = config.get("GLOBUS", "endpoint_id")
        data = {
            "site_id": compute_site_id,
            "submission_id": submission_id,
            "input_site_id": input_site_id,
            "input_endpoint": local_endpoint_id,
            "output_endpoint": local_endpoint_id,  # return to original submission site
            "output_dir": output_directory,  # jaws-writable dir to initially receive results
            "wdl_file": wdl_file,
            "json_file": json_file,
            "tag": tag,
        }
        files = {"manifest": open(staged_manifest, "r")}
        logger.debug(f"Submitting run: {data}")
        result = self._request("POST", "run", data, files)
        result["max_ram_gb"] = max_ram_gb
        del result["output_dir"]
        return result

    def inputs(self, wdl_file: str) -> dict:
        """Generate inputs template (JSON) from workflow (WDL) file."""

        from jaws_client import workflow

        if not os.path.isfile(wdl_file):
            raise JawsServiceError(f"File not found: {wdl_file}")
        stdout, stderr = workflow.womtool("inputs", wdl_file)
        if stderr:
            raise JawsServiceError(stderr)
        return stdout.strip()

    def validate(self, wdl_file: str) -> None:
        """Validate a WDL using Cromwell's WOMTool."""

        from jaws_client import workflow

        if not os.path.isfile(wdl_file):
            raise JawsServiceError(f"File not found: {wdl_file}")
        stdout, stderr = workflow.womtool("inputs", wdl_file)
        if stderr:
            raise JawsUserError(stderr)

    def get(self, run_id: int, dest: str) -> None:
        """Copy the output of a run to the specified folder."""

        from jaws_client import workflow

        result = self.status(run_id, True)
        status = result["status"]
        src = result["output_dir"]

        if status != "download complete":
            raise JawsServiceError(
                f"Run {run_id} output is not yet available; status is {status}"
            )

        if src is None:
            logger.error(f"Run {run_id} doesn't have an output_dir defined")
            raise JawsServiceError(f"Run {run_id} doesn't have an output_dir defined")

        try:
            result = workflow.rsync(
                src,
                dest,
                [
                    "-rLtq",
                    "--chmod=Du=rwx,Dg=rwx,Do=,Fu=rw,Fg=rw,Fo=",
                    "--exclude=inputs",
                    "--exclude=tmp.*",
                ],
            )
        except IOError as error:
            logger.error(f"Rsync output failed for run {run_id}: {error}")
            raise JawsServiceError(f"Error getting output for run {run_id}: {error}")
        if result.returncode != 0:
            err_msg = f"Failed to rsync {src}->{dest}: {result.stdout}; {result.stderr}"
            raise JawsServiceError(err_msg)

    def add_user(self, uid: str, email: str, name: str, admin: bool = False) -> dict:
        """Add new user and get JAWS OAuth access token (restricted)."""
        data = {"uid": uid, "email": email, "name": name, "admin": admin}
        return self._request("POST", "user", data)
