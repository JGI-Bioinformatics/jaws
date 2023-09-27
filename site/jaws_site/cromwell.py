"""
Cromwell class provides object-oriented API to Cromwell REST server.

NOTE: This class employs the verbiage employed by Cromwell, which may differ slightly from that used
elsewhere in JAWS, as clarified below:

- Cromwell "workflow": the execution of a "Workflow Specification" on a particular set of input.
  - "Workflow Specification" is colloquially referred to as a "WDL file" or "widdle"
  - a Cromwell "Workflow Specification" is called a "Workflow" in JAWS
  - a Cromwell "Workflow" is referred to as a "Run" in JAWS
  - the JAWS "cromwell_run_id" is equivalent to the Cromwell "workflow_id"

- Cromwell "task": a specific step in a Workflow Specification
  - same in JAWS, which uses the variable "cromwell_task_name" for Cromwell's "task_name"
  - JAWS cromwell_task name has expanded format: "<WDL_name>.<task_name>"; if the task is a scatter then
    the format is "<WDL_name>.<task_name>[shard_index]"

- Cromwell "call": An execution of a task
  - a Task may be atempted more than once, each with a unique "attempt" (int) and "job_id"

- Cromwell "job": the execution of a Task on the backend
  - the "job_id" is assigned by the backend; e.g. local (pid), SLURM (slurm job_id), etc.
  - same in JAWS, which uses "backend_job_id" to refer to Cromwell's "job_id"
"""


import urllib3
import requests
import requests.adapters
import logging
import os
import glob
import json
import io
import re
import boto3
import botocore


REQUEST_TIMEOUT = 120
NUMBER_OF_RETRIES = 5
BACKOFF_FACTOR = 5


logger = logging.getLogger(__package__)

aws_access_key_id = os.environ.get("AWS_ACCESS_KEY_ID", None)
aws_secret_access_key = os.environ.get("AWS_SECRET_ACCESS_KEY", None)
aws_region_name = os.environ.get("AWS_REGION_NAME", None)

session = requests.Session()
retries = urllib3.Retry(
    total=NUMBER_OF_RETRIES,
    backoff_factor=BACKOFF_FACTOR,
    status_forcelist=[500, 502, 503, 504],
)

session.mount("http://", requests.adapters.HTTPAdapter(max_retries=retries))

time_re = re.compile(r"^\s*(\d+):(\d+):(\d+)\s*$")
memory_re = re.compile(r"^\s*(\d+)\s*(\w+)\s*$")


def s3_parse_uri(full_uri):
    """
    Extract the bucket name and object key from full URI
    :param full_uri: String containing bucket and obj key, starting with "s3://"
    :ptype full_uri: str
    :return: s3 bucket name, object key
    :rtype: list
    """
    full_uri = full_uri.replace("s3://", "", 1)
    folders = full_uri.split("/")
    s3_bucket = folders.pop(0)
    obj_key = "/".join(folders)
    return s3_bucket, obj_key


def _read_file_s3(path, **kwargs):
    """
    Read the contents of a file from S3 and return it's contents.
    """
    s3_bucket, src_path = s3_parse_uri(path)
    aws_session = boto3.Session(
        aws_access_key_id=kwargs["aws_access_key_id"],
        aws_secret_access_key=kwargs["aws_secret_access_key"],
        region_name=kwargs["aws_region_name"],
    )
    s3_resource = aws_session.resource("s3")
    bucket_obj = s3_resource.Bucket(s3_bucket)
    fh = io.BytesIO()
    try:
        bucket_obj.download_fileobj(src_path, fh)
    except botocore.exceptions.ClientError as error:
        raise IOError(f"File obj not found, {src_path}: {error}")
    except Exception as error:
        raise IOError(error)
    fh.seek(0)
    data = fh.read()
    fh.close()
    contents = data.decode("utf-8")
    return contents


def _read_file_nfs(path: str):
    """
    :param path: Path to file on NFS
    :ptype path: str
    :return: contents of file
    :rtype: str
    """
    if not os.path.isfile(path):
        raise IOError(f"File not found: {path}")
    contents = None
    try:
        with open(path, "r") as file:
            contents = file.read()
    except Exception as error:
        raise (f"Error reading file, {path}: {error}")
    return contents


def _read_file(path: str, **kwargs):
    """
    Read file from NFS or S3 and return contents.
    :param path: Path to file (may be s3 item)
    :ptype path: str
    :return: contents of file
    :rtype: str
    """
    if path.startswith("s3://"):
        return _read_file_s3(path, **kwargs)
    else:
        return _read_file_nfs(path)


def _write_file_s3(path: str, content: str, **kwargs):
    s3_bucket, src_path = s3_parse_uri(path)
    aws_session = boto3.Session(
        aws_access_key_id=kwargs["aws_access_key_id"],
        aws_secret_access_key=kwargs["aws_secret_access_key"],
        region_name=kwargs["aws_region_name"],
    )
    s3_resource = aws_session.resource("s3")
    bucket_obj = s3_resource.Bucket(s3_bucket)
    bucket_obj.put(Body=content)


def _write_file_nfs(path: str, content: str):
    with open(path, "w") as fh:
        fh.write(content)


def _write_file(path: str, contents: str, **kwargs):
    """
    Write contents to NFS or S3 file.
    :param path: Path to file (may be s3 item)
    :ptype path: str
    :param contents: Contents to write
    :ptype contents: str
    """
    if path.startswith("s3://"):
        return _write_file_s3(path, contents, **kwargs)
    else:
        return _write_file_nfs(path, contents)


def sort_table(table: list, index: int):
    """Sort table by specified column value.
    :param table: table to sort
    :ptype table: list
    :param index: column index to sort by
    :ptype index: int
    :return: sorted table
    :rtype: list
    """
    table.sort(key=lambda x: x[index])
    return table


def sort_table_dict(table: list, key: str):
    """Sort table by specified column value.
    :param table: table to sort
    :ptype table: list
    :param key: column key to sort by
    :ptype key: str
    :return: sorted table
    :rtype: list
    """
    table.sort(key=lambda x: x[key])
    return table


class CallError(Exception):
    pass


class Call:
    def __init__(self, data, task_name: str):
        """
        :param data: This is the original call record from the cromwell metadata
        :ptype data: dict
        """
        if "subWorkflowMetadata" in data:
            raise CallError(f"Task {task_name} is a subworkflow, not a Call object")
        self.data = data
        self.attempt = data.get("attempt", None)
        self.shard_index = data.get("shardIndex", -1)
        self.name = task_name
        self.execution_status = data.get("executionStatus", None)
        self.result = None
        self.cached = False
        self.return_code = data.get("returnCode", None)
        self.stdout = self.data.get("stdout", None)
        self.stderr = self.data.get("stderr", None)
        self.call_root = self.data.get("callRoot", None)
        self.execution_dir = None
        self.job_id = self.data.get("jobId", None)
        self.requested_time = None
        self.requested_time_minutes = None
        self.requested_memory = None
        self.requested_memory_gb = None
        self.requested_cpu = None
        self.failure_message = None
        self.dir = None

        if self.execution_status == "Failure":
            self.result = "failed"
        elif self.execution_status == "Done":
            self.result = "succeeded"

        if "runtimeAttributes" in self.data:
            self.requested_time = self.data["runtimeAttributes"].get("time", None)
            self.requested_memory = self.data["runtimeAttributes"].get("memory", None)
            if "cpu" in self.data["runtimeAttributes"]:
                self.requested_cpu = int(self.data["runtimeAttributes"]["cpu"])

        if "failures" in self.data:
            # save last failure message only
            for failure in self.data["failures"]:
                self.failure_message = failure["message"]

        if "callCaching" in self.data and "hit" in self.data["callCaching"]:
            self.cached = self.data["callCaching"]["hit"]

    def _get_file_path(self, file_id, relpath=False):
        """
        Return the path to the specified output file, optionally replacing part of the path.
        :param file_id: either "stdout" or "stderr"
        :type file_id: str
        :param relpath: If true then make path relative to callRoot; fullpath by default
        :type relpath: bool
        :param dest: Destination root dir
        :return: Path to specified file
        :rtype: str
        """
        if file_id not in ("stdout", "stderr"):
            raise TaskError(
                f"Invalid file id, {file_id}; allowed values: stdout, stderr"
            )
        path = self.data.get(file_id)
        if relpath:
            if self.call_root is None:
                raise TaskError(
                    f"Cannot return relpath as callRoot is not defined for task {self.name}"
                )
            path = os.path.relpath(self.call_root, path)
        return path

    def stdout(self, relpath=False):
        """
        Return the path to the standard output file, optionally replacing part of the path.
        :param relpath: If true then make path relative to callRoot; fullpath by default
        :type relpath: bool
        :param dest: Destination root dir
        :return: Path to stdout file
        :rtype: str
        """
        return self._get_file_path("stdout", relpath)

    def stderr(self, relpath=False):
        """
        Return the path to the standard err file, optionally replacing part of the path.
        :param relpath: If true then make path relative to callRoot; fullpath by default
        :type relpath: bool
        :return: Path to stdout file
        :rtype: str
        """
        return self._get_file_path("stderr", relpath)

    def summary(self):
        """
        :return: Select fields
        :rtype: dict
        """
        record = {
            "name": self.name,
            "shard_index": self.shard_index,
            "attempt": self.attempt,
            "cached": self.cached,
            "job_id": self.job_id,
            "execution_status": self.execution_status,
            "result": self.result,
            "failure_message": self.failure_message,
            "call_root": self.call_root,
            "requested_time": self.requested_time,
            "requested_cpu": self.requested_cpu,
            "requested_memory": self.requested_memory,
        }
        return record

    def failed_folder(self):
        """
        Return the 'callRoot' if the attempt failed and the folder is defined; otherwise return None.
        """
        if self.execution_status == "Failed" and "callRoot" in self.data:
            return self.data["callRoot"]
        else:
            return None

    def error(self):
        """
        Errors report, if failed (else None).
        :return: errors report
        :rtype: dict
        """
        if self.execution_status != "Failed":
            return None
        result = {
            "attempt": self.attempt,
            "shardIndex": self.shard_index,
        }
        if "failures" in self.data:
            result["failures"] = self.data["failures"]
        if "jobId" in self.data:
            result["jobId"] = self.data["jobId"]
        if "returnCode" in self.data:
            result["returnCode"] = self.data["returnCode"]
        if "runtimeAttributes" in self.data:
            result["runtimeAttributes"] = self.data["runtimeAttributes"]
        if "stderr" in self.data:
            # include *contents* of stderr files, instead of file paths
            stderr_file = self.data["stderr"]
            try:
                stderr_contents = _read_file(stderr_file)
            except Exception:  # noqa
                # stderr file doesn't always exist (e.g. fails in submit step)
                result["stderrContents"] = f"File not found: {stderr_file}"
            else:
                result["stderrContents"] = stderr_contents
            stderr_submit_file = f"{stderr_file}.submit"
            try:
                stderr_submit_contents = _read_file(stderr_submit_file)
            except Exception:  # noqa
                # stderrSubmit file doesn't always exist (not used on AWS)
                result["stderrSubmitContents"] = None
            else:
                result["stderrSubmitContents"] = stderr_submit_contents
            stderr_background_file = f"{stderr_file}.background"
            try:
                stderr_background_contents = _read_file(stderr_background_file)
            except Exception:  # noqa
                result["stderrBackgroundContents"] = None
            else:
                result["stderrBackgroundContents"] = stderr_background_contents
        if "stdout" in self.data:
            # include *contents* of stdout file, instead of file path
            stdout_file = self.data["stdout"]
            try:
                stdoutContents = _read_file(stdout_file)
            except Exception:  # noqa
                result["stdoutContents"] = f"File not found: {stdout_file}"
            else:
                result["stdoutContents"] = stdoutContents
            stdout_submit_file = f"{stdout_file}.submit"
            try:
                stdout_submit_contents = _read_file(stdout_submit_file)
            except Exception:  # noqa
                # stdoutSubmit file doesn't always exist (not used on AWS)
                result["stdoutSubmitContents"] = None
            else:
                result["stdoutSubmitContents"] = stdout_submit_contents
            stdout_background_file = f"{stdout_file}.background"
            try:
                stdout_background_contents = _read_file(stdout_background_file)
            except Exception:  # noqa
                result["stdoutBackgroundContents"] = None
            else:
                result["stdoutBackgroundContents"] = stdout_background_contents
        if "callRoot" in self.data:
            result["callRoot"] = self.data["callRoot"]
            # check for existence of any task-specific *.log files.
            # callRoot is not specified for AWS/S3 files, so it works only for NFS paths.
            log_files = glob.glob(f"{self.data['callRoot']}/execution/*.log")
            for log_file_path in log_files:
                log_file_name = os.path.basename(log_file_path)
                try:
                    log_file_contents = _read_file(log_file_path)
                except Exception:  # noqa
                    result[log_file_name] = None
                else:
                    result[log_file_name] = log_file_contents
        return result


class TaskError(Exception):
    pass


class Task:
    """
    A Task may have multiple calls, corresponding to multiple attempts and/or shards.
    If a Task is a subworkflow, it's Call object will contain a Metadata object.
    """

    def __init__(self, name, data, **kwargs):
        """
        Initialize a Task object, which may contain multiple shards or a subworkflow.

        :param name: Task name
        :type name: str
        :param data: Cromwell's "calls" for a task; this always remains unchanged.
        :type data: list
        """
        self.name = name
        self.data = data
        self.url = kwargs.get("url")
        self.calls = {}
        self.subworkflows = {}
        for call_data in data:
            # shardIndex is "-1" for regular (not scattered) tasks
            shard_index = int(call_data["shardIndex"])
            # in general, there is only 1 attempt
            attempt = int(call_data["attempt"])

            # if "expandSubWorkflows" option was used, "subWorkflowMetadata" will be populated,
            # otherwise "subWorkflowId" will be present; the latter requires another GET
            if "subWorkflowId" in call_data and self.url is not None:
                sub_id = call_data["subWorkflowId"]
                cromwell = Cromwell(self.url)
                metadata = cromwell.get_metadata(sub_id)
                call_data["subWorkflowMetadata"] = metadata.data
                if shard_index not in self.subworkflows:
                    self.subworkflows[shard_index] = {}
                self.subworkflows[shard_index][attempt] = metadata
            elif "subWorkflowMetadata" in call_data:
                sub_data = call_data["subWorkflowMetadata"]
                if shard_index not in self.subworkflows:
                    self.subworkflows[shard_index] = {}
                self.subworkflows[shard_index][attempt] = Metadata(sub_data)
            else:
                if shard_index not in self.calls:
                    self.calls[shard_index] = {}
                self.calls[shard_index][attempt] = Call(call_data, self.name)

    def failed_folders(self):
        """
        Return all 'callRoot' of this task or all child tasks if this is a subworkflow.
        """
        folders = []
        for shard_index in self.calls.keys():
            for attempt in self.calls[shard_index].keys():
                call = self.calls[shard_index][attempt]
                folder = call.failed_folder()
                if folder is not None:
                    folders.append(folder)
        for shard_index in self.subworkflows.keys():
            for attempt in self.subworkflows[shard_index].keys():
                sub_meta = self.subworkflows[shard_index][attempt]
                sub_failed_folders = sub_meta.failed_folders()
                if len(sub_failed_folders) > 0:
                    folders.extend(sub_failed_folders)
        return folders

    def errors(self):
        """
        Return user friendly errors report for this task/subworkflow.
        The contents of the stderr, stderr.submit files are also added for convenience.
        :return: Errors report
        :rtype: dict
        """
        all_errors = []
        for shard_index in self.calls.keys():
            for attempt in self.calls[shard_index].keys():
                call = self.calls[shard_index][attempt]
                call_error = call.error()
                if call_error is not None:
                    all_errors.append(call_error)
        for shard_index in self.subworkflows.keys():
            for attempt in self.subworkflows[shard_index].keys():
                sub_meta = self.subworkflows[shard_index][attempt]
                sub_meta_errors = sub_meta.errors()
                if len(sub_meta_errors.keys()) > 0:
                    sub_errors = {
                        "shardIndex": shard_index,
                        "attempt": attempt,
                        "subWorkflowMetadata": sub_meta_errors,
                    }
                    all_errors.append(sub_errors)
        return all_errors

    def summary(self, **kwargs):
        """
        Generally, there is only one attempt.  The only exception is when the resubmit with more memory
        feature of Cromwell is turned on, but it's only available with "gcs" backend as of 5/25/2022.
        For simplicity, we use only the last attempt by default.
        """
        last_attempts = kwargs.get("last_attempts", True)
        result = []
        for shard_index in self.calls.keys():
            attempts = sorted(self.calls[shard_index].keys())
            if last_attempts:
                attempts = [attempts[-1]]
            for attempt in attempts:
                call = self.calls[shard_index][attempt]
                result.append(call.summary())
        for shard_index in self.subworkflows.keys():
            attempts = sorted(self.subworkflows[shard_index].keys())
            if last_attempts:
                attempts = [attempts[-1]]
            for attempt in attempts:
                sub_meta = self.subworkflows[shard_index][attempt]
                result.extend(
                    sub_meta.task_summary(last_attempts=last_attempts, relpaths=False)
                )
        return result


class CromwellError(Exception):
    """Base class for all Cromwell exceptions."""

    pass


class CromwellServiceError(CromwellError):
    """Base class for all Cromwell service errors; such errors do not indicate a problem with the run"""

    pass


class CromwellServiceConnectionError(CromwellServiceError):
    """Cromwell is not responding"""

    pass


class CromwellRunError(CromwellError):
    """Base class for all errors that are specific to a run."""

    pass


class CromwellRunNotFoundError(CromwellRunError):
    """The specified run does not exist"""

    pass


class Metadata:
    """
    An object that represents the Cromwell Run metadata record.
    A Metadata object has many Tasks.
    """

    # Commonly used params:
    # "actualWorkflowLanguage" =>  str
    # "actualWorkflowLanguageVersion" =>  str
    # "calls" =>  dict
    # "end" =>  datetime
    # "id" =>  str
    # "inputs" =>  dict
    # "labels" =>  dict
    # "outputs" =>  dict
    # "parentWorkflowId" =>  str (iff subworkflow)
    # "rootWorkflowId" =>  str (iff subworkflow)
    # "start" =>  datetime
    # "status" =>  str
    # "submission" =>  datetime
    # "submittedFiles" =>  dict
    # "workflowName" =>  str
    # "workflowProcessingEvents" =>  list
    # "workflowRoot" =>  str (path)

    def __init__(self, data, url=None):
        """
        Initialize Metadata object by retrieving metadata from Cromwell,
        including subworkflows, and initialize Task objects.
        :param workflow_id: Cromwell's UUID for the workflow (AKA run)
        :type workflow_id: str
        :param data: cromwell metadata JSON record
        :type dict:
        """
        logger = logging.getLogger(__package__)
        self.data = data
        self.url = url
        self.workflow_id = data["id"]
        self.tasks = {}  # Task objects
        if "calls" in self.data:
            calls = self.data["calls"]
            for task_name, task_data in calls.items():
                logger.debug(f"Workflow {self.workflow_id}: Init task {task_name}")
                self.tasks[task_name] = Task(task_name, task_data, url=self.url)

    def get(self, param, default=None):
        """
        Get a section of the document.
        :param param: key of parameter to get
        :type param: str
        :param default: Default value if param undefined
        """
        return self.data.get(param, default)

    def workflow_name(self):
        return self.get("workflowName", None)

    def workflow_root(self, **kwargs):
        """
        Return the base location for the Run.
        Note that workflowRoot is only defined when using NFS;
        for S3 we must infer from an output file.
        """
        workflow_root = self.get("workflowRoot", None)
        if workflow_root:
            return workflow_root

        # AWS runs don't have workflowRoot defined, so construct path
        if "executions_dir" not in kwargs:
            raise ValueError(
                "'workflowRoot' is not defined in metadata; 'executions_dir' parameter is required."
            )
        executions_dir = kwargs["executions_dir"]
        workflow_name = self.get("workflowName", None)
        if workflow_name is None:
            raise ValueError(
                "workflowName' is not defined in metadata; cannot generate workflow_root"
            )
        workflow_id = self.get("id", None)
        if workflow_id is None:
            raise ValueError(
                "id' is not defined in metadata; cannot generate workflow_root"
            )
        workflow_root = f"{executions_dir}/{workflow_name}/{workflow_id}"
        return workflow_root

    def filtered_failures(self):
        """
        Filter failures with message "Workflow failed", as those failures are duplicated in the Tasks.
        :return: list of failures
        :rtype: list
        """
        failures = self.get("failures", [])
        filtered_failures = []
        for failure in failures:
            if failure["message"].lower() != "workflow failed":
                filtered_failures.append(failure)
        return filtered_failures

    def failed_folders(self, **kwargs):
        """
        Return list of folders of failed tasks (when 'callRoot' exists).
        """
        all_folders = []
        for task_name, task in self.tasks.items():
            task_folders = task.failed_folders()
            if len(task_folders) > 0:
                all_folders.extend(task_folders)
        if len(all_folders):
            root = self.workflow_root(**kwargs)
            if not root:
                raise ValueError("The workflowRoot could not be determined")
            for i in range(len(all_folders)):
                if all_folders[i].startswith(root):
                    all_folders[i] = os.path.relpath(all_folders[i], root)
        return all_folders

    def errors(self):
        """
        Return JSON errors report.
        """
        filtered_metadata = {}
        calls = {}
        for task_name, task in self.tasks.items():
            task_errors = task.errors()
            if len(task_errors):
                calls[task_name] = task_errors
        if len(calls):
            filtered_metadata["calls"] = calls
        other_failures = self.filtered_failures()
        if len(other_failures):
            filtered_metadata["failures"] = other_failures
        return filtered_metadata

    def task_summary(self, **kwargs):
        """
        Return select info about each task, including any subworkflows.
        :return: List of task info dictionaries.
        :rtype: list
        """
        last_attempts = kwargs.get("last_attempts", True)
        relpaths = kwargs.get("relpaths", True)
        summary = []
        for task_name, task in self.tasks.items():
            summary.extend(task.summary(last_attempts=last_attempts))
        if relpaths is True:
            root = self.workflow_root() + "/"
            for i in range(0, len(summary)):
                if "call_root" in summary[i] and summary[i]["call_root"] is not None:
                    summary[i]["call_root"] = summary[i]["call_root"].removeprefix(root)
        return summary

    def task_summary_dict(self, **kwargs):
        """
        Return select info about each task, including any subworkflows.
        :return: Dictionary of task information dictionaries where relpath of call_root is key.
        :rtype: dict
        """
        workflow_root = self.workflow_root()
        summary = self.task_summary(**kwargs)
        result = {}
        for task in summary:
            call_root = task["call_root"]
            del task["call_root"]
            relpath = call_root.removeprefix(workflow_root)
            result[relpath] = task
        return result

    def job_summary(self, **kwargs) -> dict:
        """
        :return: task name for each job id
        :rtype: dict
        """
        job_ids = {}
        for task in self.task_summary():
            if "job_id" in task:
                job_id = str(task["job_id"])
                job_ids[job_id] = task.get("name", "unknown")
        return job_ids

    def outputs(self, **kwargs):
        """Returns all outputs for a workflow"""
        relpaths = kwargs.get("relpaths", False)
        outputs = self.get("outputs", {})
        if relpaths is False:
            return outputs

        # make paths relative to "{workflow_root}/{workflow_name}"
        # i.e. the path starting with the cromwell run id
        workflow_root = self.workflow_root(**kwargs)

        if not workflow_root:
            raise ValueError("The workflowRoot could not be determined.")
        elif len(outputs.keys()) == 0:
            return {}

        # convert to string, convert all paths to relpaths via replace, then convert back to dict
        output_str = json.dumps(outputs, indent=0)
        output_str = output_str.replace(workflow_root, ".")
        relpath_outputs = json.loads(output_str)
        return relpath_outputs

    def outfiles(self, **kwargs):
        """
        Return list of all output files of a run.
        Requires that the workflow_root can be determined.
        :return: Relative paths of all output files
        :rtype: list
        """
        workflow_root = self.workflow_root(**kwargs)
        outputs = self.get("outputs", {})

        if not workflow_root:
            raise ValueError("The workflowRoot could not be determined.")
        elif len(outputs.keys()) == 0:
            return []

        relpath_outputs = []

        # find all the file paths
        output_str = json.dumps(outputs, indent=0)
        regex = f"{workflow_root}\\S+"
        matches = re.findall(regex, output_str)

        # change full paths to relative paths
        for myfile in matches:
            relpath = myfile.replace(workflow_root, ".")
            relpath = re.sub(r"[\"\}\],]", "", relpath)
            relpath_outputs.append(relpath)

        return relpath_outputs


class Cromwell:
    """Class representing a Cromwell REST server."""

    workflows = {}

    def __init__(self, url: str, caller_logger=None) -> None:
        """Init object with connection info.

        :param url: url of Cromwell REST server, including port
        :type param: str
        """
        if not url.startswith("http"):
            url = f"http://{url}"
        self.url = url
        self.workflows_url = f"{url}/api/workflows/v1"
        self.engine_url = f"{url}/engine/v1/status"
        if caller_logger is not None:
            global logger
            logger = caller_logger

    def get(self, workflow_id: str, output_type: str):
        """
        Get the specified Cromwell output for a Run.
        :param workflow_id: unique ID of the run
        :ptype workflow_id: str
        :param output_type: Report name (REST endpoint)
        :ptype output_type: str
        :return: JSON report from Cromwell
        :rtype: dict
        """
        url = f"{self.workflows_url}/{workflow_id}/{output_type}"
        try:
            response = session.get(url, timeout=REQUEST_TIMEOUT)
        except requests.exceptions.ConnectionError as error:
            raise CromwellServiceError(f"Unable to reach Cromwell service: {error}")
        except requests.exceptions.RetryError:
            logger.exception(
                "Reached max retries to Cromwell server. Check cromwell server status"
            )
            raise CromwellServiceError("Max retries for Cromwell server reached")
        if response.status_code == 404:
            raise CromwellRunNotFoundError(f"Cromwell run {workflow_id} not found")
        elif response.status_code >= 400:
            raise CromwellServiceError(
                f"Error retrieving Cromwell Run report {output_type}: code {response.status_code}"
            )
        return response.json()

    def get_metadata(self, workflow_id: str, data=None):
        """Get Metadata object for a workflow-run.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :param data: optionally provide the metadata instead
        :type data: dict
        :return: Metadata object
        :rtype: cromwell.Metadata
        """
        # don't expand subworkflows because it fails for very large runs (>1M rows)
        data = self.get(workflow_id, "metadata?expandSubWorkflows=0")
        return Metadata(data, url=self.url)

    def status(self):
        """Check if Cromwell is available"""
        try:
            response = session.get(self.engine_url, timeout=REQUEST_TIMEOUT)
        except requests.exceptions.ConnectionError as error:
            raise CromwellServiceError(f"Unable to reach Cromwell service: {error}")
        if response.status_code >= 400:
            raise CromwellServiceError(
                f"Error retrieving Cromwell status: code {response.status_code}"
            )
        return True

    def abort(self, workflow_id: str):
        """Abort a run.  Raise upon error."""
        url = f"{self.workflows_url}/{workflow_id}/abort"
        try:
            response = session.post(url, timeout=REQUEST_TIMEOUT)
        except requests.exceptions.ConnectionError as error:
            raise CromwellServiceError(
                f"Workflow {workflow_id}: Abort; Unable to reach Cromwell service: {error}"
            )
        sc = response.status_code
        if sc == 200 or sc == 403:
            # 403 = finishing already; too late to cancel (do not raise)
            return response.json()
        elif sc == 400:
            raise CromwellRunError(
                f"Workflow {workflow_id}: Abort failed for malformed workflow id"
            )
        elif sc == 404:
            raise CromwellRunNotFoundError(
                f"Workflow {workflow_id}: Abort; id not found"
            )
        else:
            raise CromwellError(
                f"Workflow {workflow_id}: Abort; An unexpected Cromwell error {sc} occurred"
            )

    def _options_fh(self, **kwargs):
        """
        Provide a file handle to an in-RAM JSON file with the Cromwell options.
        """
        options = {"default_runtime_attributes": {"docker": "ubuntu:latest"}}
        if "default_container" in kwargs:
            options["default_runtime_attributes"]["docker"] = kwargs[
                "default_container"
            ]
        if "caching" in kwargs and kwargs["caching"] is False:
            options["read_from_cache"] = False
            options["write_to_cache"] = False
        else:
            options["read_from_cache"] = True
            options["write_to_cache"] = True
        fh = io.StringIO(json.dumps(options))
        fh.seek(0)
        return fh

    def submit(self, file_handles: dict, options: dict) -> int:
        """
        Submit a run to Cromwell.
        :return: Cromwell workflow uuid
        :rtype: str
        """
        files = {}
        if "wdl" in file_handles:
            files["workflowSource"] = (
                "workflowSource",
                file_handles["wdl"],
                "application/json",
            )
        else:
            raise CromwellRunError("WDL file handle required")
        if "inputs" in file_handles:
            files["workflowInputs"] = (
                "workflowInputs",
                file_handles["inputs"],
                "application/json",
            )
        else:
            raise CromwellRunError("Inputs JSON file handle required")
        if "subworkflows" in file_handles:
            files["workflowDependencies"] = (
                "workflowDependencies",
                file_handles["subworkflows"],
                "application/zip",
            )
        files["workflowOptions"] = (
            "workflowOptions",
            self._options_fh(**options),
            "application/json",
        )
        try:
            response = session.post(
                self.workflows_url, files=files, timeout=REQUEST_TIMEOUT
            )
        except requests.exceptions.ConnectionError as error:
            raise CromwellServiceError(f"Unable to reach Cromwell service: {error}")
        if response.status_code >= 400:
            raise CromwellServiceError(
                f"Error retrieving Cromwell status: code {response.status_code}"
            )
        run_id = response.json()["id"]
        return run_id

    def get_status(self, workflow_id: str):
        """
        Get the status of a workflow.
        :param workflow_id: Cromwell's workflow uuid
        :type workflow_id: str
        :return: Status of workflow
        :rtype: str
        """
        result = self.get(workflow_id, "status")
        return result["status"]
