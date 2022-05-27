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


import requests
import logging
import os
import json
import io
from datetime import datetime, timezone
from dateutil import parser
from collections import deque
import boto3
import botocore


logger = logging.getLogger(__package__)

aws_access_key_id = os.environ.get("AWS_ACCESS_KEY_ID", None)
aws_secret_access_key = os.environ.get("AWS_SECRET_ACCESS_KEY", None)
aws_region_name = os.environ.get("AWS_REGION_NAME", None)


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
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name
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
    except IOError as error:
        raise IOError(f"Error reading file, {path}: {error}")
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
        self.name = task_name
        self.attempt = data.get("attempt", None)
        self.shard_index = data.get("shardIndex", None)
        self.execution_status = data.get("executionStatus", None)
        self.result = None
        self.cached = False
        self.return_code = data.get("returnCode", None)
        self.start = parser.parse(data["start"]).strftime("%Y-%m-%d %H:%M:%S")
        self.end = None
        if "end" in data:
            self.end = parser.parse(data["end"]).strftime("%Y-%m-%d %H:%M:%S")
        self.stdout = self.data.get("stdout", None)
        self.stderr = self.data.get("stderr", None)
        self.call_root = self.data.get("callRoot", None)
        self.execution_dir = None
        self.job_id = self.data.get("jobId", None)
        self.requested_time = None
        self.requested_memory = None
        self.requested_cpu = None
        self.failure_message = None

        if self.execution_status == "Failure":
            self.result == "failed"
        elif self.execution_status == "Done":
            self.result == "succeeded"

        if "runtimeAttributes" in self.data:
            self.requested_time = self.data["runtimeAttributes"].get("time", None)
            self.requested_memory = self.data["runtimeAttributes"].get("memory", None)
            self.requested_cpu = self.data["runtimeAttributes"].get("cpu", None)

        if "failures" in self.data:
            # save last failure message only
            for failure in self.data["failures"]:
                self.failure_message = failure["message"]

        self.queue_start = None
        self.run_start = None
        self.run_end = None
        self.queue_duration = None
        self.run_duration = None
        self.dir = None

        if "callCaching" in self.data and "hit" in self.data["callCaching"]:
            self.cached = self.data["callCaching"]["hit"]

        # get queued, running, and completed times as well as durations
        if "executionEvents" in self.data:
            for event in self.data["executionEvents"]:
                if event["description"] == "RequestingExecutionToken":
                    self.queue_start = parser.parse(event["startTime"]).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
                elif event["description"] == "RunningJob":
                    self.run_start = parser.parse(event["startTime"]).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
                elif event["description"] == "CallCacheReading":
                    self.run_start = parser.parse(event["startTime"]).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
                elif event["description"] == "UpdatingJobStore":
                    self.run_end = parser.parse(event["startTime"]).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
        if self.queue_start is None:
            self.queue_start = self.start
        if self.queue_start is not None and self.run_start is not None:
            delta = parser.parse(self.run_start) - parser.parse(self.queue_start)
            self.queue_duration = str(delta)
        if self.run_start is not None and self.run_end is not None:
            delta = parser.parse(self.run_end) - parser.parse(self.run_start)
            self.run_duration = str(delta)

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

    def summary(self, **kwargs):
        """
        :return: Select fields
        :rtype: dict
        """
        real_time = True
        if "real_time" in kwargs and kwargs["real_time"] is False:
            real_time = False
        if real_time is True and self.execution_status == "Running":
            self.set_real_time_status()
        record = {
            "name": self.name,
            "shard_index": self.shard_index,
            "attempt": self.attempt,
            "cached": self.cached,
            "job_id": self.job_id,
            "execution_status": self.execution_status,
            "result": self.result,
            "failure_message": self.failure_message,
            "queue_start": self.queue_start,
            "run_start": self.run_start,
            "run_end": self.run_end,
            "queue_duration": self.queue_duration,
            "run_duration": self.run_duration,
            "call_root": self.call_root,
            "requested_time": self.requested_time,
            "requested_cpu": self.requested_cpu,
            "requested_memory": self.requested_memory,
        }
        return record

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
            result["stderrContents"] = _read_file(stderr_file)
            try:
                stderrSubmitContents = _read_file(f"{stderr_file}.submit")
            except Exception:  # noqa
                pass  # submit stderr file doesn't always exist
            else:
                result["stderrSubmitContents"] = stderrSubmitContents
        if "stdout" in self.data:
            # include *contents* of stdout file, instead of file path
            stdout_file = self.data["stdout"]
            result["stdoutContents"] = _read_file(stdout_file)
        return result

    def set_real_time_status(self) -> None:
        """
        Distinguish between "Queued" and "Running" states.
        Check the stderr mtime to see when the task started running.
        This is necessary because the executionEvents are not populated until the
        task completes execution.
        This is only used by the summary() method because a) we usually don't need
        real-time data and b) accessing the file system could be unnecessarily slow.
        """
        if self.execution_status != "Running":
            pass
        elif self.stderr is None or self.queue_start is None:
            self.execution_status = "Queued"
        elif self.stderr.startswith("s3://"):
            # we don't have access to the EBS volume; stderr will not be copied to S3 until after the
            # task completes
            pass
        else:
            try:
                stderr_mtime = os.path.getmtime(self.stderr)
            except Exception:  # noqa
                self.execution_status = "Queued"
            else:
                # task is actually "Running" so calculate the queue-wait duration
                self.run_start = datetime.fromtimestamp(
                    stderr_mtime, tz=timezone.utc
                ).strftime("%Y-%m-%d %H:%M:%S")
                delta = parser.parse(self.run_start) - parser.parse(self.queue_start)
                self.queue_duration = str(delta)


class TaskError(Exception):
    pass


class Task:
    """
    A Task may have multiple calls, corresponding to multiple attempts and/or shards.
    If a Task is a subworkflow, it's Call object will contain a Metadata object.
    """

    def __init__(self, name, data):
        """
        Initialize a Task object, which may contain multiple shards or a subworkflow.

        :param name: Task name
        :type name: str
        :param data: Cromwell's "calls" for a task; this always remains unchanged.
        :type data: list
        """
        self.name = name
        self.data = data
        self.calls = {}
        self.subworkflows = {}
        for call_data in data:
            # shardIndex is "-1" for regular (not scattered) tasks
            shard_index = int(call_data["shardIndex"])
            # in general, there is only 1 attempt
            attempt = int(call_data["attempt"])

            if "subWorkflowMetadata" in call_data:
                sub_data = call_data["subWorkflowMetadata"]
                if shard_index not in self.subworkflows:
                    self.subworkflows[shard_index] = {}
                self.subworkflows[shard_index][attempt] = Metadata(sub_data)
            else:
                if shard_index not in self.calls:
                    self.calls[shard_index] = {}
                self.calls[shard_index][attempt] = Call(call_data, self.name)

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
                sub_errors = {
                    "shardIndex": shard_index,
                    "attempt": attempt,
                    "subWorkflowMetadata": sub_meta.errors(),
                }
                all_errors.append(sub_errors)
        return all_errors

    def summary(self, **kwargs):
        """
        Generally, there is only one attempt.  The only exception is when the resubmit with more memory
        feature of Cromwell is turned on, but it's only available with "gcs" backend as of 5/25/2022.
        For simplicity, we use only the last attempt by default.
        """
        real_time = True
        if "real_time" in kwargs and kwargs["real_time"] is False:
            real_time = False
        last_attempts = False
        if "last_attempts" in kwargs and kwargs["last_attempts"] is True:
            last_attempts = True
        if last_attempts is True:
            return self.summary_last_attempts(real_time)
        else:
            return self.summary_all_attempts(real_time)

    def summary_all_attempts(self, real_time):
        result = []
        for shard_index in self.calls.keys():
            for attempt in self.calls[shard_index].keys():
                call = self.calls[shard_index][attempt]
                result.append(call.summary(real_time=real_time))
        for shard_index in self.subworkflows.keys():
            for attempt in self.subworkflows[shard_index].keys():
                sub_meta = self.subworkflows[shard_index][attempt]
                for item in sub_meta.task_summary(real_time=real_time):
                    renamed_item = item
                    name = item["name"]
                    renamed_item["name"] = f"{self.name}:{name}"
                    result.append(renamed_item)
        return result

    def summary_last_attempts(self, real_time):
        result = []
        for shard_index in self.calls.keys():
            attempts = sorted(self.calls[shard_index].keys())
            attempt = attempts[-1]
            call = self.calls[shard_index][attempt]
            result.append(call.summary(real_time=real_time))
        for shard_index in self.subworkflows.keys():
            attempts = sorted(self.subworkflows[shard_index].keys())
            attempt = attempts[-1]
            sub_meta = self.subworkflows[shard_index][attempt]
            for item in sub_meta.task_summary(real_time=real_time):
                renamed_item = item
                name = item["name"]
                renamed_item["name"] = f"{self.name}:{name}"
                result.append(renamed_item)
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

    def __init__(self, data):
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
        self.workflow_id = data["id"]
        self.tasks = {}  # Task objects
        if "calls" in self.data:
            calls = self.data["calls"]
            for task_name, task_data in calls.items():
                logger.debug(f"Workflow {self.workflow_id}: Init task {task_name}")
                self.tasks[task_name] = Task(task_name, task_data)

    def get(self, param, default=None):
        """
        Get a section of the document.
        :param param: key of parameter to get
        :type param: str
        :param default: Default value if param undefined
        """
        return self.data.get(param, default)

    def workflow_name(self):
        return self.get("workflowName", "unknown")

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
        workflow_name = self.get("workflowName", None)
        workflow_id = self.get("id", None)
        if "executions_dir" in kwargs:
            executions_dir = kwargs["executions_dir"]
            workflow_root = f"{executions_dir}/{workflow_name}/{workflow_id}"
            return workflow_root

        # for AWS, if executions dir was not provided, try to infer from an outfile
        outputs = self.get("outputs", {})
        an_outfile = None
        for key, value in outputs.items():
            if value is None:
                # skip if null (i.e. optional output was not produced)
                pass
            elif type(value) is list:
                # a sharded task may produce a list of outputs, one per shard
                for item in value:
                    if (
                        type(item) is str
                        and item is not None
                        and self.workflow_id in item
                    ):
                        an_outfile = item
                        break
            elif value is not None and type(value) is str and self.workflow_id in value:
                # a typical task produces outputs which may be a file path
                an_outfile = value
                break
        if an_outfile is not None:
            (root, partial_path) = an_outfile.split(self.workflow_id)
            workflow_root = os.path.join(root, self.workflow_id)
            return workflow_root

        # for AWS, if the run doesn't have an outfile, return path relative to cromwell-execution dir
        workflow_root = f"{workflow_name}/{workflow_id}"
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
        Return list of all tasks, including any subworkflows.
        :param real_time: Check file system to distinguish between Queued and Running states
        :ptype real_time: bool
        :return: List of task information dictionaries
        :rtype: list
        """
        summary = []
        for task_name, task in self.tasks.items():
            for item in task.summary(**kwargs):
                summary.append(item)
        return summary

    def task_log(self, **kwargs):
        """
        Return select task summary fields in table format.
        :param real_time: Check file system to distinguish between Queued and Running states
        :ptype real_time: bool
        :return: Table of tasks
        :rtype: list
        """
        table = []
        for info in self.task_summary(**kwargs):
            name = info["name"]
            if info["shard_index"] > 0:
                name = f"{name}[{info['shard_index']}]"
            row = [
                name,
                info["cached"],
                info["execution_status"],
                info["queue_start"],
                info["run_start"],
                info["run_end"],
                info["queue_duration"],
                info["run_duration"],
            ]
            table.append(row)
        return table

    def started_running(self):
        """
        :return: True if any task actually started running; false otherwise.
        :rtype: bool
        """
        for info in self.task_summary(real_time=True):
            if info["run_start"] is not None:
                return True
        return False

    def job_summary(self, **kwargs) -> dict:
        """
        :return: task name for each job id
        :rtype: dict
        """
        job_ids = {}
        for task in self.task_summary(real_time=False):
            if "job_id" in task:
                job_id = str(task["job_id"])
                job_ids[job_id] = task.get("name", "unknown")
        return job_ids

    def outputs(self, relpath=True):
        """Returns all outputs for a workflow"""
        outputs = self.get("outputs", {})
        workflowRoot = self.workflow_root()
        if relpath and workflowRoot:
            relpath_outputs = {}
            for key, value in outputs.items():
                if value is None:
                    # skip if null (i.e. optional output was not produced)
                    pass
                elif type(value) is list:
                    # a sharded task may produce a list of outputs, one per shard
                    relpath_outputs[key] = []
                    for item in value:
                        if type(item) is str and item is not None:
                            relpath_outputs[key].append(
                                item.replace(workflowRoot, ".", 1)
                            )
                elif value is not None and type(value) is str:
                    # a typical task produces outputs which may be a file path
                    relpath_outputs[key] = value.replace(workflowRoot, ".", 1)
            outputs = relpath_outputs
        return outputs

    def outfiles(self, complete=False, relpath=True):
        """
        Return list of all output files of a run.
        By default, only include files tagged as outputs for the Run.
        :param complete: All files, not just workflow outputs.
        :ptype complete: bool
        :param relpath: Convert abs paths to rel paths.
        :ptype relpath: bool
        :return: List of files
        :rtype: list
        """
        workflow_root = self.workflow_root()
        outputs = self.get("outputs", {})
        if workflow_root is None:
            return []
        elif complete is True:
            if relpath:
                return []
            else:
                return [workflow_root]
        elif len(outputs.keys()) == 0:
            return []
        full_paths = []
        for key, value in outputs.items():
            if value is None:
                # skip if null (i.e. optional output was not produced)
                pass
            elif type(value) is list:
                # a sharded task may produce a list of outputs, one per shard
                for item in value:
                    if (
                        type(item) is str
                        and item is not None
                        and item.startswith(workflow_root)
                    ):
                        full_paths.append(item)
            elif (
                value is not None
                and type(value) is str
                and value.startswith(workflow_root)
            ):
                # a typical task produces outputs which may be a file path
                full_paths.append(value)

        if relpath:
            rel_paths = []
            for path in full_paths:
                rel_paths.append(path.replace(workflow_root, ".", 1))
            return rel_paths
        else:
            return full_paths


class Cromwell:
    """Class representing a Cromwell REST server."""

    workflows = {}

    def __init__(self, url: str) -> None:
        """Init object with connection info.

        :param url: url of Cromwell REST server, including port
        :type param: str
        """
        if not url.startswith("http"):
            url = f"http://{url}"
        self.url = url
        self.workflows_url = f"{url}/api/workflows/v1"
        self.engine_url = f"{url}/engine/v1/status"

    def get_metadata(self, workflow_id: str, data=None):
        """Get Metadata object for a workflow-run.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :param data: optionally provide the metadata instead
        :type data: dict
        :return: Metadata object
        :rtype: cromwell.Metadata
        """
        url = f"{self.workflows_url}/{workflow_id}/metadata?expandSubWorkflows=1"
        try:
            response = requests.get(url)
        except requests.ConnectionError as error:
            raise error
        response.raise_for_status()
        data = response.json()
        return Metadata(data)

    def status(self):
        """Check if Cromwell is available"""
        try:
            response = requests.get(self.engine_url)
        except Exception as error:
            raise error
        response.raise_for_status()
        return True

    def abort(self, workflow_id: str):
        """Abort a run."""
        url = f"{self.workflows_url}/{workflow_id}/abort"
        try:
            response = requests.post(url)
        except Exception as error:
            raise error
        sc = response.status_code
        if sc == 200:
            return
        elif sc == 400:
            raise ValueError(f"Abort failed for malformed workflow id: {workflow_id}")
        elif sc == 403:
            return  # too late to cancel
        elif sc == 404:
            raise ValueError(f"Cannot abort; workflow not found: {workflow_id}")
        else:
            return

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
        logger = logging.getLogger(__package__)
        files = {}
        if "wdl" in file_handles:
            files["workflowSource"] = (
                "workflowSource",
                file_handles["wdl"],
                "application/json",
            )
        else:
            raise ValueError("WDL file handle required")
        if "inputs" in file_handles:
            files["workflowInputs"] = (
                "workflowInputs",
                file_handles["inputs"],
                "application/json",
            )
        else:
            raise ValueError("Inputs JSON file handle required")
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
            response = requests.post(self.workflows_url, files=files)
        except requests.ConnectionError as error:
            lines = f"{error}".splitlines()
            logger.exception(
                f"Error submitting new run: Cromwell unavailable: {lines[-1]}"
            )
            raise
        except Exception as error:
            logger.exception(f"Error submitting new run: {error}")
            raise error
        response.raise_for_status()
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
        url = f"{self.workflows_url}/{workflow_id}/status"
        try:
            response = requests.get(url)
        except Exception as error:
            raise error
        response.raise_for_status()
        result = response.json()
        return result["status"]


# TODO: Moved from deprecated tasks.py but this code has been replaced by the function below
#    def get_task_cromwell_dir_mapping(self):
#        """
#        Return a dictionary that maps the cromwell directory to the task name for each task in the run.
#        The key=cromwell_dir, the value=task_name.
#        The cromwell directory name is modified to remove the root dir.
#        Ex: /my/path/cromwell-executions/task1/execution becomes
#            cromwell-executions/task1/execution
#        """
#        tasks = self.metadata.tasks
#        cromwell_to_task_names = {}
#
#        for name in tasks:
#            for entry in tasks[name].summary():
#
#                # remove root path from cromwell-exections dir.
#                # ex: /my/path/cromwell-executions/a/b/c is converted to cromwell-executions/a/b/c
#                if len(entry) < 5:
#                    continue
#
#                cromwell_dir = entry[4]
#                idx = cromwell_dir.find("cromwell-executions")
#
#                # if cromwell-executions not found in directory name, idx=-1. skip this condition.
#                if idx >= 0:
#                    cromwell_dir = cromwell_dir[idx:]
#                    entry[4] = cromwell_dir
#                    cromwell_to_task_names[cromwell_dir] = entry[0]
#        return cromwell_to_task_names


def parse_cromwell_task_dir(task_dir):
    """
    Given path of a task, return task fields.
    """
    result = {
        "call_root": task_dir,
        "cached": False,
        "shard": -1,
    }
    if task_dir.endswith("/execution"):
        result["call_root"] = result["call_root"].rstrip("/execution")
    else:
        task_dir = f"{task_dir}/execution"
    (root_dir, subdir) = task_dir.split("cromwell-executions/")
    result["call_root_rel_path"] = subdir
    fields = deque(subdir.split("/"))
    result["wdl_name"] = fields.popleft()
    result["name"] = result["wdl_name"]
    result["cromwell_run_id"] = fields.popleft()
    if not fields[0].startswith("call-"):
        raise ValueError(f"Problem parsing {subdir}")
    result["task_name"] = fields.popleft().lstrip("call-")
    result["name"] = f"{result['name']}.{result['task_name']}"
    if fields[0].startswith("shard-"):
        result["shard"] = int(fields.popleft().lstrip("shard-"))
        result["name"] = f"{result['name']}[{result['shard']}]"
    if fields[0] == "execution":
        return result
    elif fields[0] == "cacheCopy":
        result["cached"] = True
        return result

    # subworkflow
    result["subworkflow_name"] = fields.popleft()
    if not result["subworkflow_name"].startswith("sub."):
        raise ValueError(f"Problem parsing {subdir}")
    result["subworkflow_name"] = result["subworkflow_name"].lstrip("sub.")
    result["name"] = f"{result['name']}:{result['subworkflow_name']}"
    result["subworkflow_cromwell_run_id"] = fields.popleft()
    if not fields[0].startswith("call-"):
        raise ValueError(f"Problem parsing {subdir}")
    result["sub_task_name"] = fields.popleft().lstrip("call-")
    result["name"] = f"{result['name']}.{result['sub_task_name']}"
    if fields[0].startswith("shard-"):
        result["sub_shard"] = int(fields.popleft().lstrip("shard-"))
        result["name"] = f"{result['name']}[{result['sub_shard']}]"
    if fields[0] == "execution":
        return result
    elif fields[0] == "cacheCopy":
        result["cached"] = True
        return result
    else:
        raise ValueError(f"Problem parsing {subdir}")
    return result
