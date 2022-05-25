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
import re
import json
import io
from dateutil import parser


def _read_file(path: str):
    """
    Read a file and return it's contents.
    :param path: Path to file
    :type path: str
    :return: Contents of the file if it exists, else None.
    :rtype: str
    """
    contents = None
    # TODO special handling required for files in S3  (see example in runs.py; move to common py)
    if path and os.path.isfile(path):
        with open(path, "r") as file:
            contents = file.read()
    return contents


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
        self.task_name = task_name
        self.name = task_name
        self.attempt = data.get("attempt", None)
        self.shard_index = data.get("shardIndex", None)
        self.execution_status = data.get("executionStatus", None)
        self.cached = None
        self.return_code = data.get("returnCode", None)
        self.start = data.get("start", None)
        self.end = data.get("end", None)
        self.stdout = self.data.get("stdout", None)
        self.stderr = self.data.get("stderr", None)
        self.call_root = self.data.get("callRoot", None)
        self.execution_dir = None
        self.job_id = self.data.get("jobId", None)
        self.max_time = None
        self.memory = None
        self.cpu = None

        if self.stdout is not None:
            self.execution_dir = re.sub(r"/stdout$", "", self.stdout)
        elif self.call_root is not None:
            self.execution_dir = f"{self.call_root}/execution"

        if "runtimeAttributes" in self.data:
            self.max_time = self.data["runtimeAttributes"].get("time", None)
            self.memory = self.data["runtimeAttributes"].get("memory", None)
            self.cpu = self.data["runtimeAttributes"].get("cpu", None)

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
                    self.queue_start = event["startTime"]
                elif event["description"] == "RunningJob":
                    self.run_start = event["startTime"]
                elif event["description"] == "CallCacheReading":
                    self.run_start = event["startTime"]
                elif event["description"] == "UpdatingJobStore":
                    self.run_end = event["startTime"]
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

    def summary(self, fmt="dict"):
        if fmt == "dict":
            return self.as_dict()
        elif fmt == "list":
            return self.as_list()
        else:
            raise ValueError(f"Allowed formats are 'list', 'dict'; got: ${fmt}")

    def as_dict(self):
        """
        :return: Select fields
        :rtype: dict
        """
        record = {
            "name": self.name,
            "attempt": self.attempt,
            "cached": self.cached,
            "execution_status": self.execution_status,
            "queue_start": self.queue_start,
            "run_start": self.run_start,
            "run_end": self.run_end,
            "queue_duration": self.queue_duration,
            "run_duration": self.run_duration,
            "execution_dir": self.execution_dir,
        }
        return record

    def as_list(self):
        """
        :return: Select fields
        :rtype: list
        """
        row = [
            self.name,
            self.attempt,
            self.cached,
            self.execution_status,
            self.queue_start,
            self.run_start,
            self.run_end,
            self.queue_duration,
            self.run_duration,
            self.execution_dir,
        ]
        return row

    def logs(self, fmt="list"):
        """
        :return: logs (one if task, many if subworkflow)
        :rtype: list
        """
        if self.subworkflow is None:
            return self.summary(fmt)

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
        if self.subworkflow is None:
            # simple task (not a subworkflow)
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
                result["stderrSubmitContents"] = _read_file(f"{stderr_file}.submit")
            if "stdout" in self.data:
                # include *contents* of stdout file, instead of file path
                stdout_file = self.data["stdout"]
                result["stdoutContents"] = _read_file(stdout_file)
        else:
            result["subWorkflowMetadata"] = self.subworkflow.errors()
        return result


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

    def logs(self):
        """
        Return list of all call logs (i.e. all shards, all subworkflow calls).
        :return: All call logs
        :rtype: list
        """
        all_logs = []
        for shard_index in self.calls.keys():
            for attempt in self.calls[shard_index].keys():
                call = self.calls[shard_index][attempt]
                all_logs.append(call.log())
        return all_logs
#    def logs(self, fmt="list"):
#        """
#        :return: logs (one if task, many if subworkflow)
#        :rtype: list
#        """
#        if self.subworkflow is None:
#            return [self.summary(fmt)]
#        else:
#            logs = []
#            sub_logs = self.subworkflow.logs(fmt=fmt)
#            for log in sub_logs:
#                if fmt == "dict":
#                    log["name"] = f"{self.name}:{log['name']}"
#                else:
#                    log[0] = f"{self.name}:{log[0]}"
#                logs.append(log)
#            return logs

    def errors(self):
        """
        Return user friendly errors report for this task/subworkflow.
        This is a copy of the call data with only pertinent elements of the last
        attempt included, for brevity and readability.
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
        return all_errors

    def summary(self):
        result = []
        for shard_index in self.calls.keys():
            # include only last attempt in summary
            attempt = sorted(self.calls[shard_index].keys())[-1]
            call = self.calls[shard_index][attempt]
            row = [
                call.name,
                call.job_id,
                call.cached,
                call.max_time,
                call.execution_status,
                call.call_root
            ]
            result.append(row)
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

    def workflow_root(self):
        """
        Return the base location for the Run.
        Note that workflowRoot is only defined when using NFS;
        for S3 we must infer from an output file.
        """
        workflow_root = self.get("workflowRoot", None)
        if workflow_root:
            return workflow_root

        # extract from an outfile or return empty string if no outfiles.
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
        if not an_outfile:
            return None
        (root, partial_path) = an_outfile.split(self.workflow_id)
        workflow_root = os.path.join(root, self.workflow_id)
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

    def task_summary(self):
        """
        Return table of all tasks, including any subworkflows.
        """
        summary = []
        for task_name, task in self.tasks.items():
            task_summary = task.summary()
            for (
                name,
                job_id,
                cached,
                max_time,
                execution_status,
                cromwell_dir,
            ) in task_summary:
                summary.append(
                    [name, job_id, cached, max_time, execution_status, cromwell_dir]
                )
        return summary

    def job_summary(self):
        """Return task info, organized by job_id."""
        task_summary = self.task_summary()
        job_summary = {}
        for (
            task_name,
            job_id,
            cached,
            max_time,
            execution_status,
            _,
        ) in task_summary:
            if job_id:
                job_id = str(job_id)
                job_summary[job_id] = [task_name, max_time]
        return job_summary

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


# TODO: Moved from deprecated tasks.py but this code seems wrong.
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
