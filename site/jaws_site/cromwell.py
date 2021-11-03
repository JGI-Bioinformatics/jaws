"""Cromwell class provides OO-interface to Cromwell REST server.

NOTE: This class employs the verbiage employed by Cromwell, which may differ slightly from that used
elsewhere in JAWS/JTM, as clarified below:

- Cromwell "workflow": the execution of a "Workflow Specification" on a particular set of input.
  - "Workflow Specification" is colloquially referred to as a "WDL file" or "widdle"
  - a Cromwell "Workflow Specification" is called a "Workflow" in JAWS
  - a Cromwell "Workflow" is referred to as a "Run" in JAWS
  - the JAWS "cromwell_run_id" is equivalent to the Cromwell "workflow_id"

- Cromwell "task": a specific step in a Workflow Specification
  - same in JAWS, which uses the variable "cromwell_task_name" for Cromwell's "task_name"
  - no equivalent in JTM
  - JAWS uses Cromwell's task name, in the format: "<WDL_name>.<task_name>"

- Cromwell "job": the execution of a Task
  - a Task may be atempted more than once, each with a unique "attempt" (int) and "job_id"
  - same in JAWS, which uses "cromwell_job_id" to refer to Cromwell's "job_id";
    while JTM calls this "task_id"
  - for JTM a "job" is a compute reservation and is equivalent to the scheduler job,
    so the JTM "job_id" will match that of SLURM's `squeue` or SGE's `qstat` output.
    A single JTM job (running a JTM Worker) executes many JTM tasks (Cromwell jobs).
"""


import requests
import logging
import os
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
    if path and os.path.isfile(path):
        with open(path, "r") as file:
            contents = file.read()
    return contents


class Task:
    """
    A Task may have multiple calls, corresponding to multiple attempts and/or shards.
    If a Task is a subworkflow, it will contain a Metadata object.
    """

    def __init__(self, name, data):
        """
        Initialize a Task object, which may contain a subworkflow.
        A Task may contain multiple calls corresponding to multiple shards, but only the
        last attempt (of each) is represented by the object.

        :param name: Task name
        :type name: str
        :param data: Cromwell's "calls" for a task; this always remains unchanged.
        :type data: list
        """
        self.name = name
        self.data = data
        self.subworkflows = {}
        self.max_shard_index = None
        self.num_attempts = {}  # max attempt per shard_index

        if len(data) == 0:
            return

        for call in data:
            shard_index = call["shardIndex"]
            attempt = call["attempt"]

            if self.max_shard_index is None or self.max_shard_index < shard_index:
                self.max_shard_index = shard_index

            if (
                shard_index not in self.num_attempts
                or self.num_attempts[shard_index] < attempt
            ):
                self.num_attempts[shard_index] = attempt

            if "subWorkflowMetadata" in call:
                # this task is a subworkflow, with it's own tasks
                if shard_index not in self.subworkflows:
                    self.subworkflows[shard_index] = {}
                self.subworkflows[shard_index][attempt] = Metadata(
                    call["subWorkflowMetadata"]
                )

    def get(self, key, shard_index: int = -1, attempt: int = 1, default=None):
        """
        Get an item from the call dictionary (e.g. "executionStatus")
        corresponding to the last attempt.
        :param key: Dictionary key (e.g. "executionStatus")
        :type key: str
        :param default: Default value to return if key not found (default=None).
        :return: value of the specified key, for specified attempt number.
        """
        if shard_index > self.max_shard_index:
            raise ValueError(
                f"ShardIndex {shard_index} is invalid for Task {self.name} (max. shardIndex is {self.max_shard_index}"
            )
        if (
            shard_index not in self.num_attempts
            or attempt > self.num_attempts[shard_index]
        ):
            raise ValueError(
                f"Attempt {attempt} is invalid for Task {self.name} shardIndex {shard_index}"
            )
        for call in self.data:
            if call["shardIndex"] == shard_index and call["attempt"] == attempt:
                return call.get(key, default)
        return default

    def summary(self):
        """
        Summary of calls, their jobIds, and whether or not they are cached results.
        """
        summary = []
        for call in self.data:
            name = self.name
            shard_index = call["shardIndex"]
            if shard_index > -1:
                name = f"{name}[{shard_index}]"
            attempt = call["attempt"]
            if self.num_attempts[shard_index] > 1:
                name = f"{name}.{attempt}"
            cached = False
            if "callCaching" in call and "hit" in call["callCaching"]:
                cached = call["callCaching"]["hit"]
            job_id = None
            if "jobId" in call:
                job_id = call["jobId"]
            run_time = None
            if "executionEvents" in call:
                for event in call["executionEvents"]:
                    if (
                        event["description"] == "RunningJob"
                        and "startTime" in event
                        and "endTime" in event
                    ):
                        start_time = parser.parse(event["startTime"])
                        end_time = parser.parse(event["endTime"])
                        delta = end_time - start_time
                        run_time = str(delta)
            if "subWorkflowMetadata" in call:
                subworkflow = self.subworkflows[shard_index][attempt]
                sub_task_summary = subworkflow.task_summary()
                for sub_name, sub_job_id, sub_cached, run_time in sub_task_summary:
                    summary.append(
                        [f"{name}:{sub_name}", sub_job_id, sub_cached, run_time]
                    )
            else:
                summary.append([name, job_id, cached, run_time])
        return summary

    def errors(self):
        """
        Return user friendly errors report for this task/subworkflow.
        This is a copy of the call data with only pertinent elements of the last
        attempt included, for brevity and readability.
        The contents of the stderr and stderr.submit files are also added.
        :return: Errors report (filtered task call data)
        :rtype: dict
        """
        filteredCalls = []
        for call in self.data:
            shard_index = call["shardIndex"]
            attempt = call["attempt"]
            status = call["executionStatus"]
            if status != "Failed":
                # There is only error information if this task actually failed
                continue
            filteredCall = {
                "shardIndex": shard_index,
                "attempt": attempt,
            }
            if "subWorkflowMetadata" in call:
                # include any errors from the subworkflow's tasks
                subworkflow = self.subworkflows[shard_index][attempt]
                filteredCall["subWorkflowMetadata"] = subworkflow.errors()
            else:
                # simple task (not a subworkflow)
                if "failures" in call:
                    filteredCall["failures"] = call["failures"]
                if "jobId" in call:
                    filteredCall["jobId"] = call["jobId"]
                if "returnCode" in call:
                    filteredCall["returnCode"] = call["returnCode"]
                if "runtimeAttributes" in call:
                    filteredCall["runtimeAttributes"] = call["runtimeAttributes"]
                if "stderr" in call:
                    # include *contents* of stderr files, instead of file paths
                    stderr_file = call["stderr"]
                    filteredCall["stderrContents"] = _read_file(stderr_file)
                    filteredCall["stderrSubmitContents"] = _read_file(
                        f"{stderr_file}.submit"
                    )
            filteredCalls.append(filteredCall)
        return filteredCalls

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
            raise ValueError(
                f"Invalid file id, {file_id}; allowed values: stdout, stderr"
            )
        path = self.get(file_id)
        if relpath:
            call_root = self.get("callRoot")
            path = os.path.relpath(call_root, path)
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

    def workflow_root(self):
        return self.get("workflowRoot", None)

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
            for name, job_id, cached, run_time in task_summary:
                summary.append([name, job_id, cached, run_time])
        return summary

    def outputs(self, **kwargs):
        """Returns all outputs for a workflow"""
        outputs = self.get("outputs", {})
        if "relpath" in kwargs and kwargs["relpath"] is True:
            workflowRoot = self.get("workflowRoot")
            relpath_outputs = {}
            for key, value in outputs.items():
                if value is None:
                    # skip if null (i.e. optional output was not produced)
                    pass
                elif type(value) is list:
                    # a sharded task may produce a list of outputs, one per shard
                    relpath_outputs[key] = []
                    for item in value:
                        if type(value) is str and item is not None:
                            relpath_outputs[key].append(
                                item.replace(workflowRoot, ".", 1)
                            )
                elif type(value) is str:
                    # a typical task produces outputs which may be a file path
                    relpath_outputs[key] = value.replace(workflowRoot, ".", 1)
            outputs = relpath_outputs
        return outputs


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
        response.raise_for_status()

    def submit(
        self,
        wdl_file: str,
        json_file: str,
        zip_file: str = None,
        options_file: str = None,
    ) -> int:
        """
        Submit a run to Cromwell.
        :param wdl_file: Path to WDL file
        :type wdl_file: str
        :param json_file: Path to inputs JSON file
        :type json_file: str
        :param zip_file: Path to subworkflows ZIP file (optional)
        :type zip_file: str
        :param options_file: Path to options JSON file (optional)
        :type options_file: str
        :return: Cromwell workflow uuid
        :rtype: str
        """
        logger = logging.getLogger(__package__)
        files = {}
        try:
            files["workflowInputs"] = (
                "workflowInputs",
                open(json_file, "r"),
                "application/json",
            )
        except Exception as error:
            logger.exception(f"Unable to open file, {json_file}: {error}")
            raise IOError(f"Unable to open file, {json_file}: {error}")
        try:
            files["workflowSource"] = (
                "workflowSource",
                open(wdl_file, "r"),
                "application/json",
            )
        except Exception as error:
            logger.exception(f"Unable to open file, {wdl_file}: {error}")
            raise IOError(f"Unable to open file, {wdl_file}: {error}")
        if zip_file:
            try:
                files["workflowDependencies"] = (
                    "workflowDependencies",
                    open(zip_file, "rb"),
                    "application/zip",
                )
            except Exception as error:
                raise IOError(f"Unable to open file, {zip_file}: {error}")
        if options_file:
            try:
                files["workflowOptions"] = (
                    "workflowOptions",
                    open(options_file, "r"),
                    "application/json",
                )
            except Exception as error:
                logger.exception(f"Unable to open file, {options_file}: {error}")
                raise IOError(f"Unable to open file, {options_file}: {error}")
        try:
            response = requests.post(self.workflows_url, files=files)
        except requests.ConnectionError as error:
            lines = f"{error}".splitlines()
            logger.exception(
                f"Error submitting new run: Cromwell unavailable: {lines[-1]}"
            )
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
