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
import re
import json
import io


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
            if call.get("stdout"):
                cromwell_dir = re.sub(r"/stdout$", "", call["stdout"])
            else:
                cromwell_dir = None

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
            max_time = None
            if "runtimeAttributes" in call:
                max_time = call["runtimeAttributes"].get("time", "")
            if "subWorkflowMetadata" in call:
                subworkflow = self.subworkflows[shard_index][attempt]
                sub_task_summary = subworkflow.task_summary()
                for (
                    sub_name,
                    sub_job_id,
                    sub_cached,
                    max_time,
                    cromwell_dir,
                ) in sub_task_summary:
                    summary.append(
                        [
                            f"{name}:{sub_name}",
                            sub_job_id,
                            sub_cached,
                            max_time,
                            cromwell_dir,
                        ]
                    )
            else:
                summary.append([name, job_id, cached, max_time, cromwell_dir])
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
                if "stdout" in call:
                    # include *contents* of stdout file, instead of file path
                    stdout_file = call["stdout"]
                    filteredCall["stdoutContents"] = _read_file(stdout_file)
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
            for name, job_id, cached, max_time, cromwell_dir in task_summary:
                summary.append([name, job_id, cached, max_time, cromwell_dir])
        return summary

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
        full_paths = []
        outputs = self.get("outputs", {})
        if workflow_root is None:
            return full_paths
        elif complete is True:
            return [workflow_root]
        elif len(outputs.keys()) == 0:
            return full_paths
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
        response.raise_for_status()

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
