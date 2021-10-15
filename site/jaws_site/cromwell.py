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
import json
import logging
import os


class Task:
    """A Task may have multiple calls, each with a unique job_id.
    If a Task is a subworkflow, it will contain a Metadata object.
    """

    def __init__(self, name, calls):
        """
        Initialize a Task object, which may contain a subworkflow.

        :param name: Task name
        :type name: str
        :param calls: A task's calls (AKA attempts) section of the Cromwell metadata.
        :type calls: list
        """
        self.name = name
        self.calls = calls
        self.isSubWorkflow = False

        for call in calls:
            if "subWorkflowMetadata" in call:
                self.isSubWorkflow = True
                metadata = Metadata(call["subWorkflowMetadata"])
                call["subWorkflowMetadata"] = metadata

    def get(self, key, attempt=None, default=None):
        """
        Get an item from the list of calls dictionaries (e.g. "executionStatus");
        by default, get the item from the last attempt.
        :param key: Dictionary key (e.g. "executionStatus")
        :type key: str
        :param attempt: Attempt number (first=1; default=last).
        :type attempt: int
        :param default: Default value to return if key not found (default=None).
        :return: value of the specified key, for specified attempt number.
        """
        # first attempt is 1, convert to 0-based index if in valid range
        index = None
        if attempt is None:
            index = -1  # last attempt
        else:
            attempt = int(attempt)
            if attempt == 0 or attempt > len(self.calls):
                raise ValueError("Invalid attempt; of out range")
            else:
                index = attempt - 1
        return self.calls[index].get(key, default)

    def errors(self):
        """
        Return user friendly errors report for this task/subworkflow.
        This is a copy of the call data with only pertinent elements of the last
        attempt included, for brevity and readability.
        The contents of the stderr and stderr.submit files are also added.
        :return: Errors report (filtered task call data)
        :rtype: dict
        """
        if self.get("executionStatus") != "Failed":
            # There is only error information if this task actually failed
            return []
        filteredCall = {}
        call = self.calls[-1]  # only the last attempt is relevant
        if self.isSubWorkflow:
            # include any errors from the subworkflow's tasks
            filteredCall["subWorkflowMetadata"] = call["subWorkflowMetadata"].errors()
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
                stderrFilePath = call["stderr"]
                if os.path.isfile(stderrFilePath):
                    with open(stderrFilePath, "r") as file:
                        filteredCall["stderrContents"] = file.read()
                stderrSubmitFilePath = f"{stderrFilePath}.submit"
                if os.path.isfile(stderrSubmitFilePath):
                    with open(stderrSubmitFilePath, "r") as file:
                        filteredCall["stderrSubmitContents"] = file.read()
        # the result follows the structure of the original metadata,
        # so make a list corresponding to attempts
        filteredCalls = [filteredCall]
        return filteredCalls

    def stdout(self, attempt=None, relpath=False):
        """
        Return the path to the standard output file, optionally replacing part of the path.
        :param attempt: attempt number (first is 1; default is last)
        :type attempt: int
        :param relpath: If true then make path relative to callRoot; fullpath by default
        :type relpath: bool
        :param dest: Destination root dir
        :return: Path to stdout file
        :rtype: str
        """
        path = self.get("stdout", attempt)
        if relpath:
            call_root = self.get("callRoot", attempt)
            path = os.path.relpath(call_root, path)
        return path

    def stderr(self, attempt=None, relpath=False):
        """
        Return the path to the standard err file, optionally replacing part of the path.
        :param attempt: attempt number (first is 1; default is last)
        :type attempt: int
        :param relpath: If true then make path relative to callRoot; fullpath by default
        :type relpath: bool
        :return: Path to stdout file
        :rtype: str
        """
        path = self.get("stderr", attempt)
        if relpath:
            call_root = self.get("callRoot", attempt)
            path = os.path.relpath(call_root, path)
        return path


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
            for task_name in calls.keys():
                logger.debug(f"Workflow {self.workflow_id}: Init task {task_name}")
                self.tasks[task_name] = Task(task_name, calls[task_name])

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
            # filtered_metadata["inputs"] = self.get("inputs")
        return filtered_metadata

    def task_summary(self):
        """
        Return table of all tasks, including any subworkflows.
        """
        summary = []
        for task_name, task in self.tasks.items():
            if len(task.calls):
                # include last attempt only
                call = task.calls[-1]
                if "jobId" in call:
                    summary.append([task_name, call["jobId"]])
                elif "subWorkflowMetadata" in call:
                    subworkflow = call["subWorkflowMetadata"]
                    sub_summary = subworkflow.task_summary()
                    for (sub_task_name, job_id) in sub_summary:
                        summary.append([f"{task_name}:{sub_task_name}", job_id])
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
        if "outfile" in kwargs:
            with open(kwargs["outfile"], "w") as fh:
                fh.write(json.dumps(outputs, sort_keys=True, indent=4))
        else:
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
        url = f"{self.workflows_url}/{workflow_id}/metadata"
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
