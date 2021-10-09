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

    def __init__(self, workflows_url, name, calls, cache={}):
        """
        Initialize a Task object, which may contain a subworkflow.

        :param workflows_url: The Cromwell URL to GET workflow metadata
        :type workflows_url: str
        :param name: Task name
        :type name: str
        :param calls: A task's calls (AKA attempts) section of the Cromwell metadata.
        :type calls: list
        :param cache: cached {workflow_id => metadata}; avoids need to GET from Cromwell (optional)
        :type cache: dict
        """
        logger = logging.getLogger(__package__)
        self.workflows_url = workflows_url
        self.name = name
        self.calls = calls
        self.subworkflows = {}  # subworkflow_id => Metadata obj

        for call in calls:
            if "subWorkflowId" in call:
                workflow_id = call["subWorkflowId"]
                logger.debug(
                    f"Task {self.name} is a subworkflow; getting metadata {workflow_id}"
                )
                metadata = Metadata(self.workflows_url, workflow_id, None, cache)
                # a subworkflow task has a workflow id for each attempt
                self.subworkflows[workflow_id] = metadata

    def is_subworkflow(self):
        return bool(self.subworkflows)

    def get(self, key, attempt=None, default=None):
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

    def call_root_dir(self, attempt=None):
        return self.get("callRoot", attempt)

    def execution_status(self, attempt=None):
        return self.get("executionStatus", attempt)

    def errors(self):
        """
        Return user friendly errors report for this task.
        :return: Errors report
        :rtype: dict
        """
        full_report = []
        for call in self.calls:
            if call.get("executionStatus") == "Failed":
                report = {}
                report["failures"] = call["failures"]
                if "jobId" in call:
                    report["jobId"] = call["jobId"]
                if "runtimeAttributes" in call:
                    report["runtimeAttributes"] = {}
                    runtime_attr = call["runtimeAttributes"]
                    if "memory" in runtime_attr:
                        report["runtimeAttributes"]["memory"] = runtime_attr["memory"]
                    if "cpu" in runtime_attr:
                        report["runtimeAttributes"]["cpu"] = runtime_attr["cpu"]
                    if "time" in runtime_attr:
                        report["runtimeAttributes"]["time"] = runtime_attr["time"]
                if "stderr" in call:
                    stderr_file = call["stderr"]
                    if os.path.isfile(stderr_file):
                        with open(stderr_file, "r") as file:
                            report["stderr"] = file.read()
                    stderr_submit_file = f"{stderr_file}.submit"
                    if os.path.isfile(stderr_submit_file):
                        with open(stderr_submit_file, "r") as file:
                            report["stderr.submit"] = file.read()
                full_report.append(report)
        return full_report

    def stdout(self, attempt=None, src=None, dest=None):
        """
        Return the path to the standard output file, optionally replacing part of the path.
        :param attempt: attempt number (first is 1; default is last)
        :type attempt: int
        :param src: Source root dir
        :type src: str
        :param dest: Destination root dir
        :type dest: str
        :return: Path to stdout file
        :rtype: str
        """
        path = self.get("stdout", attempt)
        if not path:
            return None
        if src and dest:
            path = os.path.join(dest, os.path.relpath(src, path))
        return path

    def stderr(self, attempt=None, src=None, dest=None):
        """
        Return the path to the standard err file, optionally replacing part of the path.
        :param attempt: attempt number (first is 1; default is last)
        :type attempt: int
        :param src: Source root dir
        :type src: str
        :param dest: Destination root dir
        :type dest: str
        :return: Path to stdout file
        :rtype: str
        """
        path = self.get("stderr", attempt)
        if not path:
            return None
        if src and dest:
            path = os.path.join(dest, os.path.relpath(src, path))
        return path

    def runtime(self, attempt=None):
        return self.get("runtimeAttributes", attempt, {})


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

    def __init__(self, workflows_url, workflow_id, data=None, cache={}):
        """
        Initialize Metadata object by retrieving metadata from Cromwell,
        including subworkflows, and initialize Task objects.
        :param workflows_url: URL for Cromwell
        :type workflows_url: str
        :param workflow_id: Cromwell's UUID for the workflow (AKA run)
        :type workflow_id: str
        :param data: optionally provide metadata JSON to avoid GET from Cromwell
        :type dict:
        :param cache: optionally provide multiple metadata to avoid GET from Cromwell
        :type cache: dict
        """
        logger = logging.getLogger(__package__)
        self.workflows_url = workflows_url
        self.workflow_id = workflow_id
        self.tasks = None
        self.data = None
        if data:
            self.data = data
        elif cache and workflow_id in cache:
            # cache dict uses workflow_id as keys and metadata as values;
            # it can be used to avoid unnecessary GETs and is currently
            # used by tests/test_cromwell.py
            self.data = cache[workflow_id]
        else:
            logger.debug(f"Get metadata for {workflow_id}")
            self._get_data()
        self._init_tasks(cache)

    def _get_data(self):
        """GET record from Cromwell REST server."""
        url = f"{self.workflows_url}/{self.workflow_id}/metadata"
        try:
            response = requests.get(url)
        except requests.ConnectionError as error:
            raise CromwellServiceConnectionError(f"Cromwell temporarily inaccessible: {error}")
        if response.status_code == 404:
            raise CromwellRunNotFoundError("Run not found")
        response.raise_for_status()
        self.data = response.json()

    def init_with_json_str(self, json_str):
        try:
            self.data = json.loads(json_str)
        except Exception as error:
            raise ValueError(f"Not valid json: {error}")

    def _init_tasks(self, cache):
        """Initialize and save Task objects.
        :param cache: Cached {workflow_id=>metadata} avoids unnecessary GET (optional)
        :type cache: dict
        """
        logger = logging.getLogger(__package__)
        self.tasks = []  # Task objects
        self.subworkflows = {}  # workflow_id => metadata obj
        if "calls" in self.data:
            calls = self.data["calls"]
            for task_name in calls.keys():
                logger.debug(f"Workflow {self.workflow_id}: Init task {task_name}")
                task = Task(self.workflows_url, task_name, calls[task_name], cache)
                self.tasks.append(task)
                if task.is_subworkflow:
                    for sub_id, sub_meta in task.subworkflows.items():
                        self.subworkflows[sub_id] = sub_meta

    def get(self, param, default=None):
        """Get a section of the document.
        :param param: key of parameter to get
        :type param: str
        :param default: Default value if param undefined
        """
        return self.data.get(param, default)

    def workflow_root(self):
        return self.get("workflowRoot", None)

    def is_subworkflow(self):
        return True if "parentWorkflowId" in self.data else False

    def execution_status(self):
        """
        Return dict of task name to execution status, for last attempt of each task.
        """
        result = {}
        for task in self.tasks:
            result[task.name] = task.execution_status()
        return result

    def failure_reason(self):
        """
        Return standard message of reason for failure, without detail.
        If there is more than one failure, only the first is returned.
        Example messages are:
        - "Workflow failed" (a Task had failde)
        - "Workflow input processing failed" (an infile was not found)
        """
        reason = None
        failures = self.get("failures")
        if failures:
            reason = failures[0]["message"]
        return reason

    def filtered_failures(self):
        """
        Filter failurs with message "Workflow failed", as those failures are duplicated in the Tasks.
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
        report = {}
        task_errors = {}
        for task in self.tasks:
            task_report = task.errors()
            if len(task_report):
                task_errors[task.name] = task_report
        if len(task_errors):
            report["calls"] = task_errors
        other_failures = self.filtered_failures()
        if len(other_failures):
            report["failures"] = self.filtered_failures()
            report["inputs"] = self.get("inputs")
        return report

    def task_summary(self):
        """
        Return table of all tasks, including any subworkflows.
        """
        summary = []
        for task in self.tasks:
            for call in task.calls:
                if "jobId" in call:
                    summary.append(
                        [self.workflow_id, task.name, call["attempt"], call["jobId"]]
                    )
                elif "subWorkflowId" in call:
                    subworkflow_id = call["subWorkflowId"]
                    subworkflow = task.subworkflows[subworkflow_id]
                    sub_summary = subworkflow.task_summary()
                    summary.extend(sub_summary)
        return summary

    def outputs(self, **kwargs):
        """Returns all outputs for a workflow"""
        outputs = self.get("outputs", {})
        if "relpath" in kwargs and kwargs["relpath"] is True:
            workflowRoot = self.get("workflowRoot")
            relpath_outputs = {}
            for key, value in outputs.items():
                if type(value) is list:
                    relpath_outputs[key] = []
                    for item in value:
                        relpath_outputs[key].append(item.replace(workflowRoot, '.', 1))
                else:
                    relpath_outputs[key] = value.replace(workflowRoot, '.', 1)
            outputs = relpath_outputs
        if "outfile" in kwargs:
            with open(kwargs["outfile"], 'w') as fh:
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

    def get_metadata(self, workflow_id: str, data=None, cache={}):
        """Get Metadata object for a workflow-run.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :param data: optionally provide the metadata instead
        :type data: dict
        :return: Metadata object
        :rtype: cromwell.Metadata
        """
        return Metadata(self.workflows_url, workflow_id, data, cache)

    def get_all_metadata(self, workflow_id: str, cache: dict = {}):
        """Get dict of all runs => metadata json for run and all subworkflows.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :return: all metadata docs { workflow_id => metadata obj }
        :rtype: dict
        """
        result = {}
        metadata = Metadata(self.workflows_url, workflow_id, None, cache)
        result[workflow_id] = metadata.data
        for sub_id, sub_meta in metadata.subworkflows.items():
            result[sub_id] = sub_meta.data
        return result

    def get_all_errors(self, workflow_id: str, cache: dict = {}):
        """Get dict of all runs => errors json for run and all subworkflows.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :return: all errors docs { workflow_id => errors obj }
        :rtype: dict
        """
        result = {}
        metadata = Metadata(self.workflows_url, workflow_id, None, cache)
        result[workflow_id] = metadata.errors()
        for sub_id, sub_meta in metadata.subworkflows.items():
            result[sub_id] = sub_meta.errors()
        return result

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
