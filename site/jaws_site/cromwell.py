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


class Task:
    """A Task may have multiple calls, each with a unique job_id.
    If a Task is a subworkflow, it will contain a Metadata object."""

    def __init__(self, workflows_url, name, calls):
        logger = logging.getLogger(__package__)
        self.workflows_url = workflows_url
        self.name = name
        self.calls = calls
        self.jobs = {}  # job_id => dict
        self.subworkflows = {}  # subworkflow_id => Metadata obj

        for call in calls:
            if "subWorkflowId" in call:
                workflow_id = call["subWorkflowId"]
                logger.debug(f"Task {self.name} is a subworkflow; getting metadata {workflow_id}")
                metadata = Metadata(self.workflows_url, workflow_id)
                self.subworkflows[workflow_id] = metadata
                # copy subworkflows' jobs
                for job_id in metadata.jobs:
                    logger.debug(f"Sub {self.name}, task {metadata.jobs[job_id]['task_name']}: job {job_id}")
                    self.jobs[job_id] = {
                        "task_name": f"{self.name}.{metadata.jobs[job_id]['task_name']}",
                        "attempt": metadata.jobs[job_id]["attempt"]
                    }
            elif "jobId" in call:
                job_id = int(call["jobId"])
                logger.debug(f"Task {self.name}: job {job_id}")
                self.jobs[job_id] = {
                    "task_name": self.name,
                    "attempt": int(call["attempt"])
                }

    def is_subworkflow(self):
        return True if len(self.subworkflows.keys()) else False

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


class CromwellException(Exception):
    """Generic exception when Cromwell does not return metadata."""

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

    def __init__(self, workflows_url, workflow_id, data=None):
        """
        Initialize Metadata object by retrieving metadata from Cromwell,
        including subworkflows, and initialize Task objects.
        """
        logger = logging.getLogger(__package__)
        self.workflows_url = workflows_url
        self.workflow_id = workflow_id
        self.tasks = None
        self.data = None
        if data:
            self.data = data
        else:
            logger.debug(f"Get metadata for {workflow_id}")
            self._get_data()
        self._init_tasks()

    def _get_data(self):
        """GET record from Cromwell REST server."""
        url = f"{self.workflows_url}/{self.workflow_id}/metadata"
        try:
            response = requests.get(url)
        except requests.ConnectionError as error:
            raise error
        response.raise_for_status()
        self.data = response.json()

    def init_with_json_str(self, json_str):
        try:
            self.data = json.loads(json_str)
        except Exception as error:
            raise ValueError(f"Not valid json: {error}")

    def _init_tasks(self):
        """Initialize and save Task objects."""
        logger = logging.getLogger(__package__)
        self.tasks = []  # Task objects
        self.jobs = {}  # job_id => dict
        self.subworkflows = {}  # workflow_id => metadata obj
        if "calls" in self.data:
            calls = self.data["calls"]
            for task_name in calls.keys():
                logger.debug(f"Workflow {self.workflow_id}: Init task {task_name}")
                task = Task(self.workflows_url, task_name, calls[task_name])
                self.tasks.append(task)
                for job_id in task.jobs:
                    job_info = task.jobs[job_id]
                    logger.debug(f"Workflow {self.workflow_id}: job {job_id} = {job_info}")
                    self.jobs[job_id] = job_info
                if task.is_subworkflow:
                    for sub_id, sub_meta in task.subworkflows.items():
                        self.subworkflows[sub_id] = sub_meta

    def get(self, param, default=None):
        """Get a section of the document."""
        return self.data.get(param, default)

    def is_subworkflow(self):
        return True if "parentWorkflowId" in self.data else False

    def get_job_info(self, job_id):
        """Get task name and attempt if found, otherwise None."""
        logger = logging.getLogger(__package__)
        logger.debug(f"Get info for job {job_id}")
        job_id = int(job_id)
        if job_id in self.jobs:
            return self.jobs[job_id]
        else:
            return None


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
        return Metadata(self.workflows_url, workflow_id, data)

    def get_all_metadata_json(self, workflow_id: str):
        """Get dict of all runs => metadata json for run and all subworkflows.

        :param workflow_id: primary key used by Cromwell
        :type workflow_id: str
        :return: all metadata docs { workflow_id => metadata json }
        :rtype: dict
        """
        result = {}
        metadata = Metadata(self.workflows_url, workflow_id)
        result[workflow_id] = metadata.data
        for sub_id, sub_meta in metadata.subworkflows.items():
            result[sub_id] = sub_meta.data

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

    def submit(self, wdl_file: str, json_file: str, zip_file: str = None) -> int:
        """
        Submit a run to Cromwell.
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
        try:
            response = requests.post(self.workflows_url, files=files)
        except Exception as error:
            logger.exception(f"Error submitting new run: {error}")
            raise error
        response.raise_for_status()
        run_id = response.json()["id"]
        return run_id

    def get_status(self, workflow_id: str):
        """
        Get the status of a workflow.
        """
        url = f"{self.workflows_url}/{workflow_id}/status"
        try:
            response = requests.get(url)
        except Exception as error:
            raise error
        response.raise_for_status()
        result = response.json()
        return result["status"]
