"""
Contains the WDLFile, WorkflowInputs, and Manifest file that constitute a user's workflow. A user
will simply specify the location of their WDL and provided inputs JSON file. A manifest file will record what
needs to be transferred over to Globus.
"""

import os
import shutil
import json
import stat
import subprocess
import re
import zipfile
import logging
import pathlib
from pathlib import Path

from jaws_client import config


def join_path(*args):
    return os.path.join(*args)


def rsync(src, dest):
    """Copy source to destination using rsync.

    :param src: Source path
    :type src: str
    :param dest: Destination path
    :type dest: str
    :return: None
    """
    cmd = f"rsync -rLptq {src} {dest}"
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode:
        raise IOError(
            f"Failed to rsync {src} to {dest}: " + stderr.decode("utf-8").strip()
        )


def convert_to_gb(mem, prefix):
    """
    Convert units to gigabytes

    :param mem: amount of memory
    :type mem: int
    :param prefix: units prefix
    :type prefix str
    :return: memory in gigabytes
    """
    conversion_table = {
        "g": lambda x: x,
        "k": lambda x: x / 1048576,
        "m": lambda x: x / 1024,
        "t": lambda x: x * 1024,
    }
    convert = conversion_table.get(prefix, lambda x: x / 1073741824)
    gb_mem = convert(mem)
    gb_mem = int(gb_mem + 0.99)
    return gb_mem


def womtool(*args):
    """
    Call WOMTool

    WOMTool is a Broad Institute tool used for validating WDLs. More information can be found here:
    https://cromwell.readthedocs.io/en/stable/WOMtool/

    :param args: WOMTool arguments
    :return: stdout and stderr of WOMTool completed process
    """
    womtool_cmd = config.conf.get("JAWS", "womtool").split()
    womtool_cmd.extend(list(args))
    proc = subprocess.run(womtool_cmd, capture_output=True, text=True)
    return proc.stdout, proc.stderr


def value_generator(val):
    """
    Helper function that recurses through the values of a dictionary

    :param val:  dictionary value
    :return:
    """
    if isinstance(val, str):
        yield val
    elif isinstance(val, dict):
        for k, v in val.items():
            yield k
            yield from value_generator(v)
    elif isinstance(val, list):
        for k in val:
            yield from value_generator(k)


def values(inputs_map):
    """
    Generates the appropriate keys depending on the object we
    are iterating over.

    The keys are going to be the name of the resource, and the values
    are going to be that actual resource. For example,

    {"jgi_sample.fasta": "/path/to/file"}

    the key is the wdl resource name, and the value is what we want to iterate over.
    This can be a string, list, dictionary, or list of lists, etc.


    :param inputs_map: a JSON dictionary
    :return: yields the values of a dictionary
    """
    for val in inputs_map.values():
        yield from value_generator(val)


def is_refdata(filepath):
    return filepath.startswith("/refdata")


def looks_like_file_path(input):
    return isinstance(input, str) and "/" in input


def is_file_accessible(filename):
    """
    Checks if the file can be access by the jaws shared user account.

    Validates files by first checking if the file is the special case `/refdata` folder.
    Returns true if that is the case since /refdata is a database that exists outside of the local
    filesystem.

    Then it checks whether the group file permissions are properly set on the group file permissions.
    Assumes that there exists a common group permission that the shared user belongs in. It will check
    whether the datasets are readable and executable (executable is needed for moving files).

    :param filename: str name of the file
    :return: True iff group perms set.
    """
    if is_refdata(filename):  # Ignores refdata since it is a special case directory
        return True
    else:
        st = os.stat(filename)
        group_permission = Path(filename).group()

        # Check if file exists or accessible
        if not os.path.exists(filename):
            return False
        # Check if group permission is correct
        if group_permission != config.conf.get("JAWS", "shared_endpoint_group"):
            return False
        # Check if file has group readable and executable
        if not bool(st.st_mode & stat.S_IRGRP) or not bool(st.st_mode & stat.S_IXGRP):
            return False
        return True


def apply_filepath_op(obj, operation):
    """
    Traverses the values of a dictionary and performs an operation on the value.

    :param obj: a python data structure (eg. list, dict, set, etc)
    :param operation: an operation to perform on the values of a python data structure
    :return: output of the operation applied
    """
    if isinstance(obj, str):
        return operation(obj)
    if isinstance(obj, int):
        return obj
    elif isinstance(obj, list):
        return [apply_filepath_op(i, operation) for i in obj]
    elif isinstance(obj, dict):
        return {
            apply_filepath_op(k, operation): apply_filepath_op(obj[k], operation)
            for k in obj
        }
    else:
        raise ValueError(
            f"cannot perform op={operation.__name__} to object of type {type(obj)}"
        )


def compress_wdls(main_wdl, output_basename, output_dir="."):
    """
    Create a new staging WDL and compress any subworkflow files into a ZIP file.

    The WDL is named based off of the submission ID and moved to a staging location that is specified in a
    configuration file. Then any subworkflow WDLs that are associated with that WDL are moved to a zip file and
    compressed.
    If there are no subworkflows, the zipfile is not produced..

    :param main_wdl: a WdlFile object
    :param output_basename: string used to generate wdl and zip output file names
    :param output_dir: path where files will be compressed to. Default is current directory.
    :return: paths to the main wdl and the compressed file (latter may be None)
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # WRITE SANITIZED MAIN WDL
    modified_wdl = main_wdl.sanitized_wdl()
    staged_wdl_filename = join_path(output_dir, f"{output_basename}.wdl")
    modified_wdl.write_to(staged_wdl_filename)

    # IF NO SUBWORKFLOWS, DONE
    if not len(main_wdl.subworkflows):
        return staged_wdl_filename, None

    # ZIP SUBWORKFLOWS
    compressed_file_format = ".zip"
    compression_dir = pathlib.Path(os.path.join(output_dir, output_basename))
    compression_dir.mkdir(parents=True, exist_ok=True)
    compressed_file = join_path(output_dir, output_basename + compressed_file_format)

    # write temporary sanitized/modified subworkflows
    for subworkflow in main_wdl.subworkflows:
        modified_sub = subworkflow.sanitized_wdl()
        new_subwf_path = join_path(compression_dir, subworkflow.name)
        modified_sub.write_to(new_subwf_path)

    # if the output zip file already exists, delete it before writing a new one
    try:
        os.remove(compressed_file)
    except FileNotFoundError:
        pass

    # compress modified subworkflow files
    if compressed_file_format == ".zip":
        with zipfile.ZipFile(compressed_file, "w") as z:
            for sub_wdl in main_wdl.subworkflows:
                staged_sub_wdl = join_path(compression_dir, sub_wdl.name)
                z.write(staged_sub_wdl, arcname=sub_wdl.name)
    else:
        raise ValueError("Only zip files are supported at this time")

    # delete the temporary modified subworkflow files
    shutil.rmtree(compression_dir)
    return staged_wdl_filename, compressed_file


class WdlFile:
    """
    A WDL object that can be queried for subworkflows, memory requirements and can also validate itself.

    The user submits a WDL file that describes their workflow. This class will validate the WDL file, and also
    keep track of any resource requirements (max_memory) required. It will also know its compressed file location.
    """

    def __init__(self, wdl_file_location, contents=None):
        """
        Constructor for the WDL file.

        :param wdl_file_location:  location where the WDL file exists
        :param staging_subdir_path:  the directory where WDLs and cromwell files are all staged
        """

        self.logger = logging.getLogger(__package__)
        self.file_location = os.path.abspath(wdl_file_location)
        self.name = self._get_wdl_name(wdl_file_location)
        self.contents = (
            contents if contents is not None else open(wdl_file_location, "r").read()
        )

        self._subworkflows = None
        self._max_ram_gb = None

    def _filter_subworkflows(self, output):
        """
        Filters the output of WOMTool inputs -l so that it can parse through and grab the paths to the sub-workflows.

        It will collect the sub-workflows to a set of WdlFiles. The parent globus_host_path, staging_dir and path
        to the ZIP file are passed down to the sub workflows.

        :param output: stdout from WOMTool validation
        :return: set of WdlFile sub-workflows
        """
        out = output.splitlines()
        filtered_out_lines = ["Success!", "List of Workflow dependencies is:", "None"]
        subworkflows = set()
        for sub in out:
            if sub.startswith("http://") or sub.startswith("https://"):
                # subworkflow file not provided, Cromwell will GET automatically
                continue
            if sub not in filtered_out_lines:
                subworkflows.add(WdlFile(sub))
        return subworkflows

    @property
    def subworkflows(self):
        """
        Property that lazily evaluates the subworkflows of a WDL file.

        Since we are calling an external subprocess that adds a lot of overhead, the self._subworkflows instance
        variable is only set when this property method is called to avoid adding that overhead when the WDL file
        object is constructed.

        :return: set of WdlFiles
        """
        if self._subworkflows is None:
            stdout, stderr = womtool("validate", "-l", self.file_location)
            self._subworkflows = self._filter_subworkflows(stdout)
        return self._subworkflows

    @staticmethod
    def _check_missing_subworkflow_msg(stderr):
        missing = set()
        m = re.search("Failed to import workflow (.+).:", stderr)
        if m:
            for sub in m.groups():
                missing.add(sub)
            raise WdlError("Subworkflows not found: " + ", ".join(missing))

    def validate(self):
        """
        Validates using WOMTool the WDL file. Any syntax errors from WDL will be raised in a WdlError
        :return:
        """
        self.logger.info(f"Validating WDL, {self.file_location}")
        _, stderr = womtool("validate", "-l", self.file_location)
        if stderr:
            self._check_missing_subworkflow_msg(stderr)
            raise WdlError(stderr)

    @staticmethod
    def _get_wdl_name(file_location):
        return os.path.basename(file_location)

    def _calculate_max_ram_gb(self):
        """
        Helper method for calculating the maximum memory in a WDL.

        :return: maximum memory specified in a WDL.
        """
        max_ram = 0
        with open(self.file_location, "r") as f:
            for line in f:
                m = re.match(
                    r"^\s+(mem|memory)\s*[:=]\s*\"?(\d+\.?\d*)([kKmMgGtT])\"?", line
                )
                if m:
                    mem = int(m.group(2))
                    prefix = m.group(3).lower()
                    mem = convert_to_gb(mem, prefix)
                    if mem > max_ram:
                        max_ram = mem
        return max_ram

    @property
    def max_ram_gb(self):
        """
        Lazily evaluates the maximum memory specified in a WDL and its imported WDL files.

        The evaluation is lazy because it also has to call WOMTool to grab the subworkflows and thus
        can be an expensive operation.

        :return: maximum RAM in GB, rounded up to the nearest GB.
        :rtype: int
        """
        if self._max_ram_gb is None:
            self._max_ram_gb = self._calculate_max_ram_gb()
            for subworkflow_file in self.subworkflows:
                if subworkflow_file.max_ram_gb > self._max_ram_gb:
                    self._max_ram_gb = subworkflow_file.max_ram_gb

        self.logger.info(f"Maximum RAM requested is {self._max_ram_gb}Gb")
        return self._max_ram_gb

    def sanitized_wdl(self):
        contents = self._remove_invalid_backends()
        return WdlFile(self.file_location, contents=contents)

    def write_to(self, destination):
        with open(destination, "w") as new_wdl:
            new_wdl.write(self.contents)

    def _remove_invalid_backends(self):
        """
        Removes the backend keyword in a WDL file.

        This sanitizes the WDL so that any declared backends are taken off. The reason for this is because we
        are ONLY using JGI Task Manager (JTM) as the backend and want to prevent any WDL from using a different
        backend.

        :param wdl: the path to a WDL file
        :return: the contents of the file without backend keyword
        """
        contents = ""
        with open(self.file_location, "r") as fo:
            for line in fo:
                m = re.match(r"^\s*backend", line)
                if not m:
                    contents += line
        return contents

    def __eq__(self, other):
        return self.file_location == other.file_location

    def __hash__(self):
        return hash(self.file_location)

    def __repr__(self):
        return repr(self.file_location)


def copy_input_files(workflow_inputs, globus_host_path, destination):
    """
    Moves the input files defined in a JSON file to a destination.

    It will make any directories that are needed in the staging path, and then either symlink them
    if they are a globus path or use rsync to move the files.

    :param workflow_inputs: JSON file where inputs are specified.
    :type workflow_inputs: str
    :param globus_host_path: The root dir of the Globus endpoint
    :type globus_host_path: str
    :type globus_host_path: str
    :param destination: path to where to moved the input files
    :type destination: str
    :return: list of the moved_files
    """
    staged_files = []

    for original_path in workflow_inputs.src_file_inputs:
        staged_path = pathlib.Path(f"{destination}{original_path}")
        staged_files.append(staged_path.as_posix())
        if os.path.isdir(original_path):
            staged_path.mkdir(mode=0o0770, parents=True, exist_ok=True)
        else:
            dirname = pathlib.Path(os.path.dirname(staged_path))
            dirname.mkdir(mode=0o0770, parents=True, exist_ok=True)

        # globus paths are accessible via symlink
        if original_path.startswith(f"{globus_host_path}"):
            if not os.path.exists(staged_path):
                os.symlink(original_path, staged_path.as_posix())
        else:
            rsync(original_path, staged_path.as_posix())
    return staged_files


class WorkflowInputs:
    """
    Represents a JSON file where the input files of a WDL are specified.

    The user will submit the WDL, and JSON files to specify the workflow to JAWS. The inputs JSON
    will store the location of the input files. These input files are then pre-pended with the staging directory
    where JAWS Central will know where to pick them up.

    """

    def __init__(self, inputs_loc, submission_id, inputs_json=None):

        self.submission_id = submission_id
        self.inputs_location = os.path.abspath(inputs_loc)
        self.basedir = os.path.dirname(self.inputs_location)

        # JSON inputs could contain relative paths, so we process the JSON file to include absolute paths
        inputs_json = (
            json.load(open(inputs_loc, "r")) if inputs_json is None else inputs_json
        )
        self.inputs_json = {}
        for k in inputs_json:
            self.inputs_json[k] = apply_filepath_op(
                inputs_json[k], self._relative_to_absolute_paths
            )

        self._src_file_inputs = None

    def prepend_paths_to_json(self, staging_dir):
        """
        Modifies the JSON file to include the adjusted paths to the JAWS Central staging area.

        :return: contents of the modified JSON file which is a WorkflowInputs object
        """
        destination_json = {}
        for k in self.inputs_json:
            destination_json[k] = apply_filepath_op(
                self.inputs_json[k], self._prepend_path(staging_dir)
            )
        return WorkflowInputs(
            self.inputs_location, self.submission_id, inputs_json=destination_json
        )

    @property
    def src_file_inputs(self):
        """
        Gathers all the input files specified in the JSON file into a set of files.

        This allows for easy traversal of all the input files rather than have to recurse down the dictionary
        and attempt to determine whether an element is a file or a keyword.

        :return: set of file paths
        """
        if self._src_file_inputs is None:
            self._src_file_inputs = self._gather_paths()
        return self._src_file_inputs

    def validate(self):
        """
        Validates all of the input files in the JSON.

        A valid file is considered one that exists and that can be read and executed by the shared Globus
        user account. If a file has the wrong permissions, or is not at the specified location,
        a WorkflowError is raised to alert the user.
        :return:
        """
        not_accessible = []
        for filepath in values(self.inputs_json):
            if looks_like_file_path(filepath):
                if not is_file_accessible(filepath):
                    not_accessible.append(filepath)
        if not_accessible:
            group_owner = config.conf.get("JAWS", "shared_endpoint_group")
            msg = "File(s) not accessible:\n" + "\n".join(not_accessible)
            msg += f"\nPlease make sure that the group owner is set to {group_owner} and the group permissions are "
            msg += "set to readable and executable."
            raise SystemExit(msg)

    def write_to(self, json_location):
        """
        Writes the modified JSON file to the specified location.

        A modified file is one that includes the JAWS central staging path prepended to all the input file paths.

        :param json_location: location to write the new JSON file
        :return:
        """
        with open(json_location, "w") as json_inputs:
            json.dump(self.inputs_json, json_inputs, indent=4)

    def _relative_to_absolute_paths(self, element):
        """
        Helper method that converts any relative paths in a JSON into absolute paths.
        It is applied to each value.
        """
        if looks_like_file_path(element) and not is_refdata(element):
            new_path = element
            if not os.path.isabs(element):
                new_path = join_path(self.basedir, element)
                new_path = os.path.abspath(new_path)
                return new_path
        return element

    @staticmethod
    def _prepend_path(path_to_prepend):
        def func(path):
            if looks_like_file_path(path) and not is_refdata(path):
                return f"{path_to_prepend}{path}"
            return path

        return func

    def _gather_paths(self):
        """
        Helper method that aggregates all the paths in a JSON file
        """
        paths = set()
        for element in values(self.inputs_json):
            if looks_like_file_path(element) and not is_refdata(element):
                paths.add(element)
        return paths


class Manifest:
    """
    A Manifest file includes all the files and their inode types that are transferred over via Globus.

    File paths are as expected by Globus transfer service; they appear as an absolute path but Globus shall
    prepend the endpoint's host_path automatically.

    e.g. A transfer task of "/mydir/myfile" will actually transfer the file "/tmp/gsharing/mydir/myfile" if
         the host_path of the endpoint is "/tmp/gsharing".

    It is a TSV (tab-separate value) file.
    """

    def __init__(self, src_host_path, src_dir, dest_host_path, dest_dir):
        """
        :param src_host_path: root directory of the source globus endpoint
        :type src_host_path: str
        :param dest_host_path: root directory of the destination globus endpoint
        :type dest_host_path: str
        """
        self.logger = logging.getLogger(__package__)
        self.src_host_path = src_host_path
        self.src_dir = src_dir
        self.dest_host_path = dest_host_path
        self.dest_dir = dest_dir
        self.manifest = []

    def add(self, *args):
        """
        Adds the specified filepath to the manifest TSV. It will include the source path, the destination path
        and its inode type (file or directory).

        :param filepath: full path of the file to add to TSV
        :return:
        """
        for filepath in args:
            if filepath is None:
                continue  # zip file may be none
            inode_type = "D" if os.path.isdir(filepath) else "F"

            # virtual paths, using host_path as root
            src_path = os.path.join("/", os.path.relpath(filepath, self.src_host_path))
            dest_path = os.path.join(
                "/",
                os.path.relpath(self.dest_dir, self.dest_host_path),
                os.path.relpath(src_path, self.src_dir),
            )

            # save
            self.manifest.append([src_path, dest_path, inode_type])

    def write_to(self, write_location):
        """
        Writes the TSV file to a specified location.

        :param write_location: a file path
        :return:
        """
        self.logger.info(f"Writing file manifest to {write_location}")
        with open(write_location, "w") as f:
            for src, dest, inode_type in self.manifest:
                f.write(f"{src}\t{dest}\t{inode_type}\n")


class WorkflowError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WdlError(Exception):
    def __init__(self, message):
        super().__init__(message)
