"""
Contains the WDLFile, WorkflowInputs, and Manifest file that constitute a user's workflow. A user
will simply specify the location of their WDL and provided inputs JSON file. A manifest file will record what
needs to be transferred over to Globus.
"""

import os
import shutil
import json
import subprocess
import re
import zipfile
import logging
import pathlib
import warnings
from jaws_client.config import Configuration
from jaws_client.wdl_runtime_validator import validate_wdl_runtime
from jaws_client.copy_progress import copy_with_progress_bar


config = Configuration()


def pretty_warning(message, category, filename, lineno, line=None):
    return f"WARNING: {message}\n"


warnings.formatwarning = pretty_warning


def join_path(*args):
    return os.path.join(*args)


def rsync(src, dest, options=["-rLtq"]):
    """Copy source to destination using rsync.

    :param src: Source path
    :type src: str
    :param dest: Destination path
    :type dest: str
    :param options: rsync options
    :type options: str
    :return: None
    """
    return subprocess.run(
        ["rsync", *options, src, dest], capture_output=True, text=True
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
    womtool_cmd = config.get("JAWS", "womtool").split()
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
    return True if isinstance(input, str) and re.match(".{0,2}/.+", input) else False


def recurse_list_type(inp_type):
    """With every recursion level crossed, the obj type is recursed as well.
    Eg. Array[Array[File]] -> Array[File] -> File.
    Eg. Array[Map[Int, String]] -> Map[Int, String]
    (will move to recurse_dict_type fuction to process this map further)"""
    return inp_type[inp_type.find("[") + 1: len(inp_type) - 1]


def recurse_dict_type(inp_type, dict_type):
    """
    Processes Map type and Pair type.
    Finds the positon of ',' that separates key, value (map) and left, right (pair)
    Returns corresponding types of key, values or left, right
    Eg. Map[Array[Array[File]], Pair[Int, String]] will return:
    key = Array[Array[File]] and value = Pair[Int, String]
    """
    start_idx = inp_type.find(dict_type) + len(dict_type)
    bracket_stack = []
    for idx in range(start_idx, len(inp_type)):
        if inp_type[idx] == "[":
            bracket_stack.append("[")
        if inp_type[idx] == "," and bracket_stack == []:
            break
        if inp_type[idx] == "]":
            bracket_stack.pop()
    key = inp_type[start_idx:idx]
    value = inp_type[idx + 2: len(inp_type) - 1]
    return key, value


def apply_filepath_op(obj, inp_type, operation):
    """
    Traverses the values of a dictionary and performs an operation on the value.

    :param obj: a python data structure (eg. list, dict, set, etc)
    :param inp_type: object type (eg. Map[Int, String], Array[Array[File]], etc)
    :param operation: an operation to perform on the values of a py data structure
    :return: output of the operation applied

    Note: Whenever operation is called a Boolean value is passed along with obj.
    This boolean values represents the flag is_file. It indicates if obj is a File.
    """
    if isinstance(obj, str):
        if "File" in inp_type:
            return operation(obj, True)
        else:
            return operation(obj, False)
    if isinstance(obj, int):
        return obj
    elif isinstance(obj, list):
        return [
            apply_filepath_op(i, recurse_list_type(inp_type), operation) for i in obj
        ]
    elif isinstance(obj, dict):
        if "Pair[" in inp_type:
            left_type, right_type = recurse_dict_type(inp_type, "Pair[")
            return {
                "Left": apply_filepath_op(obj["Left"], left_type, operation),
                "Right": apply_filepath_op(obj["Right"], right_type, operation),
            }
        if "Map[" in inp_type:
            key_type, value_type = recurse_dict_type(inp_type, "Map[")
            return {
                apply_filepath_op(k, key_type, operation): apply_filepath_op(
                    obj[k], value_type, operation
                )
                for k in obj
            }
    else:
        raise ValueError(
            f"cannot perform op={operation.__name__} to object of type {type(obj)}"
        )


class WdlFile:
    """
    A WDL object that can be queried for subworkflows, memory requirements and can also validate itself.

    The user submits a WDL file that describes their workflow. This class will validate the WDL file, and also
    keep track of any resource requirements (max_memory) required. It will also know its compressed file location.
    """

    def __init__(self, wdl_file_location, submission_id, contents=None):
        """
        Constructor for the WDL file.

        :param wdl_file_location:  location where the WDL file exists
        :param staging_subdir_path:  the directory where WDLs and cromwell files are all staged
        :param submission_id:  a uuid.uuid4() generated string
        """

        self.file_location = os.path.abspath(wdl_file_location)
        self.name = self._get_wdl_name(wdl_file_location)
        self.submission_id = submission_id
        self.contents = (
            contents if contents is not None else open(wdl_file_location, "r").read()
        )
        self._subworkflows = None
        self._max_ram_gb = None

    def _set_subworkflows(self, output):
        """
        Filters the output of WOMTool inputs -l so that it can parse through and grab the paths to the sub-workflows.

        It will collect the sub-workflows to a set of WdlFiles. The parent globus_basedir, staging_dir and path
        to the ZIP file are passed down to the sub workflows.

        Subworkflows are also checked for invalid 'backend' tag and a WdlError will be raised if found.

        :param output: stdout from WOMTool validation
        :return: set of WdlFile sub-workflows
        """
        # some stdout text include original input for
        # subprocess.run. This attempts to match and filter out that output
        # eg. ['/usr/local/anaconda3/bin/java', '-Xms512m', ...]
        command_line_regex = r"\[.*\]"
        out = output.splitlines()
        filtered_out_lines = ["Success!", "List of Workflow dependencies is:", "None"]
        subworkflows = set()
        for sub in out:
            match = re.search(command_line_regex, sub)
            if sub not in filtered_out_lines and not match:
                sub_wdl = WdlFile(sub, self.submission_id)
                sub_wdl.verify_wdl_has_no_backend_tags()
                subworkflows.add(sub_wdl)
        self._subworkflows = list(subworkflows)

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
            self.validate()
        return self._subworkflows

    @staticmethod
    def _check_missing_subworkflow_msg(stderr):
        missing = set()
        m = re.search("Failed to import workflow (.+).:", stderr)
        if m:
            for sub in m.groups():
                missing.add(sub)
            raise WdlError("Subworkflows not found: " + ", ".join(missing))

    def validate(self, compute_max_ram_gb):
        """
        Validates the WDL file using Cromwell's womtool and runtime validator.
        Any syntax errors from WDL will be raised in a WdlError.
        This is a separate method and not done automatically by the constructor because subworkflows() returns
        WdlFile objects and we wish to avoid running womtool multiple times unnecessarily.
        :param compute_max_ram_gb: maximum available ram for the site
        :ptype compute_max_ram_gb: float
        :return:
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Validating WDL, {self.file_location}")
        stdout, stderr = womtool("validate", "-l", self.file_location)
        self._set_subworkflows(stdout)
        if stderr:
            self._check_missing_subworkflow_msg(stderr)
            raise WdlError(stderr)
        self.verify_wdl_has_no_backend_tags()
        validate_wdl_runtime(self.contents, compute_max_ram_gb)

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
                    g = m.groups()
                    if g[0] == "memory":
                        mem = int(m.group(2))
                        prefix = m.group(3).lower()
                        mem = convert_to_gb(mem, prefix)
                        if mem > max_ram:
                            max_ram = mem
                    else:
                        raise WdlError(
                            "The 'mem' tag is deprecated; please use 'memory' instead"
                        )
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

        logger = logging.getLogger(__package__)
        logger.debug(f"Maximum RAM requested is {self._max_ram_gb}Gb")
        return self._max_ram_gb

    def copy_to(self, destination, permissions=0o664):
        shutil.copy(self.file_location, destination)
        os.chmod(destination, permissions)

    def verify_wdl_has_no_backend_tags(self):
        """
        Checks for disallowed "backend" keyword in a WDL file.

        Each Site has it's own configured backend; we do not allow users to specify which to
        use.  In particular, the "LOCAL" backend would cause jobs to run on the head node,
        which would cause problems if it consumed too much RAM.

        An exception is WdlError exception is raised if disallowed backend is found.

        :param wdl: the path to a WDL file
        :type wdl: str
        :return:
        """
        with open(self.file_location, "r") as fh:
            for line in fh:
                m = re.match(r"^\s*backend", line)
                if m:
                    raise WdlError("ERROR: WDLs may not contain 'backend' tag")

    def compress_wdls(self, staging_dir="."):
        """
        Create a new staging WDL and compress any subworkflow files into a ZIP file.

        The WDL is named based off of the submission ID and moved to a staging location that is specified in a
        configuration file. Then any subworkflow WDLs that are associated with that WDL are moved to a zip file and
        compressed.

        This method should be called on the main WDL object and shall include the subworkflows.
        If there are no subworkflows, the zipfile is not produced.

        :param staging_dir: path where files will be compressed to. Default is current directory.
        :type staging_dir: str
        :return: paths to the main wdl and the compressed file (latter may be None)
        """
        if not os.path.isdir(staging_dir):
            os.makedirs(staging_dir)

        # COPY MAIN WDL
        staged_wdl_filename = join_path(staging_dir, f"{self.submission_id}.wdl")
        self.copy_to(staged_wdl_filename)

        # IF NO SUBWORKFLOWS, DONE
        if not len(self.subworkflows):
            return staged_wdl_filename, None

        # ZIP SUBWORKFLOWS
        compressed_file_format = ".zip"
        compression_dir = pathlib.Path(os.path.join(staging_dir, self.submission_id))
        compression_dir.mkdir(parents=True, exist_ok=True)
        compressed_file = join_path(
            staging_dir, self.submission_id + compressed_file_format
        )
        main_wdl_dir = os.path.dirname(self.file_location)

        for subworkflow in self.subworkflows:
            sub_wdl_dir = os.path.dirname(subworkflow.file_location)
            sub_wdl_relative_path = os.path.relpath(sub_wdl_dir, start=main_wdl_dir)
            sub_filename = os.path.join(sub_wdl_relative_path, subworkflow.name)
            dirname = pathlib.Path(os.path.join(compression_dir, sub_wdl_relative_path))
            dirname.mkdir(mode=0o0770, parents=True, exist_ok=True)
            staged_sub_wdl = join_path(compression_dir, sub_filename)
            subworkflow.copy_to(staged_sub_wdl)
        try:
            os.remove(compressed_file)
        except FileNotFoundError:
            pass

        with zipfile.ZipFile(compressed_file, "w") as z:
            for sub_wdl in self.subworkflows:
                sub_wdl_dir = os.path.dirname(subworkflow.file_location)
                sub_wdl_relative_path = os.path.relpath(sub_wdl_dir, start=main_wdl_dir)
                sub_filename = os.path.join(sub_wdl_relative_path, sub_wdl.name)
                staged_sub_wdl = join_path(compression_dir, sub_filename)
                z.write(staged_sub_wdl, arcname=sub_filename)

        shutil.rmtree(compression_dir)
        return staged_wdl_filename, compressed_file

    def __eq__(self, other):
        return self.file_location == other.file_location

    def __hash__(self):
        return hash(self.file_location)

    def __repr__(self):
        return repr(self.file_location)


class WorkflowInputs:
    """
    Represents a JSON file where the input files of a WDL are specified.

    The user will submit the WDL, and JSON files to specify the workflow to JAWS. The inputs JSON
    will store the location of the input files. These input files are then pre-pended with the staging directory
    where JAWS Central will know where to pick them up.

    """

    def __init__(self, inputs_loc, submission_id, inputs_json=None, wdl_loc=None):
        """
        This class represents the inputs JSON which may contain both path and non-path elements.
        Inputs that look like relative paths and converted to absolute paths.
        Input paths are saved in a set, for ease of processing without traversing complex data structure.
        """
        self.submission_id = submission_id
        self.inputs_location = os.path.abspath(inputs_loc)
        self.basedir = os.path.dirname(self.inputs_location)
        inputs_json = (
            json.load(open(inputs_loc, "r")) if inputs_json is None else inputs_json
        )
        self.inputs_json = {}
        self.src_file_inputs = set()
        self.wdl_file_location = os.path.abspath(wdl_loc)
        self.wom_input_types = json.loads(self._get_inputs())
        for k in inputs_json:
            self.inputs_json[k] = apply_filepath_op(
                inputs_json[k], self.wom_input_types[k], self._gather_absolute_paths
            )

    def _get_inputs(self):
        """
        processes the associated wdl and returns the variable types for the workflow
        inputs through WOMTOOL's inputs command
        """
        stdout, stderr = womtool("inputs", self.wdl_file_location)
        if stderr:
            raise WdlError(stderr)
        return stdout

    def _gather_absolute_paths(self, element, is_file):
        """
        Helper method that recognizes file-like elements, converts relative to absolute paths,
        and adds them to the object's set of files.
        If element is a file and contains "http", do not get absolute paths
        """
        if is_file:
            if not element.startswith(
                ("http://", "ftp://", "https://")
            ) and not os.path.isabs(element):
                element = os.path.abspath(join_path(self.basedir, element))
            self.src_file_inputs.add(element)
        return element

    def move_input_files(self, destination, quiet=False):
        """
        Moves the input files defined in a JSON file to a destination.

        It will make any directories that are needed in the staging path and then copy the files.

        :param workflow_inputs: JSON file where inputs are specified.
        :param destination: path to where to moved the input files
        :return: list of the copied_files
        """
        copied_files = []

        for original_path in self.src_file_inputs:

            # if the path doesn't exist, it may refer to a path in a Docker container, so just
            # warn user and skip, without raising an exception
            if not os.path.exists(original_path):
                warnings.warn(f"Input path not found or inaccessible: {original_path}")
                continue

            staged_path = pathlib.Path(f"{destination}{original_path}")
            dest_path = staged_path.as_posix()
            copied_files.append(dest_path)

            if os.path.isdir(original_path):
                raise ValueError(
                    f"Invalid {original_path}; directories not supported for File types; use a tarball/zip instead."
                )

            dirname = pathlib.Path(os.path.dirname(staged_path))
            dirname.mkdir(mode=0o0770, parents=True, exist_ok=True)

            # Files must be copied in to ensure they are readable by the jaws and jtm users.  The group
            # will be set correctly as a result of the gid sticky bit and acl rules on the inputs dir.
            if quiet:
                shutil.copyfile(original_path, dest_path)
            else:
                copy_with_progress_bar(original_path, dest_path)
            os.chmod(dest_path, 0o0660)

        return copied_files

    def prepend_paths_to_json(self, staging_dir):
        """
        Modifies the paths of real files by prepending the base dir of the compute site.
        Only paths which exist are modified; other path-like strings presumably refer to files in the container.
        """
        destination_json = {}
        for k in self.inputs_json:
            destination_json[k] = apply_filepath_op(
                self.inputs_json[k],
                self.wom_input_types[k],
                self._prepend_path(staging_dir),
            )

        return WorkflowInputs(
            staging_dir,
            self.submission_id,
            destination_json,
            wdl_loc=self.wdl_file_location,
        )

    @staticmethod
    def _prepend_path(path_to_prepend):
        """
        The input files shall be transferred to the compute-site, so the paths in the inputs json file must be
        updated to reflect the new location by prepending the site's uploads basedir.
        Only paths which exist are modified because these are the only files which are transferred; the other
        paths presumably refer to files in a Docker container.
        If element is a file and contains "http", do not prepend given path to it
        """

        def func(path, is_file):
            if is_file:
                if not path.startswith(
                    ("http://", "ftp://", "https://")
                ) and os.path.exists(path):
                    return f"{path_to_prepend}{path}"
            return path

        return func

    def write_to(self, json_location):
        """
        Writes the JSON file to the specified location.

        :param json_location: location to write the new JSON file
        :return:
        """
        with open(json_location, "w") as json_inputs:
            json.dump(self.inputs_json, json_inputs, indent=4)


class Manifest:
    """
    A Manifest file includes all the files and their inode types that are transferred over via Globus.

    It is a TSV (tab-separate value) file.
    """

    def __init__(self, staging_dir, dest_dir):
        self.staging_dir = staging_dir
        self.dest_dir = dest_dir
        self.manifest = []

    def add(self, *args):
        """
        Adds the specified filepath to the manifest TSV. It will include the source path, the destination path
        and its inode type (file or directory).

        :param filepath: path of the file to add to TSV
        :return:
        """
        for filepath in args:
            if filepath is None:
                continue  # zip file may be none
            inode_type = "D" if os.path.isdir(filepath) else "F"
            src_rel_path = os.path.relpath(filepath, self.staging_dir)
            dest_rel_path = f"{self.dest_dir}/{src_rel_path}"
            self.manifest.append([filepath, dest_rel_path, inode_type])

    def write_to(self, write_location):
        """
        Writes the TSV file to a specified location.

        :param write_location: a file path
        :return:
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Writing manifest file: {write_location}")
        with open(write_location, "w") as f:
            for src, dest, inode_type in self.manifest:
                f.write(f"{src}\t{dest}\t{inode_type}\n")


class WorkflowError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WdlError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkflowInputsError(Exception):
    def __init__(self, message):
        super().__init__(message)
