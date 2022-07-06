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
import uuid
from jaws_client.config import Configuration
from jaws_client.wdl_runtime_validator import validate_wdl_runtime
from jaws_client.copy_progress import copy_with_progress_bar


config = Configuration()

REFDATA_DIR = '/refdata'  # this folder is mounted to every task's container


def pretty_warning(message, category, filename, lineno, line=None):
    return f"WARNING: {message}\n"


warnings.formatwarning = pretty_warning


def join_path(*args):
    return os.path.join(*args)


def mkdir(path, perms=0o0775):
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    if perms:
        os.chmod(path, perms)


def convert_to_gb(mem, prefix):
    """
    Convert units to gigabytes

    :param mem: amount of memory
    :type mem: int
    :param prefix: units prefix
    :type prefix str
    :return: memory in gigabytes
    :rtype: int
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
    return filepath.startswith(REFDATA_DIR)


def looks_like_file_path(input):
    return True if isinstance(input, str) and re.match(".{0,2}/.+", input) else False


def recurse_list_type(inp_type):
    """With every recursion level crossed, the obj type is recursed as well.
    Eg. Array[Array[File]] -> Array[File] -> File.
    Eg. Array[Map[Int, String]] -> Map[Int, String]
    (will move to recurse_dict_type fuction to process this map further)"""
    return inp_type[inp_type.find("[") + 1 : len(inp_type) - 1]  # noqa


def recurse_dict_type(inp_type, dict_type):
    """
    Processes Map type and Pair type.
    Finds the positon of ',' that separates key, value (map) and left, right (pair)
    Returns corresponding types of key, values or left, right
    Eg. Map[Array[Array[File]], Pair[Int, String]] will return:
    key = Array[Array[File]] and value = Pair[Int, String]
    """
    paren_map = {"{": "}", "[": "]"}
    start_idx = inp_type.find(dict_type) + len(dict_type)
    bracket_stack = []
    for idx in range(start_idx, len(inp_type)):
        element = inp_type[idx]
        if element == "[" or element == "{":
            bracket_stack.append(element)
        if element == "," and bracket_stack == []:
            break
        if (element == "]" or element == "}") and paren_map[element] == bracket_stack[
            -1
        ]:
            bracket_stack.pop()
    key = inp_type[start_idx:idx]
    value = inp_type[idx + 1 : len(inp_type) - 1]  # noqa
    return key, value


def process_wom_composite(inp_type):
    """
    Takes WomCompositeType and parses the string into a dict
    Eg. WomCompositeType{\nname->String\nbatch->Pair[Int,Int]
    \nlocations-> Array[File]\ninfo->Map[String,String]\n}
    will be parsed into
    {
        'name': 'String',
        'batch': 'Pair[Int,Int]',
        'locations': 'Array[File]',
        'info': 'Map[String,String]'
    }
    """
    dict_type = "WomCompositeType{"
    start_idx = inp_type.find(dict_type) + len(dict_type)
    elements = inp_type[start_idx : len(inp_type) - 1]  # noqa
    elements = elements.strip().split()
    struct_obj = {}
    for element in elements:
        pair = element.split("->")
        struct_obj[pair[0]] = pair[1]
    return struct_obj


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
    inp_type = inp_type.replace(" ", "")
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
        if inp_type.startswith("WomCompositeType"):
            inp_type = process_wom_composite(inp_type)
            struct_obj = {}
            for variable in inp_type:
                struct_obj[variable] = apply_filepath_op(
                    obj[variable], inp_type[variable], operation
                )
            return struct_obj
        if inp_type.startswith("Pair["):
            left_type, right_type = recurse_dict_type(inp_type, "Pair[")
            return {
                "Left": apply_filepath_op(obj["Left"], left_type, operation),
                "Right": apply_filepath_op(obj["Right"], right_type, operation),
            }
        if inp_type.startswith("Map["):
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


class Run:
    """
    A Run is comprised of a Workflow (WDL file(s)) and input files/parameters (inputs JSON file).
    """

    def __init__(self, wdl_file, json_file, input_dir, input_site_id, output_basedir, **kwargs):
        """
        A Run requires the main WDL and inputs JSON.  If the main WDL uses subworkflows,
        a ZIP archive may optionally be supplied, otherwise the object will attempt to
        create it (requires `womtool` to be installed).
        :param wdl_file: Path to main WDL file
        :ptype wdl_file: str
        :param json_file: Path to inputs JSON file
        :ptype json_file: str
        :param input_dir: Base dir where Run WDL and infiles shall be copied.
        :ptype input_dir: str
        :param input_site_id: Name of the input/submission jaws-site
        :ptype input_site_id: str
        :param output_basedir: Base dir where Run results are located.
        :ptype output_basedir: str
        """
        self.input_dir = input_dir
        self.input_site_id = input_site_id.upper()
        self.output_basedir = output_basedir

        subworkflows_zip_file = None
        if "subworkfows_zip_file" in kwargs:
            subworkflows_zip_file = kwargs["subworkflows_zip_file"]
        quiet = False
        if "quiet" in kwargs:
            quiet = True if kwargs["quiet"] else False

        # validate WDL file(s)
        self.wdl = WdlFile(wdl_file)
        self.wdl.validate()

        # validate inputs json
        self.inputs = WorkflowInputs(json_file, wdl_file)

        # The submission_id is a unique identifier used before the run_id is assigned.
        # It is used to name the temporary output directory and to name the input
        # workflow (WDL, JSON, ZIP) files.
        self.submission_id = str(uuid.uuid4())

        # This folder is where the outputs shall be returned.  Make the dir now --
        # the ACL rules will ensure the folder is writeable by `jaws`
        self.output_dir = f"{output_basedir}/{self.submission_id}"
        try:
            mkdir(self.output_dir)
        except Exception as error:
            raise IOError(f"Unable to create output dir: {error}")

        # Create optional subworkflows Zip-file and copy Run all inputs to the
        # input-dir, a folder where JAWS user can read them.  The subworkflows ZIP may be
        # supplied by the user or it will be created automatically.
        staged_wdl_file, zip_file = self.wdl.prepare_wdls(
            self.input_dir, self.submission_id, subworkflows_zip_file
        )

        # Run infiles (identified in the inputs JSON file) are copied to the input-dir,
        # so they may be read by JAWS user
        site_input_dir = os.path.join(self.input_dir, self.input_site_id)
        try:
            copied_files = self.inputs.copy_input_files(site_input_dir, quiet)
        except Exception as error:
            raise IOError(f"Unable to copy input files: {error}")

        # generate a new inputs JSON object with the paths pointing to the files copied under
        # the inputs-dir
        staged_inputs_json_file = f"{self.input_dir}/{self.submission_id}.json"
        self.inputs.prepend_paths_to_json(self.input_site_id)
        self.inputs.write_to(staged_inputs_json_file)

        # generate a (unique) list of all files to be transferred to the compute-site
        self.manifest = Manifest(input_dir)
        self.manifest.add(
            staged_inputs_json_file, staged_wdl_file, zip_file, *copied_files
        )

        # Copy workflow inputs (WDL, JSON, ZIP) to the output dir for the user's reference.
        # These copies are not used by JAWS.
        shutil.copy(json_file, self.output_dir)
        shutil.copy(wdl_file, self.output_dir)
        if subworkflows_zip_file:
            # if user-supplied zip, copy it to keep original filename
            shutil.copy(subworkflows_zip_file, self.output_dir)
        elif zip_file:
            # else copy automatically generated zip (if exists)
            shutil.copy(zip_file, f"{self.output_dir}/subworkflows.zip")


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
        :param contents: contents of file, if has already been read
        """

        self.file_location = os.path.abspath(wdl_file_location)
        self.name = self._get_wdl_name(wdl_file_location)
        self._subworkflows = None
        self.max_ram_gb = 0

        self.contents = (
            contents if contents is not None else open(wdl_file_location, "r").read()
        )
        self.lines = self.contents.split("\n")

        # the max RAM required by any task must be supplied to JAWS with the Run
        self.set_max_ram_gb()

        # specifying the Cromwell backend to use (e.g. LOCAL) is disallowed
        self.verify_wdl_has_no_backend_tags()

        # verify required tags are set (e.g. memory)
        validate_wdl_runtime(self.contents)

    def _set_subworkflows(self, output):
        """
        Filters the output of WOMTool inputs -l so that it can parse through and grab the paths to the sub-workflows.

        It will collect the sub-workflows to a set of WdlFiles. The parent basedir, staging_dir and path
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
            if (
                sub not in filtered_out_lines
                and not match
                and not sub.startswith("http")
            ):
                sub_wdl = WdlFile(sub)
                subworkflows.add(sub_wdl)
                if sub_wdl.max_ram_gb > self.max_ram_gb:
                    self.max_ram_gb = sub_wdl.max_ram_gb
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
    def missing_subworkflows_error_msg(stderr):
        missing = set()
        m = re.search("Failed to import workflow (.+).:", stderr)
        if m:
            for sub in m.groups():
                missing.add(sub)
            raise WdlError("Subworkflows not found: " + ", ".join(missing))

    def validate(self):
        """
        Validates the WDL file using Cromwell's womtool.
        Any syntax errors from WDL will be raised in a WdlError.
        Any subworkflows will also be identified and initialized.
        This is a separate method and not done automatically by the constructor because subworkflows() returns
        WdlFile objects and we wish to avoid running womtool multiple times unnecessarily.
        :return:
        """
        logger = logging.getLogger(__package__)
        logger.debug(f"Validating WDL, {self.file_location}")

        # validate WDL using womtool
        stdout, stderr = womtool("validate", "-l", self.file_location)

        # define subworkflow WDLs, if any
        self._set_subworkflows(stdout)
        if stderr:
            self.missing_subworkflows_error_msg(stderr)
            raise WdlError(stderr)

    @staticmethod
    def _get_wdl_name(file_location):
        return os.path.basename(file_location)

    def set_max_ram_gb(self):
        """
        Helper method for calculating the maximum memory in a WDL.

        :return: maximum memory specified in a WDL.
        """
        for line in self.lines:
            m = re.match(r"^\s+memory\s*[:=]\s*\"?(\d+\.?\d*)([kKmMgGtT])\"?", line)
            if m:
                g = m.groups()
                mem = int(g[0])
                prefix = g[1].lower()
                mem = convert_to_gb(mem, prefix)
                if mem > self.max_ram_gb:
                    self.max_ram_gb = mem

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
        for line in self.lines:
            m = re.match(r"^\s*backend", line)
            if m:
                raise WdlError("ERROR: WDLs may not contain 'backend' tag")

    def prepare_wdls(self, staging_dir, submission_id, user_zip=None):
        """
        Create a new staging WDL and compress any subworkflow files into a ZIP file.

        The WDL is named based off of the submission ID and moved to a staging location that is specified in a
        configuration file. Then any subworkflow WDLs that are associated with that WDL are moved to a zip file and
        compressed.

        This method should be called on the main WDL object and shall include the subworkflows.
        If there are no subworkflows, the zipfile is not produced.

        The user may supply a prepared zip file, which is copied instead of producing one.

        :param staging_dir: path where files will be compressed to. Default is current directory.
        :type staging_dir: str
        :param user_zip: path to user-supplied subworkflows zip file (optional)
        :type user_zip: str
        :return: paths to the main wdl and the compressed file (latter may be None)
        """
        if not os.path.isdir(staging_dir):
            os.makedirs(staging_dir)

        # the main wdl must be copied to submission path and named by submission id
        staged_wdl_filename = join_path(staging_dir, f"{submission_id}.wdl")
        self.copy_to(staged_wdl_filename)

        # any subworkflows must be zipped and also copied to submission path and are
        # also named using the submission id
        compressed_file = join_path(staging_dir, submission_id + ".zip")
        if not len(self.subworkflows):
            return staged_wdl_filename, None
        elif user_zip:
            shutil.copy(user_zip, compressed_file)
            os.chmod(compressed_file, 0o660)
        else:
            self.compress_subworkflows(compressed_file, staging_dir, submission_id)
        return staged_wdl_filename, compressed_file

    def compress_subworkflows(
        self, compressed_file: str, staging_dir: str, submission_id: str
    ):
        """
        Generate ZIP file of required subworkflows.
        """
        compression_dir = pathlib.Path(os.path.join(staging_dir, submission_id))
        compression_dir.mkdir(parents=True, exist_ok=True)
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

    def __init__(self, inputs_loc, wdl_loc, inputs_json=None):
        """
        This class represents the inputs JSON which may contain both path and non-path elements.
        Inputs that look like relative paths are converted to absolute paths.
        Input paths are saved in a set, for ease of processing without traversing complex data structure.
        """
        self.inputs_location = os.path.abspath(inputs_loc)
        self.basedir = os.path.dirname(self.inputs_location)
        if not inputs_json:
            try:
                inputs_json = json.load(open(inputs_loc, "r"))
            except json.JSONDecodeError as error:
                raise WorkflowInputsError(f"Your inputs JSON is invalid: {error}")
            except Exception as error:
                raise (f"Unable to read JSON infile: {error}")
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

    def copy_input_files(self, destination, quiet=False):
        """
        Copies the input files defined in a JSON file to a destination.

        It will make any directories that are needed in the staging path and then copy the files.

        :param workflow_inputs: JSON file where inputs are specified.
        :param destination: path to where to moved the input files
        :return: list of the copied_files
        """
        copied_files = []

        for original_path in self.src_file_inputs:

            if not os.path.exists(original_path):
                raise ValueError(f"Input path not found or inaccessible: {original_path}")

            staged_path = pathlib.Path(f"{destination}{original_path}")
            dest_path = staged_path.as_posix()
            copied_files.append(dest_path)

            if os.path.isdir(original_path):
                raise ValueError(
                    f"Invalid {original_path}; directories not supported for File types; use a tarball/zip instead."
                )

            dirname = pathlib.Path(os.path.dirname(staged_path))
            dirname.mkdir(mode=0o0770, parents=True, exist_ok=True)

            # Files must be copied in to ensure they are readable by the "jaws" linux user.  The group
            # will be set correctly as a result of the gid sticky bit and acl rules on the inputs dir.
            copy_with_progress_bar(original_path, dest_path, quiet=quiet)
            try:
                os.chmod(dest_path, 0o0666)
            except Exception as error:
                raise IOError(f"Error chmod {dest_path}: {error}")

        return copied_files

    def prepend_paths_to_json(self, staging_dir):
        """
        Modifies the paths of real files by prepending the inputs dir to where they were copied.
        Thus the inputs-json will have file paths which point to the new (copied) location that JAWS can access.
        Only paths which exist are modified; other path-like strings presumably refer to files in the container.
        """
        new_json = {}
        for k in self.inputs_json:
            new_json[k] = apply_filepath_op(
                self.inputs_json[k],
                self.wom_input_types[k],
                self._prepend_path(staging_dir),
            )
        self.inputs_json = new_json
        return new_json

    @staticmethod
    def _prepend_path(path_to_prepend):
        """
        The input files shall be transferred to the compute-site, so the paths in the inputs json file must be
        updated to reflect the new location by prepending the site's inputs basedir.
        Only paths which exist are modified because these are the only files which are transferred; the other
        paths presumably refer to files in a Docker container.
        If element is a file and is an URL or in /refdata, do not prepend given path to it
        """

        def func(path, is_file):
            if is_file:
                if not path.startswith(
                    ("http://", "ftp://", "https://", REFDATA_DIR)
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
            json.dump(self.inputs_json, json_inputs, indent=2)


class Manifest:
    """
    The Manifest is a list of all the Run's input files which must be transferred to the compute jaws-site.
    Paths stored are relative to the specified base directory.
    """

    def __init__(self, basedir):
        self.basedir = basedir

    def add(self, *files):
        # use a set to ensure files are unique
        uniq_files = set()
        for filepath in files:
            if filepath and os.path.isfile(filepath):
                abs_path = os.path.abspath(filepath)
                rel_path = os.path.relpath(abs_path, self.basedir)
                uniq_files.add(rel_path)
        self.files = list(uniq_files)


class WorkflowError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WdlError(Exception):
    def __init__(self, message):
        super().__init__(message)


class WorkflowInputsError(Exception):
    def __init__(self, message):
        super().__init__(message)
