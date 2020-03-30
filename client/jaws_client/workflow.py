"""
Class for WDL (WOCON) and inputs (JSON).
"""

import os
import shutil
import json
import subprocess
import re
import zipfile
import uuid
import logging
from jaws_client import config


def rsync(src, dest):
    """Copy source to destination using rsync.

    :param src: Source path
    :type src: str
    :param dest: Destination path
    :type dest: str
    :return:
    """
    cmd = f"rsync -rLptq {src} {dest}"
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode:
        raise IOError(f"Failed to rsync {src} to {dest}: " + stderr.strip())


class WorkflowError(Exception):
    def __init__(self, message):
        super().__init__(message)


class Workflow:
    """The Workflow class is used to validate a WDL and manipulate the inputs JSON"""

    wdl_file = None
    json_file = None
    subworkflows = None
    inputs_file = None
    inputs_dict = None
    manifest_file = None
    max_ram_gb = 0

    def __init__(self, wdl_file: str, inputs_file: str) -> None:
        """Constructor will validate the WDL and inputs and raise an exception if invalid.

        :param wdl_file: Path to the main workflow specification (WDL) file
        :type wdl_file: str
        :param inputs_file: Path to the inputs (JSON) file
        :type inputs_file:
        :return:
        """
        self.logger = logging.getLogger(__package__)
        if not os.path.exists(wdl_file):
            raise IOError(f"wdl_file not found: {wdl_file}")
        if not os.path.exists(inputs_file):
            raise IOError(f"inputs_file not found: {inputs_file}")
        self.wdl_file = wdl_file
        self.inputs_file = inputs_file

    def validate(self):
        """Validate a workflow; raise on error"""
        self.validate_wdls()
        self.validate_inputs()

    def validate_wdls(self):
        """Validate main WDL and subworkflows."""
        self.logger.info(f"Validating WDL, {self.wdl_file}")
        proc = subprocess.run(
            [
                "java",
                "-jar",
                config.conf.get("JAWS", "womtool_jar"),
                "validate",
                "-l",
                self.wdl_file,
            ],
            capture_output=True,
            text=True,
        )
        if proc.stderr:
            raise WorkflowError(proc.stderr)
        self.subworkflows = set(proc.stdout.splitlines())
        if "Success!" in self.subworkflows:
            self.subworkflows.remove("Success!")
        if "List of Workflow dependencies is:" in self.subworkflows:
            self.subworkflows.remove("List of Workflow dependencies is:")
        if "None" in self.subworkflows:
            self.subworkflows.remove("None")
        missing = set()
        for line in proc.stderr.splitlines():
            m = re.match("Failed to import workflow (.+).:", line)
            if m:
                for sub in m.groups():
                    missing.add(sub)
        if missing:
            raise WorkflowError("Subworkflows not found: " + ", ".join(missing))

    def max_ram_gb(self):
        """Determine the maximum amount of RAM requested by the workflow.

        :return: maximum RAM in GB, rounded up to the nearest GB.
        :rtype: int
        """
        self.max_ram_gb = 0
        self._max_ram_gb(self.wdl_file)
        for subworkflow_file in self.subworkflows:
            self._max_ram_gb(subworkflow_file)
        self.logger.info(f"Maximum RAM requested is {self.max_ram_gb}Gb")
        return self.max_ram_gb

    def _max_ram_gb(self, infile):
        """Check a WDL file and update the self.max_ram_gb value if new max."""
        with open(infile, "r") as f:
            lines = f.readlines()
        for line in lines:
            m = re.match(r"^\s+memory\s*[:=]\s*\"?(\d+\.?\d*)([kKmMgGtT])\"?", line)
            if m:
                mem = m.group(1)
                prefix = m.group(2)
                prefix = prefix.lower()
                if prefix == "k":
                    mem = int(mem / 1048576)
                elif prefix == "m":
                    mem = int(mem / 1024)
                elif prefix == "t":
                    mem = int(mem * 1024)
                else:
                    mem = int(mem / 1073741824)
                mem = int(mem + 0.99)  # round up
                if mem > self.max_ram_gb:
                    self.max_ram_gb = mem

    def validate_inputs(self):
        """Reads inputs JSON, verify files exist, and save their absolute paths in object."""
        self.logger.info(f"Validating inputs file, {self.inputs_file}")
        self.inputs_dict = {}  # parameter key => value
        self.input_files = set()  # set of absolute paths of input files

        # READ JSON FILE
        self.logger.debug("Validating inputs JSON file")
        try:
            with open(self.inputs_file, "r") as f:
                self.inputs_dict = json.load(f)
        except Exception as e:
            raise WorkflowError(f"Unable to load inputs json, {self.inputs_file}: {e}")

        # IDENTIFY FILES, RETURN IF NONE
        self._identify_file_parameters()
        if not self.file_items:
            return

        # CONVERT RELATIVE PATHS TO ABSOLUTE
        json_dirname = os.path.dirname(
            self.inputs_file
        )  # paths are relative to inputs json file
        for key in self.inputs_dict.keys():
            value = self.inputs_dict[key]
            if key in self.file_items:
                t = self.file_items[key]
                if t == "File":
                    fo = value
                    if not os.path.isabs(fo):
                        fo = os.path.join(json_dirname, fo)
                    fa = os.path.abspath(fo)
                    self.input_files.add(fa)
                    self.inputs_dict[key] = fa
                elif t == "List":
                    new_value = []
                    for fo in value:
                        if not os.path.isabs(fo):
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        self.input_files.add(fa)
                        new_value.append(fa)
                    self.inputs_dict[key] = new_value
                elif t == "File:File":
                    new_value = {}
                    for fo1 in value:
                        if not os.path.isabs(fo1):
                            fo1 = os.path.join(json_dirname, fo1)
                        fa1 = os.path.abspath(fo1)
                        self.input_files.add(fa1)
                        fo2 = value[fo1]
                        if not os.path.isabs(fo2):
                            fo2 = os.path.join(json_dirname, fo2)
                        fa2 = os.path.abspath(fo2)
                        self.input_files.add(fa2)
                        new_value[fa1] = fa2
                    self.inputs_dict[key] = new_value
                elif t == "File:Other":
                    new_value = {}
                    for fo in value:
                        if not os.path.isabs(fo):
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        self.input_files.add(fa)
                        new_value[fa] = value[fo]
                    self.inputs_dict[key] = new_value
                elif t == "Other:File":
                    new_value = {}
                    for o in value:
                        fo = value[o]
                        if not os.path.isabs(fo):
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        self.input_files.add(fa)
                        new_value[o] = fa
                    self.inputs_dict[key] = new_value

        # VALIDATE INPUT FILES
        bad_input = []
        for a_file in self.input_files:
            if not os.path.exists(a_file) or not os.access(a_file, os.R_OK):
                bad_input.append(a_file)
        if len(bad_input):
            raise WorkflowError(
                "Missing/unreadable input files: " + ", ".join(bad_input)
            )

    def _identify_file_parameters(self):
        """Determine which input parameters are of the File type."""
        self.file_items = {}
        proc = subprocess.Popen(
            [
                "java",
                "-jar",
                config.conf.get("JAWS", "womtool_jar"),
                "inputs",
                self.wdl_file,
            ],
            stdout=subprocess.PIPE,
        )
        stdout = proc.communicate()[0].decode("utf8")
        wom_output = json.loads(stdout)
        for key in wom_output.keys():
            value = wom_output[key]
            if value.startswith("File"):
                self.file_items[key] = "File"
            elif value == "Array[File]":
                self.file_items[key] = "List"
            elif value == "Map[File, File]":
                self.file_items[key] = "File:File"
            elif value.startswith("Map[File,"):
                self.file_items[key] = "File:Other"
            elif value.startswith("Map[") and value.endswith(", File]"):
                self.file_items[key] = "Other:File"

    def prepare_submission(self, dest_basedir, dest_staging_subdir):
        """Prepare a workflow for submission to a JAWS-Site."""
        self.basename = str(uuid.uuid4())
        self.globus_basedir = config.conf.get("GLOBUS", "basedir")
        self.staging_dir = os.path.join(
            self.globus_basedir, config.conf.get("JAWS", "staging_subdir"),
        )
        if not os.path.isdir(self.staging_dir):
            os.makedirs(self.staging_dir)
        self.dest_basedir = dest_basedir
        self.dest_staging_subdir = dest_staging_subdir
        self.dest_staging_dir = os.path.join(dest_basedir, dest_staging_subdir)
        self._prepare_wdls()
        self._prepare_infiles()
        self._prepare_inputs_json()
        self._write_manifest()
        return self.basename

    def _prepare_wdls(self):
        """Filter WDLs (including subworkflows) and write to staging dir."""
        self.logger.info(f"Staging WDLs to {self.staging_dir}")
        self.manifest = []  # [[src_abs_path, dest_rel_path, inode_type], ...]

        # FILTER MAIN WDL
        self.staged_wdl_file = os.path.join(self.staging_dir, f"{self.basename}.wdl")
        self.__filter_wdl(self.wdl_file, self.staged_wdl_file)
        self.__add_to_manifest(self.staged_wdl_file)

        # FILTER AND ZIP SUBWORKFLOWS
        if not self.subworkflows:
            self.staged_zip_file = None
            return
        tmpdir = os.path.join(self.staging_dir, self.basename)
        os.mkdir(tmpdir, 0o0777)
        self.staged_zip_file = os.path.join(self.staging_dir, f"{self.basename}.zip")
        self.staged_subworkflows = set()
        for subworkflow_file in self.subworkflows:
            subworkflow_basename = os.path.basename(subworkflow_file)
            staged_subworkflow_file = os.path.join(tmpdir, subworkflow_basename)
            self.__filter_wdl(subworkflow_file, staged_subworkflow_file)
            self.staged_subworkflows.add(staged_subworkflow_file)
        if os.path.exists(self.staged_zip_file):
            os.remove(self.staged_zip_file)
        try:
            with zipfile.ZipFile(self.staged_zip_file, "w") as z:
                for sub_wdl in self.staged_subworkflows:
                    z.write(sub_wdl, arcname=os.path.basename(sub_wdl))
        except Exception as e:
            raise WorkflowError(
                f"Error writing subworkflows zipfile, {self.staged_zip_file}: {e}"
            )
        shutil.rmtree(tmpdir)
        self.__add_to_manifest(self.staged_zip_file)

    def __filter_wdl(self, infile: str, outfile: str) -> None:
        """Remove any disallowed "backend" tags from the WDL file and write.

        :param infile: Path to source WDL file
        :type infile: str
        :param outfile: Path to result WDL file
        :type outfile: str
        :return:
        """
        with open(outfile, "w") as fo:
            with open(infile, "r") as fi:
                for line in fi:
                    m = re.match(r"^\s*backend", line)
                    if not m:
                        fo.write(line)

    def __add_to_manifest(self, src_abs_path):
        """Add a file to the manifest list, used for Globus transfers.
        Destination paths are relative to the endpoint's basedir (may be root)."""
        inode_type = "D" if os.path.isdir(src_abs_path) else "F"
        src_rel_path = os.path.relpath(src_abs_path, self.staging_dir)
        dest_rel_path = f"{self.dest_staging_subdir}/{src_rel_path}"
        self.manifest.append([src_abs_path, dest_rel_path, inode_type])

    def _prepare_infiles(self):
        """Stage infiles and write new inputs JSON file"""
        self.infiles_map = {}  # { source abs path => dest abs path }
        self.dest_input_files = set()
        source_site_id = config.conf.get("JAWS", "site_id")
        data_dir = os.path.join(self.staging_dir, source_site_id)
        self.logger.info(f"Staging infiles to {data_dir}")

        # STAGE ALL INPUTS
        for orig_path in self.input_files:
            staged_path = f"{data_dir}{orig_path}"  # os.path.join won't work because "orig_path" is an abs path
            # MKDIRS AS NECESSARY
            if os.path.isdir(orig_path):
                if not os.path.isdir(staged_path):
                    os.makedirs(staged_path, 0o0700)
            else:
                dirname = os.path.dirname(staged_path)
                if not os.path.isdir(dirname):
                    os.makedirs(dirname, 0o0700)
            # PATHS UNDER BASEDIR ARE ACCESSIBLE BY GLOBUS VIA SYMLINK; OTHERWISE MAKE A COPY
            if orig_path.startswith(self.globus_basedir):
                if not os.path.exists(staged_path):
                    os.symlink(orig_path, staged_path)
            else:
                rsync(orig_path, staged_path)
            self.__add_to_manifest(staged_path)
            dest_path = staged_path.replace(self.staging_dir, self.dest_staging_dir)
            self.infiles_map[
                orig_path
            ] = dest_path  # mapping dict for updating inputs JSON
            self.dest_input_files.add(dest_path)  # add to set of self.infiles_map

    def _prepare_inputs_json(self):
        """Write inputs JSON file with dest paths to staging area."""
        self.staged_json_file = os.path.join(self.staging_dir, f"{self.basename}.json")
        self.dest_inputs_dict = {}
        self.logger.debug(f"Writing prepared inputs JSON to {self.staged_json_file}")
        for key, value in self.inputs_dict.items():
            if key in self.file_items:
                t = self.file_items[key]
                if t == "File":
                    source = value
                    if source in self.infiles_map:
                        self.dest_inputs_dict[key] = self.infiles_map[source]
                elif t == "List":
                    new_value = []
                    for source in value:
                        new_value.append(self.infiles_map[source])
                    self.dest_inputs_dict[key] = new_value
                elif t == "File:File":
                    new_value = {}
                    for source1 in value:
                        dest1 = self.infiles_map[source1]
                        source2 = value[source1]
                        dest2 = self.infiles_map[source2]
                        new_value[dest1] = dest2
                    self.dest_inputs_dict[key] = new_value
                elif t == "File:Other":
                    new_value = {}
                    for source in value:
                        dest = self.infiles_map[source]
                        new_value[dest] = value[source]
                    self.dest_inputs_dict[key] = new_value
                elif t == "Other:File":
                    new_value = {}
                    for o in value:
                        source = value[o]
                        dest = self.infiles_map[source]
                        new_value[o] = dest
                    self.dest_inputs_dict[key] = new_value
            else:
                self.dest_inputs_dict[key] = value
        try:
            with open(self.staged_json_file, "w") as f:
                json.dump(self.dest_inputs_dict, f, indent=4)
        except Exception as e:
            raise WorkflowError(
                f"Failed to write staged JSON, {self.staged_json_file}: {e}"
            )
        self.__add_to_manifest(self.staged_json_file)

    def _write_manifest(self):
        """Write manifest.tsv file"""
        self.manifest_file = os.path.join(self.staging_dir, f"{self.basename}.tsv")
        self.logger.info(f"Writing file manifest to {self.manifest_file}")
        self.logger.debug(f"Writing manifest file: {self.manifest_file}")
        try:
            with open(self.manifest_file, "w") as f:
                for src, dest, inode_type in self.manifest:
                    f.write(f"{src}\t{dest}\t{inode_type}\n")
        except Exception as e:
            raise WorkflowError(
                f"Error writing manifest file, {self.manifest_file}: {e}"
            )
