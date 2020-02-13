"""
Class for WDL (WOCON) and inputs (JSON).
"""

import sys
import os
import shutil
import json
import subprocess
from subprocess import Popen, PIPE, STDOUT
import re
import zipfile

DEBUG = False
if "JAWS_DEBUG" in os.environ: DEBUG = True

if "WOMTOOL" not in os.environ: sys.exit("\$WOMTOOL env var not defined")
WOMTOOL = os.environ["WOMTOOL"]

class Workflow:

    wdl_file = None
    json_file = None
    zip_file = None
    subworkflows = None
    inputs_file = None
    inputs_dict = None
    manifest_file = None
    max_ram_gb = 0

    def __init__(self, wdl_file, inputs_file):
        if not os.path.exists(wdl_file): sys.exit("wdl_file not found: %s" % (wdl_file,))
        if not os.path.exists(inputs_file): sys.exit("inputs_file not found: %s" % (inputs_file,))
        self.wdl_file = wdl_file
        self.inputs_file = inputs_file

    def validate(self):
        """
        Perform all tests.  Returns True upon success; False otherwise.
        """
        if not self.validate_workflow(): return False
        if not self.validate_inputs(): return False
        return True


    #########################
    ## WORKFLOW (WDL) METHODS
    ## including subworkflows

    def validate_workflow(self):
        """
        Validate main WDL file and identify any subworkflows.  Returns True on success, False otherwise.
        """
        # VALIDATE MAIN WDL
        if DEBUG: print("Validating workflow")
        proc = subprocess.run(["java", "-jar", WOMTOOL, "validate", "-l", self.wdl_file], capture_output=True, text=True)
        subworkflows = set(proc.stdout.splitlines())
        if "Success!" in subworkflows: subworkflows.remove("Success!")
        if "List of Workflow dependencies is:" in subworkflows: subworkflows.remove("List of Workflow dependencies is:")
        if "None" in subworkflows: subworkflows.remove("None")
        self.subworkflows = subworkflows

        # ANY MISSING SUBWORKFLOWS?
        missing = set()
        for line in proc.stderr.splitlines():
            m = re.match("Failed to import workflow (.+).:", line)
            if m:
                for sub in m.groups(): missing.add(sub)
        if missing:
            for sub in missing: sys.stderr.write("Subworkflow file not found: %s\n" % (sub,))
            return False
        return True

    def _identify_file_parameters(self):
        """
        Validate the WDL using Cromwell's womtool and define the set of parameters of "File" type as self.file_items.  Used by load_inputs_json()
        """
        if DEBUG: print("Identifying \"File\" parameters")
        items = {}
        proc = subprocess.Popen(["java", "-jar", WOMTOOL, "inputs", self.wdl_file], stdout=subprocess.PIPE)
        stdout = proc.communicate()[0].decode('utf8')
        wom_output = json.loads(stdout)
        for key in wom_output.keys():
            value = wom_output[key]
            if value.startswith('File'):
                items[key] = 'File'
            elif value == 'Array[File]':
                items[key] = 'List'
            elif value == 'Map[File, File]':
                items[key] = 'File:File'
            elif value.startswith('Map[File,'):
                items[key] = 'File:Other'
            elif value.startswith('Map[') and value.endswith(', File]'):
                items[key] = 'Other:File'
        self.file_items = items
        if not items: print("[WARNING] No input files specified")
        return self.file_items


    def _filter_wdl(self, infile, outfile):
        """
        Removes any "backend" tags in the WDL file.  Also updates the self.max_ram_gb value.  Returns True on success; False otherwise.
        We disallow "backend" tags in the WDL since they override the server defaults and may overload the local server.  Note that this tag must be entirely contained within a single line as the WDL is not fully parsed with python and the java womtool doesn't have this feature.  This tag is currently only required when running on AWS.
        We determine the maximum amount of requested memory (in gigabytes) across all WDLs (i.e. including subworkflows) to ensure the cluster has nodes capable of running the workflow (otherwise the job shall never be run by the scheduler).
        """
        assert(infile)
        assert(outfile)

        # READ WDL
        if DEBUG: print("Reading WDL: %s" % (infile,))
        lines = None
        try:
            with open(infile, 'r') as f:
                lines = f.readlines()
        except:
            sys.stderr.write("Error reading WDL file: %s" % (infile,))
            return False
        new_wdl = ""
        for line in lines:
            # DETERMINE MAX REQUESTED RAM
            m = re.match(r"^\s+memory\s*[:=]\s*\"?(\d+\.?\d*)([kKmMgGtT])\"?", line)
            if m:
                mem = m.group(1)
                prefix = m.group(2)
                prefix = prefix.lower()
                if prefix == "k":
                    mem = int(mem / 1048576)
                elif prefix == "m":
                    mem = int(mem / 1024)
                elif prefix == "g":
                    mem = int(mem+0.99) # round up
                elif prefix == "t":
                    mem = int(mem * 1024)
                else:
                    mem = int(mem / 1073741824)
                if mem > self.max_ram_gb: self.max_ram_gb = mem
            # FILTER "backend" TAGS
            m = re.match(r"^\s*backend", line)
            if not m:
                new_wdl = new_wdl + line

        # WRITE (FILTERED) WDL
        if DEBUG: print("Writing WDL: %s" % (outfile,))
        try:
            with open(outfile, 'w') as f:
                f.write(new_wdl)
        except:
            sys.stderr.write("Error writing WDL file: %s" % (outfile,))
        return True


    def prepare_wdls(self, staging_dir, submission_id):
        """
        Copy and filter WDLs (including subworkflows) to new directory. This changes the wdl paths in the object.
        """
        # MAIN WDL
        assert(self.wdl_file)
        assert(staging_dir)
        assert(submission_id)
        main_wdl_file = os.path.join(staging_dir, "%s.wdl" % (submission_id,))
        self._filter_wdl(self.wdl_file, main_wdl_file)
        self.wdl_file = main_wdl_file

        # PREPARE SUBWORKFLOWS
        if self.subworkflows:
            # CREATE FOLDER FOR (FILTERED) SUBWORKFLOWS
            dest = os.path.join(staging_dir, submission_id)
            os.mkdir(dest, 0o0777)
            self.zip_file = os.path.join(staging_dir, "%s.zip" % (submission_id,))
            new_subworkflows = set()
            for sub_wdl_src in self.subworkflows:
                sub_name = os.path.basename(sub_wdl_src)
                sub_wdl_dest = os.path.join(dest, sub_name)
                self._filter_wdl(sub_wdl_src, sub_wdl_dest)
                new_subworkflows.add(sub_wdl_dest)
            self.subworkflows = new_subworkflows 
            if not self.zip_subworkflows(): return False
            shutil.rmtree(dest)
        return True

    def zip_subworkflows(self):
        """
        Create a zip file of the subworkflows, as required by Cromwell.
        """
        if not self.subworkflows: return True
        assert(self.zip_file)
        if os.path.exists(self.zip_file): os.remove(self.zip_file)
        print("Zipping subworkflows to %s" % (self.zip_file,))
        try:
            with zipfile.ZipFile(self.zip_file, 'w') as z:
                for sub_wdl in self.subworkflows:
                    z.write(sub_wdl, arcname=os.path.basename(sub_wdl))
        except:
            sys.stderr.write("Error writing subworkflows zipfile: %s" % (self.zip_file,))
            if os.path.exists(self.zip_file): os.remove(self.zip_file)
            return False
        return True


    ########################
    ## INPUTS (JSON) METHODS

    def validate_inputs(self):
        """
        Load inputs json file, converts all paths to absolute, and verifies they exist.
        This method populates:
            self.inputs_dict  (parameter key => value)
            self.source_files (set of absolute paths)
        Returns True if successful; False otherwise.  Does not write any outfile.
        """
        if DEBUG: print("Validating inputs JSON file")
        try:
            with open(self.inputs_file, 'r') as f:
                self.inputs_dict = json.load(f)
        except:
            if DEBUG: sys.stderr.write("Unable to load inputs json\n")
            return False

        # CONVERT RELATIVE PATHS TO ABSOLUTE PATHS
        source_files = set() 
        if not self._identify_file_parameters():
            # no "File" parameters
            self.source_files = source_files
            return True
        items = self.file_items
        if DEBUG: print("Making all paths absolute")
        inputs_dict = self.inputs_dict
        has_rel_paths = False
        json_dirname = os.path.dirname(self.inputs_file)
        for key in inputs_dict.keys():
            value = inputs_dict[key]
            if key in items:
                t = items[key]
                if t == 'File':
                    fo = value
                    if not os.path.isabs(fo):
                        has_rel_paths = True
                        fo = os.path.join(json_dirname, fo)
                    fa = os.path.abspath(fo)
                    source_files.add(fa)
                    inputs_dict[key] = fa
                elif t == 'List':
                    new_value = []
                    for fo in value:
                        if not os.path.isabs(fo):
                            has_rel_paths = True
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        source_files.add(fa)
                        new_value.append(fa)
                    inputs_dict[key] = new_value
                elif t == 'File:File':
                    new_value = {}
                    for fo1 in value:
                        if not os.path.isabs(fo1):
                            has_rel_paths = True
                            fo1 = os.path.join(json_dirname, fo1)
                        fa1 = os.path.abspath(fo1)
                        source_files.add(fa1)
                        fo2 = value[fo1]
                        if not os.path.isabs(fo2):
                            has_rel_paths = True
                            fo2 = os.path.join(json_dirname, fo2)
                        fa2 = os.path.abspath(fo2)
                        source_files.add(fa2)
                        new_value[fa1]=fa2
                    inputs_dict[key] = new_value
                elif t == 'File:Other':
                    new_value = {}
                    for fo in value:
                        if not os.path.isabs(fo):
                            has_rel_paths = True
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        source_files.add(fa)
                        new_value[fa] = value[fo]
                    inputs_dict[key] = new_value
                elif t == 'Other:File':
                    new_value = {}
                    for o in value:
                        fo = value[o]
                        if not os.path.isabs(fo):
                            has_rel_paths = True
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        source_files.add(fa)
                        new_value[o] = fa
                    inputs_dict[key] = new_value
        self.source_files = source_files
        return self._validate_infiles()
        

    def _validate_infiles(self):
        """
        Verify that all infiles exist, are readable, and are files rather than folders.
        All files are checked, even if there are errors (doesn't quit on first error).
        Returns True if OK, False otherwise.
        """
        is_okay = True
        if DEBUG: print("Validating input files")
        if len(self.source_files) == 0:
            if DEBUG: print("[WARNING] Workflow doesn't have any input files")
        for a_file in self.source_files: # abs paths
            if not os.path.exists(a_file):
                sys.stderr.write("[ERROR] File not found: %s\n" % (a_file,))
                is_okay = False
            elif not os.access(a_file, os.R_OK):
                sys.stderr.write("[ERROR] File not readable: %s\n" % (a_file,))
                is_okay = False
            elif os.path.isdir(a_file):
                sys.stderr.write("[WARNING] Is a dir, not a file: %s\n" % (a_file,))
            #elif DEBUG: print("\t%s" % (a_file,))
        return is_okay


    def write_inputs_json(self):
        """
        Write inputs JSON to specified outfile (e.g. after manipulating paths).
        Returns True on success; False otherwise.
        """
        assert(self.json_file)
        if DEBUG: print("Writing inputs JSON: %s" % (self.json_file,))
        if self.inputs_dict is None:
            if DEBUG: sys.stderr.write("\t[ERROR] Inputs not defined\n")
            return False
        try:
            with open(self.json_file, 'w') as f:
                json.dump(self.inputs_dict, f, indent=4)
        except:
            sys.stderr.write("[ERROR] Failed to write inputs JSON file\n")
            return False
        return True


    def _update_input_paths(self, infiles):
        """
        Update self.input_dict using mapping information provided in infiles dictionary (source => dest).
        """
        items = self.file_items
        if DEBUG: print("Updating paths in inputs JSON")
        inputs_dict = self.inputs_dict
        for key in inputs_dict.keys():
            value = inputs_dict[key]
            if key in items:
                t = items[key]
                if t == 'File':
                    source = value
                    if source in infiles:
                        inputs_dict[key] = infiles[source]
                elif t == 'List':
                    new_value = []
                    for source in value:
                        new_value.append(infiles[source])
                    inputs_dict[key] = new_value
                elif t == 'File:File':
                    new_value = {}
                    for source1 in value:
                        dest1 =infiles[source1]
                        source2 = value[source1]
                        dest2 = infiles[source2]
                        new_value[dest1]=dest2
                    inputs_dict[key] = new_value
                elif t == 'File:Other':
                    new_value = {}
                    for source in value:
                        dest = infiles[source]
                        new_value[dest] = value[source]
                    inputs_dict[key] = new_value
                elif t == 'Other:File':
                    new_value = {}
                    for o in value:
                        source = value[o]
                        dest = infiles[source]
                        new_value[o] = dest
                    inputs_dict[key] = new_value


    def prepare_inputs(self, globus_basedir, staging_subdir, site_name, submission_id):
        """
        Copy or symlink (depending on path) all infiles to the staging directory so they may be transferred via Globus.
        Also calculates total gigabytes to be transferred and generates the file manifest for sending via Globus.
        """
        # VALIDATE FOLDERS
        staging_dir = os.path.join(globus_basedir, staging_subdir)
        data_dir = os.path.join(staging_dir, site_name)
        if not os.path.isdir(data_dir): os.makedirs(data_dir, 0o0777)
        if not os.path.isdir(staging_dir): os.makedirs(staging_dir, 0o0777)

        # COPY OR SYMLINK FILES
        if DEBUG: print("Staging infiles")
        infiles = {} # source => dest
        new_source_files = set()
        transfer_mb = 0
        manifest = [] # files to xfer
        for source in self.source_files:
            new_source = source
# DEPRECATED
#            # MIRRORED REFERENCE DATA FOLDERS NO LONGER NEEDED -- CACHING ALL FILES NOW
#            if local_ref_dir and remote_ref_dir and source.startswith(local_ref_dir):
#                new_source = source.replace(local_ref_dir, remote_ref_dir, 1)
#                if DEBUG: print("\t%s" % (source,))
# /DEPRECATED
            dest = "%s/%s" % (data_dir, source) # os.path.join won't work because "source" is an abs path
            if DEBUG: print("\t%s" % (source,))
            dirname = os.path.dirname(dest)
            if not os.path.isdir(dirname): os.makedirs(dirname, 0o0777)
            if source.startswith(globus_basedir):
                if not os.path.exists(dest): os.symlink(source, dest)
            else:
                if not os.path.exists(dest) or int(os.path.getmtime(source)) != int(os.path.getmtime(dest)):
                    shutil.copy2(source, dest, follow_symlinks=True)
            new_source = "./%s%s" % (site_name, source)
            remote_dest = "%s%s" % (site_name, source)
            transfer_mb = transfer_mb + int(os.path.getsize(os.path.realpath(source))/1048576 + 0.5)
            manifest.append([dest, remote_dest])
            infiles[source] = new_source # mapping dict for updating inputs JSON
            new_source_files.add(new_source) # add to set of infiles
        self.source_files = new_source_files
        self.transfer_gb = int(transfer_mb/1024 + 0.5)
        self.manifest = manifest

        # CHANGE INFILE PATHS TO NEW LOCATION
        self._update_input_paths(infiles)

        # WRITE JSON INPUTS FILE
        self.json_file = os.path.join(staging_dir, "%s.json" % (submission_id,))
        if not self.write_inputs_json(): sys.exit("FAILED to write new inputs json")

        # ADD WORKFLOW FILES TO MANIFEST
        wdl_basename = os.path.basename(self.wdl_file)
        self.manifest.append([self.wdl_file, wdl_basename])
        if self.zip_file:
            zip_basename = os.path.basename(self.zip_file)
            self.manifest.append([self.zip_file, zip_basename])
        json_basename = os.path.basename(self.json_file)
        self.manifest.append([self.json_file, json_basename])


    def write_manifest(self, staging_dir, submission_id):
        """
        Write manifest.tsv file
        """
        outfile = self.manifest_file = os.path.join(staging_dir, "%s.tsv" % (submission_id,))
        if DEBUG: print("Writing manifest file: %s" % (outfile,))
        try:
            with open(outfile, 'w') as f:
                for src, dest in self.manifest:
                   f.write("%s\t%s\n" % (src, dest)) 
        except:
            sys.stderr.write("Error writing WDL file: %s" % (outfile,))
