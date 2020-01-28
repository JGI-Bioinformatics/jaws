"""
Class for WDL (WOCON) and inputs (JSON).
"""

import sys
import os
import json
import subprocess
import re
import zipfile

class Workflow:

    verbose = True
    womtool = None
    wdl_file = None
    inputs_file = None
    inputs_json = None

    def __init__(self, womtool, wdl_file, inputs_file, verbose=True):
        if not os.path.exists(womtool): sys.exit("womtool not found: %s" % (womtool,))
        if not os.path.exists(wdl_file): sys.exit("wdl_file not found: %s" % (wdl_file,))
        if not os.path.exists(inputs_file): sys.exit("inputs_file not found: %s" % (inputs_file,))
        self.womtool = womtool
        self.wdl_file = wdl_file
        self.inputs_file = inputs_file

    def identify_file_parameters(self):
        """
        Validate the WDL using Cromwell's womtool and define the set of parameters of "File" type.
        """
        if self.verbose: print("Identifying File parameters")
        items = {}
        proc = subprocess.Popen(["java", "-jar", self.womtool, "inputs", self.wdl_file], stdout=subprocess.PIPE)
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
        if self.verbose:
            if len(items) == 0:
                print("\tNo files specified")
            else:
                for param in items:
                    print("\t%s" % (param,))
        return self.file_items

    def load_inputs_json(self):
        """
        Load inputs json file and save in object.
        """
        if self.verbose: print("Loading inputs json file")
        try:
            with open(self.inputs_file, 'r') as f:
                self.inputs_json = json.load(f)
            if self.verbose: print("\tOK")
            return True
        except:
            if self.verbose: sys.stderr.write("Unable to load inputs json\n")
            return False
        
    def write_inputs_json(self, outfile=None):
        """
        Overwrite inputs json file (e.g. use when file paths changed from relative to absolute).
        """
        if self.verbose: print("Writing inputs json file")
        if self.inputs_json is None:
            if self.verbose: sys.stderr.write("\tInputs not defined\n")
            return False
        if outfile is None:
            outfile = self.inputs_file
        try:
            with open(outfile, 'w') as f:
                json.dump(self.inputs_json, f, indent=4)
            if self.verbose: print("\tOK")
            return True
        except:
            sys.stderr.write("Error overwriting inputs json file\n")
            return False

    def rel_to_abs(self, overwrite=True):
        """
        Convert all file paths to absolute ones.  By default, if the original contains any relative paths, the file with overwritten after the transformation.  Returns True on success, False on failure, None if there are no file parameters.
        """
        items = self.file_items

        if self.verbose: print("Converting relative input paths to absolute")

        inputs_json = None  # contents of inputs json file
        source_files = set() 
        has_rel_paths = False

        inputs_json = self.inputs_json
        json_dirname = os.path.dirname(self.inputs_file)
        for key in inputs_json.keys():
            value = inputs_json[key]
            if key in items:
                t = items[key]
                if t == 'File':
                    fo = value
                    if not os.path.isabs(fo):
                        has_rel_paths = True
                        fo = os.path.join(json_dirname, fo)
                    fa = os.path.abspath(fo)
                    source_files.add(fa)
                    inputs_json[key] = fa
                elif t == 'List':
                    new_value = []
                    for fo in value:
                        if not os.path.isabs(fo):
                            has_rel_paths = True
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        source_files.add(fa)
                        new_value.append(fa)
                    inputs_json[key] = new_value
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
                    inputs_json[key] = new_value
                elif t == 'File:Other':
                    new_value = {}
                    for fo in value:
                        if not os.path.isabs(fo):
                            has_rel_paths = True
                            fo = os.path.join(json_dirname, fo)
                        fa = os.path.abspath(fo)
                        source_files.add(fa)
                        new_value[fa] = value[fo]
                    inputs_json[key] = new_value
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
                    inputs_json[key] = new_value
        self.source_files = source_files
        if len(source_files) == 0:
            if self.verbose: print("\tNo files specified")
            return None
        elif has_rel_paths:
            if self.verbose: print("\tRelative paths found")
            if overwrite is False:
                self.inputs_file = self.inputs_file + ".abs"
            return self.write_inputs_json()
        elif self.verbose:
            print("\tNo relative paths found")
            return True

    def validate_infiles(self):
        """
        Verify that all infiles exist, are readable, and are files (not folders).  Returns True if pass, None if there aren't any infiles, and False if any bad infiles found.
        """
        if self.verbose: print("Validating input files")
        if len(self.source_files) == 0:
            if self.verbose: print("\tWorkflow doesn't have input files")
            return None
        result = True
        for fa in self.source_files:
            if not os.path.exists(fa):
                sys.stderr.write("File not found: %s\n" % (fa,))
                result = False
            elif not os.access(fa, os.R_OK):
                sys.stderr.write("File unreadable: %s\n" % (fa,))
                result = False
            elif os.path.isdir(fa):
                sys.stderr.write("Is a dir, not a file: %s\n" % (fa,))
            elif self.verbose:
                print("\t%s" % (fa,))
        return result

    def validate_subworkflows(self):
        """
        Use Cromwell's "womtool" to determine subworkflow paths. Returns None if there are no subworkflows, True if all subworkflow files exist, False otherwise.  Paths are stored in self.subworkflows.
        """
        if self.verbose: print("Validating subworkflows")
        proc = subprocess.run(["java", "-jar", self.womtool, "validate", "-l", self.wdl_file], capture_output=True, text=True)
        subworkflows = set(proc.stdout.splitlines())
        if "Success!" in subworkflows:
            subworkflows.remove("Success!")
        if "List of Workflow dependencies is:" in subworkflows:
            subworkflows.remove("List of Workflow dependencies is:")
        if "None" in subworkflows:
            subworkflows.remove("None")
        self.subworkflows = subworkflows
        if self.verbose:
            for sub in subworkflows:
                print("\t%s" % (sub,))

        missing = set()
        for line in proc.stderr.splitlines():
            m = re.match("Failed to import workflow (.+).:", line)
            if m:
                for sub in m.groups():
                    missing.add(sub)
        if len(missing):
            for sub in missing:
                sys.stderr.write("Subworkflow file not found: %s\n" % (sub,))
            return False
        elif len(subworkflows) == 0:
            if self.verbose: print("\tNo subworkflows specified")
            return None
        else:
            return True

    def zip_subworkflows(self, zip_file=None):
        """
        Create a zip archive with all of the workflow's required subworkflows.  These WDLs are all saved in the `cwd` of the archive, regardless of their original location.  Returns True upon success, False if failed, None if there aren't any subworkflows.  If the outfile exists, it is overwritten.
        """
        if not self.subworkflows and len(self.subworkflows):
            return None
        if zip_file is None: zip_file = re.sub('\.json$', '', self.wdl_file) + ".sub.zip"
        if os.path.exists(zip_file): os.remove(zip_file)
        if self.verbose: print("Zipping subworkflows")
        try:
            with zipfile.ZipFile(zip_file, 'w') as z:
                for s in self.subworkflows:
                    z.write(s, arcname=os.path.basename(s))
            if self.verbose: print("\t%s" % (zip_file,))
            return True
        except:
            sys.stderr.write("Error writing subworkflows zipfile: %s" % (zip_file,))
            if os.path.exists(zip_file): os.remove(zip_file)
            return False

    def validate(self):
        """
        Perform all tests but do not overwrite inputs json with absolute file paths.  Exit upon error; return otherwise.
        """
        result = self.validate_subworkflows()
        if result is False:
            sys.exit("FAIL")
        elif result is True:
            if self.zip_subworkflows() is False: sys.exit("FAIL")
        if self.identify_file_parameters() is False: sys.exit("FAIL")
        if self.load_inputs_json() is False: sys.exit("FAIL")
        if self.rel_to_abs(overwrite=False) is False: sys.exit("FAIL")
        if self.validate_infiles() is False: sys.exit("FAIL")
        print("PASS")
