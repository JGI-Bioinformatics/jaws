"""
Functions for WDL files
"""

# HISTORY:
# - 2019 : created by Jeff Froula

import sys
import os
import re
import subprocess
import json
import collections
import stat

WOMTOOL = os.environ["WOMTOOL"]

def womtool(wdl):
    cmd="java -jar {} validate {}".format(WOMTOOL, wdl),

    process = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    outerr = process.stderr.decode('utf-8')

    if (not outerr):
        return 1
    else:
        print("WDL validation failed because.\n{}".format(outerr))
        return 0


def _flatten(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(_flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def _checkPathPermissions(filepath):
    '''
    checks if a file's parent directories have world-executable permissions and 
    the file itself has world-readable permissions.  It returns a list
    of directories (in the file's path) that don't
    '''
    previous_dir = ''

    print("-- testing file: {}".format(filepath))
    dir_list = filepath.split('/')
    for d in dir_list:
        if (d):
            previous_dir += '/' + d
            if (os.path.isdir(previous_dir)):
                # use stat to check if a directory has wold readable permissions
                st = os.stat(previous_dir)
                if (not bool(st.st_mode & stat.S_IXOTH)):
                    print("     Permissions failed for the following directory: {}".format(previous_dir))
                    return 0
            if (os.path.isfile(previous_dir)):
                st = os.stat(filepath)
                if (not bool(st.st_mode & stat.S_IROTH)):
                    print("     Permissions failed for the following file: {}".format(previous_dir))
                    return 0
    if (os.path.exists(filepath)):
        print("     OK")
    else:
        print("     File doesn't exist")

    return 1

def inputFilePermissions(wdl,input_json):
    '''
    This finds a user's file path from the inputs.json and checks that it is readable by everyone.
    Each parent directory must also have appropriate world-executable permissions.
    File paths are detected in the inputs.json by the presence of back slashes.
    Assumption: the file paths must be full paths.
    '''

    # lets first get a list of "Files" from input.json
    with open(input_json) as j:
        data = json.load(j)

    flat_dict = _flatten(data)
    files_to_check = []
    for k,v in flat_dict.items():
        if isinstance(v, list):
            for i in v:
                if '/' in i:
                    files_to_check.append(i)
        if isinstance(v, tuple):
            for i in v:
                if '/' in i:
                    files_to_check.append(i)
        if isinstance(v, str):
            if '/' in v:
                files_to_check.append(v)


    failed=0
    for filepath in files_to_check:
        filepath=filepath.strip()
        if(_checkPathPermissions(filepath)):
            failed=1

    if (failed):
        return 0
    else:
        return 1


def searchCommandShifter(wdl):
    '''
    This tests that the word shifter or docker is not being used to run a command. 
    The backend.providers sections of the config file will determine whether its a shifter of docker call.
    '''
    f = open(wdl, 'r')
    text=f.read()
    search_for_shifter = re.search(r"^\s*shifter", text, re.MULTILINE)
    if search_for_shifter:
        print("Error: You should not call shifter in the command stanza. "
              "You should instead set the image in a runtime stanza like this:\n"
              "runtime { docker: <your_image> } ")
        return 0


    print("testing if docker is part of command ...")
    search_for_docker = re.search(r"^\s*docker", text, re.MULTILINE)
    if search_for_docker:
        print("Error: You should not call docker in the command stanza. "
              "You should instead set the image to use in a runtime stanza like follows:\n"
              "runtime { docker: <your_image> } ")
        return 0

def checkRuntimeFormat(wdl):
    '''
    This checks that the runtime attribute section of the wdl has excepted key words 
    '''
    accepted_keys = ["mem","cpu","time","poolname","poolsize","cluster","nwpn"]
    f = open(wdl, 'r')
    text=f.read()

    # matches the whole runtime stanza
    runtime_matches = re.findall(r"runtime\s*{(.*?)(?:}\s*}|output|command)", text, re.DOTALL)
    for match in runtime_matches:
        # print(match)

        # test if just an emtpy runtime stanza which is OK
        found_something = re.search(r"\S", match, re.DOTALL)
        if not found_something:
            print("Warning: There was a runtime stanza that was empty. Will continue anyway.")
            continue

        # this part we'll be testing that the key names conform to JGI standards(i.e. they must
        # be in the "accepted_keys" array.
        get_keys = re.findall(r"(\w+):\s", match)

        # if runtime stanza is not emtpy and not valid key structure found
        if not get_keys:
            print("There doesn't seem to be any valid variables specified in the runtime stanza. "
                  "You need to use colon after key word.\n{}".format(match))
            return 0

        # check that key names are correct
        for key in get_keys:
            if key not in accepted_keys:
                print("error: {} is not an excepted key. please make sure "
                      "you use only the following keys in the runtime stanza:\n {}".format(key,', '.join(accepted_keys)))
                return 0

    return 1


