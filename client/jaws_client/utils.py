"""
Miscellaneous, stand-alone utility functions.  These don't interact with any JAWS systems.
"""

import sys
import os
import click
import re
import subprocess
import pathlib
import json
import pprint
import requests

from jaws_client import user
from jaws_client import wdl_functions as w

JAWS_URL = os.environ["JAWS_URL"]
#CROMWELL = os.environ["CROMWELL"]
#WOMTOOL = os.environ["WOMTOOL"]

@click.group()
def util():
    """
    Workflow developer utilities
    """
    pass

@util.command()
def status():
    """
    Current system status
    """
    current_user = user.User()
    url = "%s/status" % (JAWS_URL,)
    try:
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("JAWS Central is DOWN")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@util.command()
@click.argument("wdl")
@click.argument("input_json")
def validate(wdl, input_json):
    """
    Validate your WDL
    """
    if not os.path.exists(wdl):
        print("Error: %s does not exist" % wdl)
        sys.exit(1)
    if not os.path.exists(input_json):
        print("Error: %s does not exist" % input_json)
        sys.exit(1)

    # run womtool.jar for basic wdl formatting check
    print("### womtool.jar validation ###")
    if (w.womtool(wdl)):
        print("success")
    else:
        print("Error: Womtool.jar validate failed")
    print("#################\n")

    # check that inputfiles all have read permissions (i.e. r-- for all)
    print("### Input file permissions check ###")
    if (w.inputFilePermissions(wdl,input_json)):
        print("success")
    else:
        print("Error: Some files don't have world readable permissions or their upstream directories don't have world-executable permissions")
    print("#################\n")

    # make sure shifter or docker is not called to run commands
    #print("### check shifter/docker not used to call command ###")
    #if (w.searchCommandShifter(wdl)):
        # print("### success ###\n")
    #else:
    #    print("Error: shifter or docker should not be used to call commands")
    #print("#################\n")

    # test that runtime stanza if formatted correctly
    print("### check runtime format ###")
    if (w.checkRuntimeFormat(wdl)):
        print("success")
    else:
        print("Error: runtime keyword check failed")
    print("#################\n")

@util.command()
@click.argument("wdl")
def inputs(wdl):
    """
    Generate inputs template from WDL file.
    """
    cmd = [ "java", "-jar", config["CROMWELL"]["womtool"], "inputs", wdl ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout)
    if stderr:
        sys.stderr.write(stderr)
    sys.exit(process.returncode)

# TODO
#@util.command()
#@click.argument('wdl')
#@click.argument('inputs')
#def run(inputs, wdl):
#    """
#    Run a workflow with local cromwell
#    """
#    sys.exit("Not yet implemented")

# 20190315 - initial version
# 20190805 - when calling a sub-workflow, use the root's task name instead of the sub workflow's task name.
@util.command()
@click.argument("cromwelldir")
@click.argument("outdir")
def wfcopy(cromwelldir, outdir, verbose=False):
    """
    Copy cromwell output to specified dir.
    """
    logdir = os.path.join(outdir, "log")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)

    cromwellFilesToSkip = ['stdout.background', 'stderr.background', 'script.background', 'script.submit']
    taskname = None
    rcfile = os.path.join(logdir, "workflow.rc")
    
    if os.path.exists(rcfile):
        os.remove(rcfile)
    with open(rcfile, 'w') as fh:
        fh.write("#ExitCode\tTask\n")

    for rootdir, subdirs, files in os.walk(cromwelldir):
        if os.path.basename(rootdir).startswith('call-') and rootdir == os.path.join(cromwelldir, os.path.basename(rootdir)):
            taskname = re.sub(r'^call-', '', os.path.basename(rootdir))
            
        if rootdir.endswith('execution'):
            parentdir = str(pathlib.Path(rootdir).parent)
            shardname = None

            if re.search(r'shard-', parentdir):
                shardname = os.path.basename(parentdir)

            taskdir = os.path.join(outdir, taskname)
            if not os.path.exists(taskdir):
                os.makedirs(taskdir)

            for dname in subdirs:
                rsync(os.path.join(rootdir, dname), taskdir, verbose=verbose)

            for fname in files:
                fullname = os.path.join(rootdir, fname)
                outname = "%s-%s"%(taskname, shardname) if shardname else taskname

                if fname == 'stdout':
                    rsync(fullname, os.path.join(logdir, "%s.stdout"%outname), verbose=verbose)
                elif fname == 'stderr':
                    rsync(fullname, os.path.join(logdir, "%s.stderr"%outname), verbose=verbose)
                elif fname == 'script':
                    rsync(fullname, os.path.join(logdir, "%s.script"%outname), verbose=verbose)
                elif fname == 'rc':
                    with open(fullname, 'r') as fh:
                        exitcode = fh.readlines()[0].strip()
                    with open(rcfile, 'a') as fh:
                        fh.write("%s\t%s\n"%(exitcode, outname))
                elif fname not in cromwellFilesToSkip:
                    rsync(fullname, taskdir, verbose=verbose)
                    
def runCommand(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.strip(), stderr.strip(), process.returncode

def rsync(src, dest, verbose=False):
    cmd = "rsync -a %s %s"%(src, dest)
    if verbose:
        print(cmd)
    stdout, stderr, exitcode = runCommand(cmd)
    if exitcode:
        sys.stderr.write(stderr)
        sys.exit(exitcode)

if __name__ == "__main__":
    util()
