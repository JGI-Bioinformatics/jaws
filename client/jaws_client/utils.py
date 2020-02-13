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

from jaws_client import wfcopy as wfc
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


# CONTRIBUTED/MAINTAINED BY JEFF FROULA
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


# CONTRIBUTED/MAINTAINED BY STEPHAN TRONG
@util.command()
@click.argument("cromwelldir")
@click.argument("outdir")
@click.option('--flattenShardDir', default=False)
def wfcopy(cromwelldir, outdir, flattensharddir):
    """
    Copy cromwell output to specified dir.
    """
    wfc.wfcopy(cromwelldir, outdir, flattenShardDir=flattensharddir)


if __name__ == "__main__":
    util()
