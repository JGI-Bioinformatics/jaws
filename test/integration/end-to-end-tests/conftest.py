# content of conftest.py
import json
import pytest
import smtplib
from subprocess import Popen, PIPE



def run(cmd):
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             cwd="/global/cscratch1/sd/jaws/jfroula", universal_newlines=True)

    stdout,stderr=output.communicate()
    rc=output.returncode

    return rc,stdout,stderr

@pytest.fixture()
def sleep_little_baby():
    cmd='./go.sh > tmp.txt'
    rc,stdout,stderr = run(cmd)
    return stdout

@pytest.fixture(scope="module")
def submit_wdl():
    """
    This is a fixture that will submit a wdl for all functions to use.  
    This function returns the output of a wdl submission. 
    """
    
    cmd = "jaws run submit /global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/TestsWDLs/fq_count.wdl /global/cscratch1/sd/jfroula/JAWS/jaws/test/integration/end-to-end-tests/TestsWDLs/fq_count.json fq_count_out nersc"
    rc,stdout,stderr = run(cmd)
    print(stderr)

    assert rc == 0
    data = json.loads(stdout)
    return data

@pytest.fixture(scope="module")
def submit_subworkflow():
    """
    This is a fixture that will submit a subworkflow wdl for all functions to use.  
    This function returns the output of the wdl submission. 
    """
    
    cmd = "jaws run submit TestsWDLs/jaws_alignment_example.wdl TestsWDLs/jaws_alignment_example.json alignment_out nersc"
    rc,stdout,stderr = run(cmd)
    print(stderr)

    assert rc == 0
    data = json.loads(stdout)
    return data

#@pytest.fixture(autouse=True)
#def env_setup(monkeypatch):
#    """set jaws-test variables from jaws-test.env"""
#    monkeypatch.setenv('JAWS_CLIENT_CONFIG', '/global/u2/j/jfroula/jaws-test.conf')
#    monkeypatch.setenv('JAWS_CLIENT_LOG', '/global/cscratch1/sd/jaws/jfroula/jaws-test.log')
