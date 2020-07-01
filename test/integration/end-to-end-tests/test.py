#!/usr/bin/env python
from subprocess import Popen, PIPE
from time import sleep
import os

def run(cmd):
    output = Popen(cmd, stdout=PIPE,
             stderr=PIPE, shell=True,
             universal_newlines=True)
             #cwd="/global/cscratch1/sd/jaws/jfroula", universal_newlines=True)

    stdout,stderr=output.communicate()
    rc=output.returncode

    return rc,stdout,stderr

def test_sleep(sleep_little_baby):
    stdout=0
    while (int(stdout) < 4):
        if os.path.exists("tmp.txt"):
            sleep(5)
            cmd='cat tmp.txt'
            rc,stdout,stderr = run(cmd)
            print(f"test_sleep {stdout}")
        else:
            sleep(5)

def mytest_next_sleep(sleep_little_baby):
    cmd='cat tmp.txt'
    rc,stdout,stderr = run(cmd)
    print(f"test_next_sleep {stdout}")
