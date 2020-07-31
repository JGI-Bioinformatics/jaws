#!/usr/bin/env python

from subprocess import Popen, PIPE

from click.testing import CliRunner
from testcase1 import status

def test_hello_world():
  runner = CliRunner()
  result = runner.invoke(status)
  assert result.exit_code == 0
  #assert result.output == 'Hello Peter!\n'
  #print(f"output: {result.output}")

def run(cmd):
  output = Popen(cmd, stdout=PIPE, 
    stderr=PIPE, shell=True, 
    cwd="/global/cscratch1/sd/jaws/jfroula", universal_newlines=True)

  stdout,stderr=output.communicate()
  rc=output.returncode
  return rc,stdout,stderr



if __name__ == "__main__":
  test_hello_world()
  rc,stdout,stderr = run('ls -l')
  print(f"stdout: {stdout}")
  print(f"stderr: {stderr}")
  print(f"errcode: {rc}")
