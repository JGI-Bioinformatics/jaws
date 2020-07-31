#!/usr/bin/env python

from click.testing import CliRunner
from sync import cli

def test_sync():
  runner = CliRunner()
  result = runner.invoke(cli, ['--debug', 'sync', 'we are Syncing'])
  assert result.exit_code == 0
  assert 'Debug mode is on' in result.output
  assert 'Syncing' in result.output

  print(f"result: {result} ")

if __name__ == "__main__":
    test_sync()
