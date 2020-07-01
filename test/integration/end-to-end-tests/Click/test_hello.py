#!/usr/bin/env python

from click.testing import CliRunner
from hello import hello

def test_hello_world():
  runner = CliRunner()
  result = runner.invoke(hello, ['Peter'])
  assert result.exit_code == 0
  assert result.output == 'Hello Peter!\n'
  print(f"result: {result}")

if __name__ == "__main__":
  test_hello_world()
