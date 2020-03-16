#!/usr/bin/env python
import click

from jaws_jtm.lib.config import COMPUTE_RESOURCES


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug):
    print("Debug mode is %s" % ("on" if debug else "off"))


@cli.command()
@click.option(
    "--task-json-file", "-f", help="user task in json format",
)
@click.option(
    "--cluster",
    "-cl",
    help="cluster (site) name to run task",
    type=click.Choice(COMPUTE_RESOURCES, case_sensitive=False),
)
def submit(task_json_file, cluster):
    print("submit")
    print(task_json_file)
    print(cluster)


@cli.command()
def kill():
    print("kill")


@cli.command()
def isalive():
    print("isalive")


@cli.command()
def status():
    print("status")


# if __name__ == '__main__':
def jtm():
    cli()
