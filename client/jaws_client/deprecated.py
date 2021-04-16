"""Deprecated commands print useful error message"""


import warnings


import click
@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def run():
    """JAWS Run-Workflows Commands (DEPRECATED)"""
    pass


@run.command()
def queue() -> None:
    """"DEPRECATED: please try 'jaws queue' or 'jaws --help' instead"""
    warnings.warn(queue.__doc__)


@run.command()
def history():
    """"DEPRECATED: please try 'jaws history' or 'jaws --help' instead"""
    warnings.warn(history.__doc__)


@run.command()
def status():
    """"DEPRECATED: please try 'jaws status' or 'jaws --help' instead"""
    warnings.warn(status.__doc__)


@run.command()
def task_status():
    """"DEPRECATED: please try 'jaws task-status' or 'jaws --help' instead"""
    warnings.warn(task_status.__doc__)


@run.command()
def metadata():
    """"DEPRECATED: please try 'jaws metadata' or 'jaws --help' instead"""
    warnings.warn(metadata.__doc__)


@run.command()
def log():
    """"DEPRECATED: please try 'jaws log' or 'jaws --help' instead"""
    warnings.warn(log.__doc__)


@run.command()
def task_log():
    """"DEPRECATED: please try 'jaws task-log' or 'jaws --help' instead"""
    warnings.warn(task_log.__doc__)


@run.command()
def errors():
    """"DEPRECATED: please try 'jaws errors' or 'jaws --help' instead"""
    warnings.warn(errors.__doc__)


@run.command()
def cancel():
    """"DEPRECATED: please try 'jaws cancel' or 'jaws --help' instead"""
    warnings.warn(cancel.__doc__)


@run.command()
def list_sites():
    """"DEPRECATED: please try 'jaws list-sites' or 'jaws --help' instead"""
    warnings.warn(list_sites.__doc__)


@run.command()
def submit():
    """"DEPRECATED: please try 'jaws submit' or 'jaws --help' instead"""
    warnings.warn(submit.__doc__)


@run.command()
def inputs():
    """"DEPRECATED: please try 'jaws inputs' or 'jaws --help' instead"""
    warnings.warn(inputs.__doc__)


@run.command()
def validate(wdl_file: str) -> None:
    """"DEPRECATED: please try 'jaws validate' or 'jaws --help' instead"""
    warnings.warn(validate.__doc__)


@run.command()
def get(run_id: int, dest: str) -> None:
    """"DEPRECATED: please try 'jaws get' or 'jaws --help' instead"""
    warnings.warn(get.__doc__)
