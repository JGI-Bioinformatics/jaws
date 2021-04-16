"""
Deprecated commands receive useful warnings.
"""

import click
import warnings


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def run():
    """JAWS Run Commands (DEPRECATED)"""
    pass


@run.command()
def queue():
    """DEPRECATED: try "jaws queue" instead"""
    warnings.warn(queue.__doc__)


@run.command()
def history():
    """DEPRECATED: try "jaws history" instead"""
    warnings.warn(history.__doc__)


@run.command()
def status():
    """DEPRECATED: try "jaws status" instead"""
    warnings.warn(status.__doc__)


@run.command()
def task_status():
    """DEPRECATED: try "jaws task-status" instead"""
    warnings.warn(task_status.__doc__)


@run.command()
def metadata():
    """DEPRECATED: try "jaws metadata" instead"""
    warnings.warn(metadata.__doc__)


@run.command()
def log():
    """DEPRECATED: try "jaws log" instead"""
    warnings.warn(log.__doc__)


@run.command()
def task_log():
    """DEPRECATED: try "jaws task-log" instead"""
    warnings.warn(task_log.__doc__)


@run.command()
def errors():
    """DEPRECATED: try "jaws errors" instead"""
    warnings.warn(errors.__doc__)


@run.command()
def cancel():
    """DEPRECATED: try "jaws cancel" instead"""
    warnings.warn(cancel.__doc__)


@run.command()
def list_sites():
    """DEPRECATED: try "jaws list_sites" instead"""
    warnings.warn(list_sites.__doc__)


@run.command()
def submit():
    """DEPRECATED: try "jaws submit" instead"""
    warnings.warn(submit.__doc__)


@run.command()
def inputs():
    """DEPRECATED: try "jaws inputs" instead"""
    warnings.warn(inputs.__doc__)


@run.command()
def validate():
    """DEPRECATED: try "jaws validate" instead"""
    warnings.warn(validate.__doc__)


@run.command()
def get():
    """DEPRECATED: try "jaws get" instead"""
    warnings.warn(get.__doc__)
