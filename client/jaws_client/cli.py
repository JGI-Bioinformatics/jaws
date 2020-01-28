"""
Top-level Click 'group' for JAWS command-line client
"""

import click
from jaws_client import analysis,catalog,utils

@click.group()
def cli():
    """
    JGI Analysis Workflows Service (JAWS)
    """
    pass

def jaws():
    cli.add_command(analysis.run)
    cli.add_command(catalog.wdl)
    cli.add_command(utils.util)
    cli()
