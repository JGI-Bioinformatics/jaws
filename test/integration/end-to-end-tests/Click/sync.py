#!/usr/bin/env python

import click

@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
   click.echo('Debug mode is %s' % ('on' if debug else 'off'))

@cli.command()
@click.argument('name')
def sync(name):
   click.echo(name)
