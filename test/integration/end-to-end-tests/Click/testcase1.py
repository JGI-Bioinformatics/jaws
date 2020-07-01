#!/usr/bin/env python

import click

@click.command()
def status():
   click.echo('running status')

