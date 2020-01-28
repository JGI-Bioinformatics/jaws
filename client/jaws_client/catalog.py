"""
JAWS Workflows Repository
"""

import sys
import os
import json
import getpass
import requests
#import html2text
from bs4 import BeautifulSoup
import csv
import time
import click
from jaws_client import user

# JAWS CONFIG
JAWS_URL = os.environ["JAWS_URL"]
jaws_catalog = "%s/wdl" % (JAWS_URL,)

@click.group()
def wdl():
    """
    Workflows Catalog
    """
    pass

@wdl.command()
def list():
    """
    List shared workflows
    """
    current_user = user.User()
    try:
        r = requests.get(jaws_catalog, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
def versions(name):
    """
    List available versions of a workflow
    """
    current_user = user.User()
    try:
        url = "%s/%s" % (jaws_catalog, name)
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
@click.argument('version')
def about(name, version):
    """
    Return README document for a workflow.
    """
    url = "%s/%s/%s/doc" % (jaws_catalog, name, version)
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.text
    result = result.rstrip('"\n')
    result = result.lstrip('"')
    result = result.replace("\\n", "\n")
    soup = BeautifulSoup(result, features="html.parser")
    print(soup.get_text())

@wdl.command()
@click.argument('name')
@click.argument('version')
def get(name, version):
    """
    Get WDL specification for a workflow.
    """
    url = "%s/%s/%s" % (jaws_catalog, name, version)
    current_user = user.User()
    try:
        r = requests.get(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.text
    print(result)

@wdl.command()
@click.argument('name')
@click.argument('version')
@click.argument('wdl_file')
@click.argument('md_file')
def add(name, version, wdl_file, md_file):
    """
    Add a workflow to the catalog
    """
    url = "%s/%s/%s" % (jaws_catalog, name, version)
    files = {
         'wdl_file': open(wdl_file,'r'),
         'md_file': open(md_file,'r')
    }
    current_user = user.User()
    try:
        r = requests.post(url, files=files, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
@click.argument('version')
def release(name, version):
    """
    Mark a version as immutable production release.
    """
    data = { "release" : True }
    url = "%s/%s/%s" % (jaws_catalog, name, version)
    current_user = user.User()
    try:
        r = requests.put(url, data=data, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
@click.argument('version')
@click.argument('wdl_file')
def update_wdl(name, version, wdl_file):
    """
    Update a workflow WDL in the catalog
    """
    url = "%s/%s/%s" % (jaws_catalog, name, version)
    files = {
         'wdl_file': open(wdl_file,'r')
    }
    current_user = user.User()
    try:
        r = requests.put(url, files=files, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
@click.argument('version')
@click.argument('md_file')
def update_doc(name, version, md_file):
    """
    Update a workflow README in the catalog
    """
    url = "%s/%s/%s/doc" % (jaws_catalog, name, version)
    files = {
         'md_file': open(md_file,'r')
    }
    current_user = user.User()
    try:
        r = requests.put(url, files=files, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

@wdl.command()
@click.argument('name')
@click.argument('version')
def delete(name, version):
    """
    Remove a workflow from the catalog.
    """
    url = "%s/%s/%s" % (jaws_catalog, name, version)
    current_user = user.User()
    try:
        r = requests.delete(url, headers=current_user.header())
    except:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))

#@wdl.command()
#@click.argument('name')
#@click.argument('version')
#def owners(name, version):
#    """
#    Show a workflow's owners
#    """
#    url = "%s/%s/%s/owners" % (jaws_catalog, name, version)
#    current_user = user.User()
#    try:
#        r = requests.get(url, headers=current_user.header())
#    except:
#        sys.exit("Unable to communicate with JAWS server")
#    if r.status_code != 200:
#        sys.exit(r.text)
#    result = r.json()
#    print(json.dumps(result, indent=4, sort_keys=True))

#@wdl.command()
#@click.argument('name')
#@click.argument('version')
#@click.argument('username')
#def add_owner(name, version, username):
#    """
#    Add another user to a workflow's owners' list.
#    """
#    data = { "new_user" : username }
#    url = "%s/%s/%s/owners" % (jaws_catalog, name, version)
#    current_user = user.User()
#    try:
#        r = requests.put(url, data=data, headers=current_user.header())
#    except:
#        sys.exit("Unable to communicate with JAWS server")
#    if r.status_code != 200:
#        sys.exit(r.text)
#    result = r.json()
#    print(json.dumps(result, indent=4, sort_keys=True))
