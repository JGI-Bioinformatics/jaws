"""
JAWS Workflows Repository
"""

import sys
import json
import requests
from bs4 import BeautifulSoup
import click
from . import config, user


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def wdl() -> None:
    """JAWS Workflows Catalog"""
    pass


@wdl.command()
def list() -> None:
    """List available workflows in the Catalog.

    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
def versions(name: str) -> None:
    """List available versions of specified workflow.

    :param name: The name of the workflow
    :type name: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
@click.argument("version")
def about(name: str, version: str) -> None:
    """Return README document for a workflow.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}/doc'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
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
@click.argument("name")
@click.argument("version")
def get(name: str, version: str) -> None:
    """Get workflow specification (WDL) for a workflow.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}'
    try:
        r = requests.get(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.text
    print(result)


@wdl.command()
@click.argument("name")
@click.argument("version")
@click.argument("wdl_file")
@click.argument("md_file")
def add(name: str, version: str, wdl_file: str, md_file: str) -> None:
    """Add a new workflow to the catalog.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :param wdl_file: Path to the workflow specification (WDL) file
    :type wdl_file: str
    :param md_file: Path to the README file in markdown format
    :type md_file: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}'
    files = {"wdl_file": open(wdl_file, "r"), "md_file": open(md_file, "r")}
    try:
        r = requests.post(url, files=files, headers=current_user.header())
    except Exception:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
@click.argument("version")
def release(name: str, version: str) -> None:
    """ Mark a version as released, which makes it's WDL immutable.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :return:
    """
    current_user = user.User()
    data = {"release": True}
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}'
    try:
        r = requests.put(url, data=data, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
@click.argument("version")
@click.argument("wdl_file")
def update_wdl(name: str, version: str, wdl_file: str) -> None:
    """Update a workflow's WDL in the catalog.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :param wdl_file: Path to the workflow specification (WDL) file
    :type wdl_file: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}'
    files = {"wdl_file": open(wdl_file, "r")}
    try:
        r = requests.put(url, files=files, headers=current_user.header())
    except Exception:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
@click.argument("version")
@click.argument("md_file")
def update_doc(name: str, version: str, md_file: str) -> None:
    """Update a workflow's README in the catalog.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :param md_file: Path to the README file in markdown format
    :type md_file: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}/doc'
    files = {"md_file": open(md_file, "r")}
    try:
        r = requests.put(url, files=files, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@wdl.command()
@click.argument("name")
@click.argument("version")
def delete(name: str, version: str) -> None:
    """Remove a workflow from the catalog.

    :param name: The name of the workflow
    :type name: str
    :param version: The version of the workflow
    :type version: str
    :return:
    """
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/wdl/{name}/{version}'
    try:
        r = requests.delete(url, headers=current_user.header())
    except requests.ConnectionError:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))
