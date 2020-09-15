"""
Workflows Catalog REST endpoints.
"""

import os
import click
import markdown
from typing import Tuple
from jaws_central import catalog


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def wdl() -> None:
    """JAWS Workflows Catalog"""
    pass


@wdl.command()
def list_wdls(user: str) -> Tuple[dict, int]:
    """Retrieve workflows from database.

    :param user: Current user's ID
    :type user: str
    :return:
    """
    result = catalog.list_wdls()
    print(result)


@wdl.command()
@click.argument("user")
@click.argument("name")
def get_versions(user: str, name: str) -> Tuple[dict, int]:
    """Returns a list of all (non-deprecated) versions of a particular workflow.

    :param user: Current user's ID
    :type user: str
    :param name: Name of the workflow
    :param str:
    :return:
    """
    result = catalog.get_versions(name)
    print(result)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
def get_doc(user: str, name: str, version: str) -> Tuple[str, int]:
    """Returns a workflow's README (stored in md format) as HTML.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    doc = catalog.get_doc(name, version)
    print(markdown.markdown(doc))


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
def get_wdl(user: str, name: str, version: str) -> Tuple[str, int]:
    """Returns a workflow's WDL

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    wdl = catalog.get_wdl(name, version)
    print(wdl)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
def release_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Tag a workflow as "released", which makes it's WDL immutable.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    catalog.release_wdl(user, name, version)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
def del_wdl(user: str, name: str, version: str) -> Tuple[dict, int]:
    """Delete a workflow.  If it was "released", then tags as "deprecated", rather than being purged from db.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :return:
    """
    catalog.del_wdl(user, name, version)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
@click.argument("wdl_file")
def update_wdl(user: str, name: str, version: str, wdl_file: str) -> Tuple[dict, int]:
    """Update the WDL file for a workflow.  This cannot be updated after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param wdl_file: path to new WDL document
    :type wdl_file: str
    :return:
    """
    new_wdl = wdl_file.read()
    os.remove(wdl_file)
    catalog.update_wdl(user, name, version, new_wdl)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
@click.argument("doc_file")
def update_doc(user: str, name: str, version: str, doc_file: str) -> Tuple[dict, int]:
    """Update the doc file for a workflow.  This can be updated even after release.

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param doc_file: Path to new README document
    :type doc_file: str
    :return:
    """
    new_doc = doc_file.read()
    os.remove(doc_file)
    catalog.update_doc(user, name, version, new_doc)


@wdl.command()
@click.argument("user")
@click.argument("name")
@click.argument("version")
@click.argument("wdl_file")
@click.argument("doc_file")
def add_wdl(user: str, name: str, version: str, wdl_file: str, doc_file: str) -> Tuple[dict, int]:
    """Add a new workflow to the catalog

    :param user: Current user's ID
    :type user: str
    :param name: Name of a workflow
    :type name: str
    :param version: Version of a workflow
    :type version: str
    :param wdl_file: path to new WDL document
    :type wdl_file: str
    :param wdl_file: path to new WDL document
    :type wdl_file: str
    :return:
    """
    new_wdl = wdl_file.read()
    os.remove(wdl_file)
    new_doc = doc_file.read()
    os.remove(doc_file)
    catalog.add_wdl(user, name, version, new_wdl, new_doc)
