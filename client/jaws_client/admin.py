"""
JAWS administrator functions
"""

import json
import requests
import click
import logging

from jaws_client import config, user


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def admin():
    """JAWS Administrator Commands"""
    pass


@admin.command()
@click.argument("uid")
@click.argument("email")
@click.argument("name")
@click.option("--admin", default=False, help="Grant admin access")
def log(uid: str, email: str, name: str, admin: bool) -> None:
    """Add new user and get JAWS OAuth access token.

    :param uid: New user's ID
    :type id: str
    :param email: New user's email address
    :type id: str
    :param name: New user's full name
    :type id: str
    :param admin: Grant new user administrator access
    :type id: bool
    :return: JSON response contains new user's access token
    :rtype: dict
    """
    logger = logging.getLogger(__package__)
    current_user = user.User()
    data = {
        "uid": uid,
        "email": email,
        "name": name,
        "admin": admin
    }
    url = f'{config.conf.get("JAWS", "url")}/user'
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except Exception as error:
        logger.error(f"Unable to add user: {error}")
        raise SystemExit(f"Unable to add user: {error}")
    result = r.json()
    if r.status_code != 201:
        raise SystemExit(result["detail"])
    logger.info(f"Add new user: {uid} ({email})")
    print(json.dumps(result, indent=4, sort_keys=True))
