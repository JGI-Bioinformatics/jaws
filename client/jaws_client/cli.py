"""
JAWS CLI
"""

import sys
import click
import os
import requests
import json
import logging
import globus_sdk
from jaws_client import log, config, analysis, catalog, user


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """JGI Analysis Workflows Service"""
    pass


def find_config_file(env_var: str, config_file: str) -> str:
    """Find config file by checking env var, cwd, and home, in that order.

    :param config_file: filename of configuration file
    :type config_file: str
    :param env_var: name of config file environment variable
    :type env_var: str
    :return: absolute path if found, None otherwise.
    :rtype: str
    """
    if env_var in os.environ:
        cf = os.environ[env_var]
        if os.path.isfile(cf):
            return cf
        else:
            raise IOError(f"File not found: {cf}")

    cf = os.path.abspath(config_file)
    if os.path.isfile(cf):
        return cf

    cf = os.path.join(os.environ["HOME"], config_file)
    if os.path.isfile(cf):
        return cf
    raise IOError("Unable to find config file")


@cli.command()
def status() -> None:
    """Current system status."""
    url = f'{config.conf.get("JAWS", "url")}/status'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        sys.exit("JAWS Central is DOWN")
    if r.status_code != 200:
        sys.exit(r.text)
    result = r.json()
    print(json.dumps(result, indent=4, sort_keys=True))


@cli.command()
def login() -> None:
    """Authenticate using Globus OAuth2 and provide tokens to JAWS-Central"""
    logger = logging.getLogger(__package__)
    requested_scopes = (
        "openid profile email urn:globus:auth:scope:transfer.api.globus.org:all"
    )
    globus_client = globus_sdk.NativeAppAuthClient(
        config.conf.get("GLOBUS", "client_id")
    )
    globus_client.oauth2_start_flow(
        requested_scopes=requested_scopes, refresh_tokens=True
    )
    authorize_url = globus_client.oauth2_get_authorize_url()
    print(f"Please go to this URL and log in:\n{authorize_url}")
    auth_code = input("Paste the authorization code here: ").strip()
    logger.debug("Authenticate with Globus")
    try:
        token_response = globus_client.oauth2_exchange_code_for_tokens(auth_code)
    except globus_sdk.GlobusError:
        sys.exit("Authentication failed")

    # POST TOKENS TO JAWS CENTRAL
    data = {}
    auth_service_name = config.conf.get("GLOBUS", "auth_service_name")
    auth = token_response.by_resource_server[auth_service_name]
    data["auth_refresh_token"] = auth["refresh_token"]
    transfer_service_name = config.conf.get("GLOBUS", "transfer_service_name")
    transfer = token_response.by_resource_server[transfer_service_name]
    data["transfer_refresh_token"] = transfer["refresh_token"]
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/auth'
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except requests.exceptions.RequestException:
        sys.exit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        sys.exit(r.text)


def jaws():
    """Entrypoint for jaws-client app."""
    log_file = os.path.join(os.environ["HOME"], "jaws.log")
    logger = log.setup_logger(__package__, log_file)

    config_file = find_config_file("JAWS_CLIENT_CONFIG", "jaws-client.ini")
    if not config_file:
        logger.critical("Unable to find jaws config file")
    config.conf = config.JawsConfig(config_file)

    cli.add_command(analysis.run)
    cli.add_command(catalog.wdl)
    cli()
