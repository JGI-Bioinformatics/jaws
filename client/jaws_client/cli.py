"""
JAWS CLI
"""

import click
import os
import requests
import json
import logging
import globus_sdk
from jaws_client import log, config, analysis, catalog, user

JAWS_CLIENT_LOG = "JAWS_CLIENT_LOG"
JAWS_DEFAULT_LOG = "~/jaws.log"
JAWS_CLIENT_CONFIG = "JAWS_CLIENT_CONFIG"
JAWS_DEFAULT_CONFIG = "~/jaws.conf"


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", default=None, help="Config INI file")
@click.option("--log", "log_file", default=None, help="Log file")
def cli(config_file: str, log_file: str):
    """JGI Analysis Workflows Service."""
    if log_file is None:
        if JAWS_CLIENT_LOG in os.environ:
            log_file = os.environ[JAWS_CLIENT_LOG]
        else:
            log_file = os.path.expanduser(JAWS_DEFAULT_LOG)
    logger = log.setup_logger(__package__, log_file)

    if config_file is None:
        if JAWS_CLIENT_CONFIG in os.environ:
            config_file = os.environ[JAWS_CLIENT_CONFIG]
        else:
            config_file = os.path.expanduser(JAWS_DEFAULT_CONFIG)
    conf = config.Configuration(config_file)
    jaws = conf.get("JAWS", "name")
    logger.debug(f"Using {jaws} : {config_file}")


@cli.command()
def status() -> None:
    """Current system status."""
    url = f'{config.conf.get("JAWS", "url")}/status'
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException:
        raise SystemExit("JAWS Central is DOWN")
    if r.status_code != 200:
        raise SystemExit(r.text)
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
        raise SystemExit("Authentication failed")

    # POST TOKENS TO JAWS CENTRAL
    data = {}
    auth_service_name = "auth.globus.org"
    auth = token_response.by_resource_server[auth_service_name]
    data["auth_refresh_token"] = auth["refresh_token"]
    transfer_service_name = "transfer.api.globus.org"
    transfer = token_response.by_resource_server[transfer_service_name]
    data["transfer_refresh_token"] = transfer["refresh_token"]
    current_user = user.User()
    url = f'{config.conf.get("JAWS", "url")}/auth'
    try:
        r = requests.post(url, data=data, headers=current_user.header())
    except requests.exceptions.RequestException:
        raise SystemExit("Unable to communicate with JAWS server")
    if r.status_code != 200:
        raise SystemExit(r.text)


def jaws():
    """Entrypoint for jaws-client app."""
    cli.add_command(analysis.run)
    cli.add_command(catalog.wdl)
    cli()
