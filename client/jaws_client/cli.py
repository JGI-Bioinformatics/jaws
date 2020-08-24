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

JAWS_LOG_ENV = "JAWS_CLIENT_LOG"
JAWS_USER_LOG = os.path.expanduser("~/jaws.log")
JAWS_CONFIG_ENV = "JAWS_CLIENT_CONFIG"
JAWS_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-client.conf")
JAWS_USER_CONFIG_ENV = "JAWS_USER_CONFIG"
JAWS_USER_CONFIG_DEFAULT_FILE = os.path.expanduser("~/jaws-user.conf")


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "jaws_config_file", default=None, help="JAWS config file")
@click.option("--user", "user_config_file", default=None, help="User config file")
@click.option("--log", "log_file", default=None, help="Log file")
@click.option("--log-level", "log_level", default="INFO", help="Logging level")
def cli(jaws_config_file: str, user_config_file: str, log_file: str, log_level: str):
    """JGI Analysis Workflows Service"""
    if log_file is None:
        log_file = (
            os.environ[JAWS_LOG_ENV] if JAWS_LOG_ENV in os.environ else JAWS_USER_LOG
        )
    logger = log.setup_logger(__package__, log_file, log_level)
    if jaws_config_file is None:
        jaws_config_file = (
            os.environ[JAWS_CONFIG_ENV]
            if JAWS_CONFIG_ENV in os.environ
            else JAWS_CONFIG_DEFAULT_FILE
        )
    if user_config_file is None:
        user_config_file = (
            os.environ[JAWS_USER_CONFIG_ENV]
            if JAWS_USER_CONFIG_ENV in os.environ
            else JAWS_USER_CONFIG_DEFAULT_FILE
        )
    conf = config.Configuration(jaws_config_file, user_config_file)
    if conf:
        logger.debug(f"Config using {jaws_config_file}, {user_config_file}")


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
def info() -> None:
    """JAWS version and info."""
    url = f'{config.conf.get("JAWS", "url")}/info'
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
