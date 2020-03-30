import os
import sys
import logging
from jaws_client import config


class User:
    """User class"""
    token_file = None
    access_token = None

    def __init__(self) -> None:
        """Load user's access token from file"""
        logger = logging.getLogger(__package__)
        if "HOME" not in os.environ:
            sys.exit('Env var "HOME" not defined')
        self.token_file = os.path.join(os.environ["HOME"], f'.{config.conf.get("JAWS", "name")}')
        logger.debug(f"Reading user access token from {self.token_file}")
        if not os.path.isfile(self.token_file):
            sys.exit(f"Access token file not found: {self.token_file}.  Get yours from a JAWS Admin.")
        with open(self.token_file, "r") as token_file:
            self.access_token = token_file.read().strip()
        if not self.access_token:
            sys.exit(f"Access token file not found: {self.token_file}.  Get yours from a JAWS Admin.")

    def header(self) -> str:
        """Return HTTP OAuth2 header containing the authentication token"""
        header = {"Authorization": f"Bearer {self.access_token}"}
        return header
