from jaws_client import config


class User:
    """User class"""
    token_file = None
    access_token = None

    def __init__(self) -> None:
        """Load user's access token from file"""
        self.access_token = config.conf.get("USER", "token")
        if not self.access_token:
            raise SystemExit("User access token required; an contact admin to get yours.")

    def header(self) -> str:
        """Return HTTP OAuth2 header containing the authentication token"""
        header = {"Authorization": f"Bearer {self.access_token}"}
        return header
