import logging
from types import new_class
from typing import Dict

from jaws_site import config
from jaws_site.utils import run_sh_command

logger = logging.getLogger(__package__)


class PoolManager:
    """Class representing a single Run"""

    def __init__(self, **kwargs):
        self.site = config.conf.get("SITE", "id"),

    def add_workers():
        logger.info("Checking to add workers to pool")
        return None

    def rm_workers():
        logger.info("Checking to remove workers to pool")
        return None
