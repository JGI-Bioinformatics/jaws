"""
JAWS Site Pubsub
~~~~~~~~~~~~~~~~

A basic publish/subscribe library for JAWS site.
"""

__path__ = __import__("pkgutil").extend_path(__path__, __name__)

from .core import JawsPubsub

__all__ =  ["JawsPubsub"]
