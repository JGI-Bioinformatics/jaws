#!/usr/bin/env python

"""
Miscellaneous utilities
"""

import sys
import os
from flask import make_response, abort, Flask, request, redirect, url_for, current_app


def status():
    """
    Check system health
    """
    #config = current_app.config
    # TODO: check all systems, not just this one
    result = { "JAWS": "UP" }
    return result, 200

