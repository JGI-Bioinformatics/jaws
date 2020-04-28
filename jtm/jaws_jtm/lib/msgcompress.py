#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Seung-Jin Sul (ssul@lbl.gov)
"""
compress and pickle msg
"""
import pickle as cPickle
import zlib


# -------------------------------------------------------------------------------
def zdumps(obj):
    """
    Dumps pickleable object into zlib compressed string
    :param obj: pickle it and compress
    :return:
    """
    # 1 is fastest and produces the least compression,
    # 9 is slowest and produces the most.
    # 0 is no compression
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 5)  # py3


# -------------------------------------------------------------------------------
def zloads(zstr):
    """
    Loads pickleable object from zlib compressed string
    :param zstr: compress msg
    :return:
    """
    return cPickle.loads(zlib.decompress(zstr))  # py3
