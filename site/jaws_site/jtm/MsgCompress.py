#! /usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Copyright 2018-2019 (C) DOE JGI, LBL
# Seung-Jin Sul (ssul@lbl.gov)

"""
compress and pickle msg
"""

import sys
if sys.version_info[0] < 3:
    import cPickle
else:  # py3
    # import _pickle as cPickle
    import pickle as cPickle
import zlib


#-------------------------------------------------------------------------------
def zdumps(obj):
    """
    Dumps pickleable object into zlib compressed string
    :param obj: pickle it and compress
    :return:
    """
    ## 1 is fastest and produces the least compression,
    ## 9 is slowest and produces the most.
    ## 0 is no compression
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 5)
    # Todo: match the protocol # among manager, worker, and the JTM CLI
    #  ref) https://stackoverflow.com/questions/38466523/python-how-to-handle-the-valueerror-unsupported-pickle-protocol-4-error
    #  Py2 ==> max=2, Py3 ==> max=4
    # return zlib.compress(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL), 5)  # py3


#-------------------------------------------------------------------------------
def zloads(zstr):
    """
    Loads pickleable object from zlib compressed string
    :param zstr: compress msg
    :return:
    """
    return cPickle.loads(zlib.decompress(zstr))
    # return pickle.loads(zlib.decompress(zstr))  # py3


## To do
## json or msgpack can be used for
##
##  pickle -- If you have no desire to support any language other than
    #Python, then using the pickle encoding will gain you the support of all
    #built-in Python data types (except class instances), smaller messages when
    #sending binary files, and a slight speedup over JSON processing. .
#
#import json
#
#msg = {
#        "version": "1.0",
#        "host": "www1",
#        "short_message": "Short message",
#        "full_message": "Backtrace here and more stuff",
#}
#zmessage = zlib.compress(str(msg))
#data_json = json.dumps(msg)
#
#connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
#channel = connection.channel()
#
#print str(msg)
#channel.queue_declare(queue='nodeLogs')
#channel.basic_publish(exchange='', routing_key='nodeLogs', body=data_json)
#print " [x] Sent 'Hello World!'"
#connection.close()