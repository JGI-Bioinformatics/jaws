#! /usr/bin/env python
#
# JtmInterface returns
#   0 if ready
#   1 if queued
#   2 if running
#   4 if successfully done
#   -1, -2, -3, -4: failed
#
# jtm-isalive exits with 0 if it's in ['ready', 'queued', 'running'] status
#                       1 if it's done successfully or failed
#
import sys

from jaws_jtm.lib.common import eprint
from jaws_jtm.lib.jtminterface import JtmInterface


def isalive():
    assert len(sys.argv) == 2, "USAGE: jtm-isalive <task_id>"
    ret = int(
        JtmInterface("status", info_tag=sys.argv[1]).call(task_id=int(sys.argv[1]))
    )
    if ret == -88:
        eprint("jtm-isalive: command timeout.")
        sys.exit(-1)
    elif ret in [0, 1, 2]:
        print("yes")
        sys.exit(0)
    else:
        print("no")
        sys.exit(1)
