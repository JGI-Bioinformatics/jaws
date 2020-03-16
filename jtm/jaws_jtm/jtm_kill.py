#! /usr/bin/env python
#
# JtmInterface returns code
#   0: terminated successfully
#  -1: no task found
#   4: the task has already been completed
#
# jtm-kill exits with 0 if it's successfully terminated.#
#                     1 else
#
import sys
from jaws_jtm.lib.jtminterface import JtmInterface
from jaws_jtm.lib.run import eprint


def kill():
    assert len(sys.argv) == 2, "USAGE: jtm-kill <task_id>"
    taskID = int(sys.argv[1])
    ret = int(JtmInterface('kill',
                           info_tag=taskID).call(task_id=taskID))
    if ret == -88:
        eprint("jtm-kill: command timeout.")
        sys.exit(-1)
    elif ret == -5:
        eprint("jtm-kill: task id not found.")
        sys.exit(-1)
    sys.exit(0) if ret == 0 else sys.exit(1)
    # if ret != 0:
    #     print "jtm-kill failed with task id %d" % taskID
    # sys.exit(0)
