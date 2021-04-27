#!/usr/bin/env python
"""
This is a script to validate the runtime sections of a WDL.
Some basic validation of the runtime values are done in wdl_parser_functions.py => checkValueSyntax()
when the runtime dictionary is created.
"""

import sys
import re


class WdlRuntimeError(Exception):
    pass


class WdlRuntimeMemoryError(WdlRuntimeError):
    pass


class WdlRuntimeTimeError(WdlRuntimeError):
    pass


class WDLStanzas:
    """
    This class contains functions to parse and manipulate the
    different stanzas of a WDL file.

    Attributes
    ----------
    f:  filehandler
        filehandler of the wdl file
    stanza_dict: dictionary
        this is the main dictionary to which all stanzas will be added.
    """

    def __init__(self, wdl):
        """This just opens a WDL file and saves the filehandler in a variable

        Input: string
             the path to the wdl
        Output: filehandler
             saves a filehandler (f) of the wdl
        """

        self.wdl = wdl

    def checkValueSyntax(self, task, key, value):
        """
        make sure the values for the runtime params are a string or int and have the correct formatting.

        time: "00:30:00"
        memory: "250G"
        poolname: "some-unique-name"
        shared: 0
        node: 1
        nwpn: 1
        constraint: "skylake"
        qos: "jgi_shared"
        account: "fungalp"
        cpu: 12
        cluster: cori
        """

        # time: "00:30:00"
        if key == "time":
            time = value.split(":")
            if len(time) != 3:
                sys.exit(
                    'Error in task %s. Time is not in the correct format of "hh:mm:ss"'
                    % task
                )
            z = 0
            for i in time:
                i = int(i)
                z = z + i
            if z == 0:
                sys.exit(
                    "Error in task %s. No time was specified (everything was zero). You had %s"
                    % (task, value)
                )

        # memory: "250G"
        if key == "mem":
            sys.exit(
                'Error in task %s. "mem" has been deprecated in favor of "memory"'
                % (task)
            )

        if key == "memory":
            if re.search(r"[gG]", value):
                value = value.upper()
                return value
            else:
                sys.exit(
                    'Error in task %s. "memory" should have a "G" to specify gigabytes. You had %s'
                    % (task, value)
                )

            if value < 1:
                sys.exit(
                    'Error in task %s. "memory" should have a value greater than 1. You had %s'
                    % (task, value)
                )

        # shared: 0
        if key == "shared":
            if int(value) != 0 and int(value) != 1:
                sys.exit(
                    'Error in task %s. "shared" has to be an int, 1 or 0. You had %s'
                    % (task, value)
                )

        # if value is string, then make sure it's lower case
        # otherwise it is suppose to be an int
        if re.search(r"[A-Za-z:]", value):
            value = value.lower()
        else:
            value = int(value)

        return value

    def loadStanza(self, stanza):
        """Load a stanza type into a dictionary

        Input: string
                the name of a stanza that you want to parse.
        Output: dictionary
                a dictionary of the stanzas that were found, where key is task name and values are in a list.
        Strategy: go through each line of the wdl document and when the specified stanza block is found
        start keeping track of '{' and '}' brackets.  When the total number of '}' closing brackets
        matches '{' opening brackets, then add all the contents (i.e. all lines from
        when the stanza was seen) to a dictionary.
        """

        accepted_task_stanza_keys = [ "input", "command", "meta", "runtime", "output", "parameter_meta", ] # noqa
        accepted_keys = [ "mem", "cpu", "time", "poolname", "poolsize", "cluster", "nwpn", ] # noqa
        start_count = 0
        left_bracket_count = 0
        right_bracket_count = 0
        task = ""  # the name of the task
        self.stanza_dict = (
            {}
        )  # this is the main dictionary to which all stanzas will be added.

        # look at each line of the WDL
        for line in self.wdl.split("\n"):
            line = line.strip()

            # parse the task name
            if m := re.match(r"task\s+(\w+)", line):
                task = m.groups()[0]
                stanza_content = []

            # The stanza is found, so start keeping track of brackets.
            # If the stanza is "command", we can expect either {} or <<< >>> surrounding the block.
            # Therefore, this expression will match either case (?:{|<<<).
            s = stanza + r"\s*((?:{|<<<))"
            p = re.compile(s)
            if m := p.match(line):
                start_count = 1

            if start_count:
                stanza_content.append(line)
                if left_bracket := re.findall("(?:{|<<<)", line):
                    left_bracket_count += len(left_bracket)

                if right_bracket := re.findall("(?:}|>>>)", line):
                    right_bracket_count += len(right_bracket)

                # save results when we've reached end of stanza
                if right_bracket_count == left_bracket_count:

                    # remove instance of "<stanza> {" (i.e. "runtime {") and last "}" bracket.
                    stanza_str = (":::").join(stanza_content)

                    # bracket =
                    # if stanza == command and bracket == '{':

                    # replace_me = stanza + "\s*" + bracket
                    # new_set = {x.replace(replace_me, '') for x in stanza_content}
                    text_after = p.sub(
                        "", stanza_str
                    )  # remove unessessary stanza name (i.e. runtime {)
                    text_after = re.sub(
                        r"(?:}|>>>)", "", text_after
                    )  # remove right-most bracket '}' or '>>>'
                    text_after = re.sub(
                        r"^\s*", "", text_after
                    )  # remove left-most space
                    text_after = re.sub(
                        r"\s*$", "", text_after
                    )  # remove right-most space

                    # create string back into a list
                    stanza_list = text_after.split(":::")
                    stanza_list = [
                        x for x in stanza_list if x
                    ]  # remove all empty elements in list created by regex

                    # make stanza_list into a dictionary (so we'll have a dictionary of a dictionary)
                    mydict = {}
                    for pair in stanza_list:
                        pair_arry = pair.split(":", 1)
                        if len(pair_arry) < 2:
                            sys.exit(
                                'Error for task %s. This parameter "%s" may have an empty value (or one that couldn\'t be parsed).' # noqa: E501,E261
                                % (task, pair)
                            )
                        mk = pair_arry[0]
                        mv = pair_arry[1]
                        # cleanup junk from value
                        mv = mv.replace('"', "").strip()
                        mv = self.checkValueSyntax(task, mk, mv)
                        mydict[mk] = mv

                    self.stanza_dict[task] = mydict
                    start_count = 0


def allRequiredParams(task_name, task_dict):
    """Check that the user has included the minimum runtime parameters and that other params
    are acceptable values
    """
    if "time" not in task_dict:
        raise WdlRuntimeTimeError(
            "Task: %s Test: allRequiredParams. %s is a required parameter for runtime"
            % (task_name, "time")
        )

    if "memory" not in task_dict:
        raise WdlRuntimeMemoryError(
            "Task: %s Test: allRequiredParams. %s is a required parameter for runtime"
            % (task_name, "memory")
        )

    # check that qos is an acceptable value.
    accepted_qos = ["jgi_shared", "jgi_exvivo", "regular"]
    if "qos" in task_dict and task_dict["qos"] not in accepted_qos:
        raise WdlRuntimeError(
            "Task: %s Test: memoryParam. The value for qos is not a recognized value for skylake. You had: %s."
            % (task_name, task_dict["qos"])
        )

    # check that constraint is an acceptable value.
    accepted_constraint = [
        "haswell",
        "knl",
        "skylake",
        "lr3_c32,jgi_m256",
        "lr3_c32,jgi_m512",
    ]
    if (
        "constraint" in task_dict
        and task_dict["constraint"].lower() not in accepted_constraint
    ):
        raise WdlRuntimeError(
            'Task: %s Test: timeParam. "constraint" must be one of the following values: %s. We found "%s"'
            % (task_name, accepted_constraint, task_dict["constraint"].lower())
        )


def timeParam(task_name, task_dict):
    """
    Check that if constraint: is set to skylake or knl then the max time is respected.
    If constraint: is anything else, then check that max time is 72hrs.
    Also check that constraint is not set to some non-recognized value.
    """
    if "time" not in task_dict:
        raise WdlRuntimeError(
            'Task: %s Test: timeParam. %s is a required parameter for runtime so "timeParam" test has been skipped.'
            % (task_name, "time")
        )
        return

    hours = task_dict["time"].split(":")[0]

    if "constraint" in task_dict:
        myconstraint = task_dict["constraint"].lower()

        # check skylake mem
        if myconstraint == "skylake":
            if int(hours) > 168:
                raise WdlRuntimeError(
                    "Task: %s Test: timeParam. You are limited to 168hrs where constraint=%s"
                    % (task_name, myconstraint)
                )
        # check knl mem
        elif myconstraint == "knl":
            if int(hours) > 48:
                raise WdlRuntimeError(
                    "Task: %s Test: timeParam. You are limited to 48hrs where constraint=%s"
                    % (task_name, myconstraint)
                )
        # check haswell or other mem
        else:
            # if constraint exists but is not skylake or knl, then it is haswell or jgi?, so limit 72hrs.
            if int(hours) > 72:
                raise WdlRuntimeError(
                    "Task: %s Test: timeParam. You are limited to 72hrs where constraint=%s"
                    % (task_name, myconstraint)
                )
    else:
        # if constraint is not included in the runtime, the default is haswell, so 72hrs limit
        if int(hours) > 72:
            raise WdlRuntimeError(
                'Task: %s Test: timeParam. You are limited to 72hrs when constraint is the default value("%s")'
                % (task_name, "haswell")
            )


def memoryParam(task_name, task_dict):
    """
    check the user hasn't requested too much memory for the specified resource
    memory: '5G'  # You have get '115G' for all machines j|'500G']

    For JGI
     lr3                              316        64     32      64      72
     lr3        lr3_c32,jgi_m256       32       256     32      64      72
     lr3        lr3_c32,jgi_m512        8       512     32      64      72
     jgi                               40       256     32      64      72

    For skylake
    memory: '250G' if qos: 'jgi_shared'
    memory: '758G' if qos: 'jgi_exvivo'

    For knl
    memory: 96

    For haswell
    memory: 128
    """

    # TODO: make checks for jgi

    # we've already verified that memory: is a string that include a "G" for gigabytes when the task_dict
    # was created, so grab just the int.
    if "memory" not in task_dict:
        raise WdlRuntimeMemoryError(
            'Task: %s Test: memoryParam. %s is a required parameter for runtime. The "memoryParam" test has been skipped.' # noqa: E501,E261
            % (task_name, "memory")
        )
        return

    mem = int(re.sub("[gG]", "", task_dict["memory"]))

    # if no constraint, then default is "haswell" which has 128G max
    if "constraint" not in task_dict:
        if mem > 128:
            raise WdlRuntimeMemoryError(
                "Task: %s Test: memoryParam. You are limited to 128G for haswell, which is the default constraint when non is specified. You had %s" # noqa: E501,E261
                % (task_name, task_dict["memory"])
            )

    # skylake
    # TODO: can skylake have other values? What about qos in general?
    if "constraint" in task_dict and task_dict["constraint"] == "skylake":
        if "account" not in task_dict:
            raise WdlRuntimeMemoryError(
                "Task: %s Test: memoryParam. If you are using skylake, you must have account: set to fungalp." # noqa: E501,E261
                % (task_name)
            )
        if "qos" not in task_dict:
            raise WdlRuntimeMemoryError(
                "Task: %s Test: memoryParam. If you are using skylake, you must have qos: set to jgi_exvivo(250G) or jgi_shared(758G)" # noqa: E501,E261
                % (task_name)
            )

        if "qos" in task_dict and re.search(r"[gG]", task_dict["memory"]):
            if task_dict["qos"] not in ["jgi_shared", "jgi_exvivo"]:
                raise WdlRuntimeMemoryError(
                    "Task: %s Test: memoryParam. If you are using skylake, you must have qos: set to jgi_exvivo(250G) or jgi_shared(758G): Your value was %s" # noqa: E501,E261
                    % (task_name, task_dict["qos"])
                )

            mem = re.sub("[gG]", "", task_dict["memory"])
            if task_dict["qos"] == "jgi_shared" and int(mem) > 250:
                raise WdlRuntimeMemoryError(
                    "Task: %s Test: memoryParam. You are limited to 250G for qos: jgi_shared. You had %s"
                    % (task_name, task_dict["memory"])
                )
            elif task_dict["qos"] == "jgi_exvivo" and int(mem) > 758:
                raise WdlRuntimeMemoryError(
                    "Task: %s Test: memoryParam. You are limited to 758G for qos: jgi_exvivo. You had %s"
                    % (task_name, task_dict["memory"])
                )

    # knl
    if "constraint" in task_dict and task_dict["constraint"] == "knl":
        if mem > 96:
            raise WdlRuntimeMemoryError(
                "Task: %s Test: memoryParam. You are limited to 96G for knl. You had %s"
                % (task_name, task_dict["memory"])
            )


def runtimeCombinations(task_name, task_dict):
    """
    # this is handled in the memoryParam function
      memory: "250G"
      constraint: "skylake"
      qos: "jgi_shared"
      account: "fungalp"

    # this is handled in the memoryParam function
      memory: "758G"
      constraint: "skylake"
      qos: "jgi_exvivo"
      account: "fungalp"

    # this is handled below
      Using non-priority queue ("genepool")
      memory: "10G"
      qos: "regular"
      account: "m342"
    """

    # Using non-priority queue ("genepool")
    if "qos" in task_dict and task_dict["qos"] == "regular":
        if "account" not in task_dict:
            raise WdlRuntimeError(
                "Task %s Test: runtimeCombinations. 'account' is required when qos: 'regular'."
                % (task_name)
            )
        if task_dict["account"] != "m342":
            raise WdlRuntimeError(
                "Task %s Test: runtimeCombinations. 'account' needs to be set to 'm342' when qos: 'regular'."
                % (task_name)
            )

    # the skylake combinations have been checked in the memoryParam function.


def validate_wdl_runtime(wdl: str) -> None:
    # Run the validations
    #
    # A dictionary of all the runtime stanzas is created here.
    # The keys are the task names and the value is a dictionary of each parameter like:
    # assy': {'docker: container','time: 03:00:00'}
    # The parenthesis are removed.
    doc = WDLStanzas(wdl)
    doc.loadStanza("runtime")

    for task_name in doc.stanza_dict.keys():
        task_dict = doc.stanza_dict[task_name]

        # Check that we have minimum runtime params (i.e. time & memory)
        allRequiredParams(task_name, task_dict)

        # Check that if constraint: is set to skylake or knl then the max time is respected.
        # If constraint: is anything else, then check that max time is 72hrs.
        # Also check that constraint is not set to some non-recognized value.
        timeParam(task_name, task_dict)

        # check the user hasn't requested too much memory for the specified resource
        memoryParam(task_name, task_dict)

        # Some runtime params require other params to be set. Check that the combinations
        # of runtime parameters are correct.
        runtimeCombinations(task_name, task_dict)
