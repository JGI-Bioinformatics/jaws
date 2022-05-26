#!/usr/bin/env python
"""
This is a script to validate the runtime sections of a WDL.
Some basic validation of the runtime values are done in wdl_parser_functions.py => checkValueSyntax()
when the runtime dictionary is created.
"""

import sys
import re
import warnings


def pretty_warning(message, category, filename, lineno, line=None):
    return f"WARNING: {message}\n"


warnings.formatwarning = pretty_warning


class WdlRuntimeError(Exception):
    pass


class WdlRuntimeMemoryError(WdlRuntimeError):
    pass


class WdlRuntimeTimeError(WdlRuntimeError):
    pass


class WdlRuntimeCpuError(WdlRuntimeError):
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
            if value[0].isalpha():
                # if the value starts with a letter, it is a variable name
                return None
            elif re.search(r".+:.+:.+", value):
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
            else:
                sys.exit(
                    'Error in in the runtime section of task %s. "time" should have the format: "hh:mm:ss". You had %s'
                    % (task, value)
                )

        # memory: "250G"
        if key == "mem":
            sys.exit(
                'Error in task %s. "mem" has been deprecated in favor of "memory"'
                % (task)
            )

        if key == "memory":
            if value[0].isalpha():
                # if the value starts with a letter, it is a variable name
                return None
            elif re.search(r"[gG]", value):
                value = value.upper()
                return value
            else:
                sys.exit(
                    'Error in the runtime section of task %s. "memory" should have '
                    'a "G" to specify gigabytes. You had %s' % (task, value)
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
            m = re.match(r"task\s+(\w+)", line)
            if m:
                task = m.groups()[0]
                stanza_content = []

            # s = stanza + r"\s*(?:{|<<<)"
            # match the stanza
            s = stanza
            p = re.compile(s)
            m = p.match(line)
            if m:
                # The stanza is found, so start keeping track of brackets.
                # If the stanza is "command", we can expect either {} or <<< >>> surrounding the block.
                # Therefore, this expression will match either case (?:{|<<<).
                start_count = 1
                left_bracket_count = 0
                right_bracket_count = 0

                # make sure there is a starting bracket
                if "{" not in line and "<<<" not in line:
                    raise WdlRuntimeError(
                        f'Error for task {task}. The stanza: "{line}" seems to be missing an opening bracket.'
                    )

            if start_count:
                stanza_content.append(line)
                left_bracket = re.findall("(?:{|<<<)", line)
                if left_bracket:
                    left_bracket_count += len(left_bracket)

                right_bracket = re.findall("(?:}|>>>)", line)
                if right_bracket:
                    right_bracket_count += len(right_bracket)

                # save results when we've reached end of stanza
                if right_bracket_count == left_bracket_count:

                    # remove instance of "<stanza> {" (i.e. "runtime {") and last "}" bracket.
                    stanza_str = (":::").join(stanza_content)

                    s = stanza + r"\s*(?:{|<<<)"
                    text_after = re.sub(
                        s, "", stanza_str
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
                    for param_key in stanza_list:
                        if param_key.startswith("#"):
                            continue

                        next_stanza = re.search(r"(?:{|<<<)", param_key)

                        # make sure the last bracket is not missing, or we'll bleed into the next stanza
                        if next_stanza:
                            raise WdlRuntimeError(
                                f'Error for task {task}. "{param_key}" seems to be the next stanza and not a '
                                "parameter key name. Maybe there was a missing bracket."
                            )

                        param_key_array = param_key.split(":", 1)
                        if len(param_key_array) < 2:
                            raise WdlRuntimeError(
                                f"Error for task {task}. This parameter {param_key} may have an "
                                "empty value or one that couldn't be parsed"
                            )
                        mk = param_key_array[0]
                        mv = param_key_array[1]
                        # cleanup junk from value
                        mv = mv.replace('"', "").strip()
                        mv = self.checkValueSyntax(task, mk, mv)
                        mydict[mk] = mv

                    self.stanza_dict[task] = mydict
                    start_count = 0


def spellCheck(task_name, task_dict):
    accepted_param_names = ["docker", "time", "memory", "cpu"]
    for key in task_dict.keys():
        if key not in accepted_param_names:
            warnings.warn(
                "%s is not a known runtime parameter and will be ignored." % key,
                SyntaxWarning,
            )


def allRequiredParams(task_name, task_dict):
    """Check that the user has included the minimum runtime parameters and that other params
    are acceptable values
    """
    if "cpu" not in task_dict:
        warnings.warn(
            "Task: %s allRequiredParams. %s is a required parameter for runtime"
            % (task_name, "cpu")
        )

    if "time" not in task_dict:
        warnings.warn(
            "Task: %s allRequiredParams. %s is a required parameter for runtime"
            % (task_name, "cpu")
        )

    if "memory" not in task_dict:
        raise WdlRuntimeMemoryError(
            "Task: %s allRequiredParams. %s is a required parameter for runtime"
            % (task_name, "memory")
        )


def memoryParam(task_name, task_dict, compute_max_ram_gb=None):
    """Check that the user hasn't requested too much memory for the specified resource."""

    if "memory" not in task_dict:
        raise WdlRuntimeMemoryError(
            "Task: %s memoryParam. %s is a required parameter for runtime."
            % (task_name, "memory")
        )
        return
    elif task_dict["memory"] is None:
        # memory is a variable name (i.e. memory shall be calculated); nothing to do
        return

    # memory is a string that include a "G" for gigabytes, so grab just the int
    mem = int(re.sub("[gG]", "", task_dict["memory"]))

    if compute_max_ram_gb and mem > compute_max_ram_gb:
        raise WdlRuntimeMemoryError(
            f"Task: {task_name} memoryParam. You are limited to {compute_max_ram_gb} for the requested site, you had {mem}"  # noqa
        )


def validate_wdl_runtime(wdl: str, compute_max_ram_gb: float = None) -> None:
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

        # Check the spelling of the runtime param keys. Give a warning if a paramemter is not in the list
        spellCheck(task_name, task_dict)

        # Check that we have minimum runtime params (i.e. time & memory)
        allRequiredParams(task_name, task_dict)

        # verify memory has been defined
        memoryParam(task_name, task_dict, compute_max_ram_gb)
