import configparser
import os
import re


class EnvInterpolation(configparser.BasicInterpolation):
    """
    Interpolation which expands environment variables in values and also sets the vars param
    of get() so that it pulls from a dictionary of environment variables to override any
    settings that have a key in the dictionary/environment.

    Cribbed shamelessly from https://gist.github.com/malexer/ee2f93b1973120925e8beb3f36b184b8
    Uses normal interpolation syntax for configparser
    https://docs.python.org/3/library/configparser.html#interpolation-of-values
    """

    def before_get(self, parser, section, option, value, defaults):
        value = super().before_get(parser, section, option, value, defaults)
        return os.path.expandvars(value)


class JAWSConfigParser(configparser.ConfigParser):
    """
    Subclass of the normal configparser that allows config options to be totally overriden by an
    environment variable. ignoring any other settings in the config file.

    A new parameter is introduced to the init() function, env_override, which is a prefix for
    any environment variables that should be used as overrides for the get() function. Must be
    at least 3 characters long and no longer than 20 characters long. This is basically to avoid
    prefixes that are too short colliding with normal environment variables, and also to avoid
    prefixes that are too long messing up config names that are longish. Any environment variables
    that match the prefix have the prefix removed from the name, and the name/value are added to
    the _vars dict. The environment is only checked when the config object is instantiated, no
    changes to env vars after instantiation are picked up.

    So, if env_override="ENV__" is passed to the constructor, then any environment variables like
    ENV__host and ENV__port are set, then _vars['host'] and _vars['port'] will be set to the values
    found, and entries in the config file (irrespective of section!) with the name "host" will be
    set to the value originally from ENV__host, and port will be from ENV__port. Whatever value is
    on the right hand side of the assignment in the config file are totally ignored.

    Note that env_override should be set as a positional parameter, not as a named parameter
    because the singleton metaclass constructor for jaws_central.config.Configuration will throw
    an error on a named parameter.

    After some review, we realized that section names will be necessary. If the section that comes
    after the env_override prefix starts with all caps and under underscores, then everything up
    the last underscore will be treated as a section name, and a nested dictionary will be created
    where the first level is the section name, and the second level is a dictionary with that sections
    settings. Section names remain uppercased, the trailing underscore is stripped.

    It is not acceptable to have some environment variables with sections and others without,
    this will cause an exception to be thrown.

    For example, if env_override is "_ENV_", and the environment contains the following variables:
    _ENV_username=myself
    _ENV_password=mypassword
    _ENV_url=http://jaws.jgi.gov:5000/

    The _vars will contain:
    { "username" : "myself",
      "password" : "mypassword",
      "url" : "http://jaws.jgi.gov:5000/"
    }

    On the other hand, if the environment contains these variables:
    _ENV_DB_username=myself
    _ENV_DB_password=mypassword
    _ENV_DB_host=mysql:9000
    _ENV_API_username=someone
    _ENV_API_password=somewpassword
    _ENV_API_url=http://jaws.jgi.gov:5000/

    Then _vars will contain:
    { "DB" : { "username" : "myself",
               "password" : "mypassword",
               "host" : "mysql:9000"
             },
      "API": { "username" : "someone",
               "password" : "somepassword",
               "url" : "http://jaws.jgi.gov:5000/"
             }
    }


    Finally, if the environment contains this:
    _ENV_username=myself
    _ENV_DB_password=mypassword
    _ENV_DB_host=mysql:9000
    _ENV_API_username=someone
    _ENV_API_password=somewpassword
    _ENV_API_url=http://jaws.jgi.gov:5000/

    Then a KeyError will be raised because all but the first environment variable starting with _ENV_ have section
    names, but the first does not.

    When the get() method is called, if _vars is a nest dictionary, then the section passed in will
    be used to look for a matching section in _vars, and the matching dict() will be passed to the
    parent method. No sections match then None will be passed for the vars parameter. If _vars is
    not a nest dictionary, then the contents of _vars will be passed in the vars param as is.
    """

    def __init__(self, env_override=None, **kwargs):
        if env_override is not None:
            # Dictionary to be used as vars param in get(), populated from env vars by looking
            # for the prefix env_override in the name, and then strip off the prefix for the
            # name of the setting.
            strlen = len(env_override)
            if strlen < 3 or strlen > 20:  # The prefix must be 3-20 characters long
                raise ValueError(
                    "env_override prefix must be from 3-20 characters in length"
                )
            # should change to use removeprefix() for Python 3.9+ instead of slicing off prefix
            basevars = {
                k[strlen:]: os.environ[k]
                for k in os.environ.keys()
                if k.startswith(env_override)
            }
            if len(basevars) > 0:
                section_pre = re.compile(
                    "^([A-Z][:0-9A-Z_]*)_(.+)"
                )  # match all caps section name, ending in _
                temp = dict()
                for key, val in basevars.items():
                    res = section_pre.match(key)
                    if res is not None:
                        sec = res.group(1)
                        opt = res.group(2)
                        if sec not in temp:
                            temp[sec] = {opt: val}
                        else:
                            temp[sec][opt] = val
                templen = sum(len(v) for v in temp.values())
                if templen == len(basevars):
                    self._vars = temp  # All environment variables options contained section names
                elif templen == 0:
                    self._vars = basevars  # No environment variables options contained section names
                else:
                    raise KeyError(
                        "Cannot mix options with and without section prefixes"
                    )
            else:
                self._vars = dict()
        else:
            self._vars = None
        configparser.ConfigParser.__init__(self, **kwargs)

    def get(self, section, option, **kwargs):
        """
        Override the standard getter by forcing the vars argument to use the dictionary derived
        from the environment variables with the prefix in the env_override param. Call the parent
        getter with _vars set to use the existing functionality for overriding values.
        """
        if self._vars is not None:
            # Force vars to self._vars if it has been set. This means env_override will also override
            # any calls to get() with an explicitly set vars

            # Uppercase the section name so that it matches any section names in _vars, should preclude
            # accidental matches against option name
            usection = section.upper()
            if usection in self._vars:
                kwargs["vars"] = self._vars[usection]
            elif len(self._vars) > 0:
                # Section didn't match anything, then either we don't have a matching section
                # and we use an empty dict() for vars
                # or else _vars is a simple 1 level dictionary and we just pass along all of _vars
                key1 = next(iter(self._vars))
                if type(self._vars[key1]) is dict:
                    kwargs["vars"] = dict()
                else:
                    kwargs["vars"] = self._vars
        value = super().get(section, option, **kwargs)
        return value
