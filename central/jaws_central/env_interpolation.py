import configparser
import os


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

    """

    def __init__(self,  env_override=None, **kwargs):
        if env_override is not None:
            # Dictionary to be used as vars param in get(), populated from env vars by looking
            # for the prefix env_override in the name, and then strip off the prefix for the
            # name of the setting.
            strlen = len(env_override)
            if strlen < 3 or strlen > 20:  # The prefix must be 3-20 characters long
                raise(ValueError("env_override prefix must be from 3-20 characters in length"))
            # should change to use removeprefix() for Python 3.9+ instead of slicing off prefix
            self._vars = {k[strlen:]: os.environ[k] for k in os.environ.keys() if k.startswith(env_override)}
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
            # any calls to get() with an explcitly set vars
            kwargs['vars'] = self._vars
        value = super().get(section, option, **kwargs)
        return(value)
