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

    A new parameter is introduced to the init() function, env_var_override, which is a prefix for
    any environment variables that should be used as overrides for the get() function.
    """

    def __init__(self,  env_var_override=None, **kwargs):
        if env_var_override is not None:
            self._vars = dict()  # Dictionary to be used as vars param in get(), populated from env vars
            # ToDo: take the value of env_var_override as a prefix for any environment variables meant
            # to be passed as a dictionary to the get() vars argument, which overrides any settings.
        else:
            self._vars = None
        configparser.BasicInterpolation.__init__(self, **kwargs)

    def before_get(self, parser, section, option, value, defaults):
        value = super().before_get(parser, section, option, value, defaults)
        return os.path.expandvars(value)

    def get(self, section, option, **kwargs):
        """
        Override the standard getter by forcing the vars argument to use the dictionary derived
        from the environment variables with the prefix in the env_var_override param to init
        """
        value = super().get(self, section, option, vars=self._vars, **kwargs)
        return(value)
