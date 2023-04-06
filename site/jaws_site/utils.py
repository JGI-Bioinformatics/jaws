import json
import os


class Singleton(type):
    """
    Metaclass for creating singleton objects.
    """

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    def _destructor(cls):
        del cls._instances[cls]


def write_json_file(outfile: str, contents: dict):
    with open(outfile, "w") as fh:
        fh.write(json.dumps(contents, sort_keys=True, indent=4))
    os.chmod(outfile, 0o0664)
