import os

# from jaws_site import config

abspath = os.path.abspath(__file__)
dirpath = os.path.dirname(abspath)
conf = config.Configuration(os.path.join(dirpath, "..", "jaws-site.ini"))
