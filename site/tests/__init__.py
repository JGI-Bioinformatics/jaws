import os

# from jaws_site import config

abspath = os.path.abspath(__file__)
dirpath = os.path.dirname(abspath)
conf = jaws_site.config.Configuration(os.path.join(dirpath, "..", "jaws-site.ini"))
