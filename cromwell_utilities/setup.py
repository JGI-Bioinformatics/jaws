from setuptools import setup
import os

setup(
    name="cromwell_utilities",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="Cromwell utilities",
    url="https://gitlab.com/jgi-dsi/aa/jaws/cromwell_utilities",
    author="The JAWS Team",
    packages=["cromwell_utilities"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["cromwell-utils = cromwell_utilities.cromwell_utils:cromwell", ]},
    zip_safe=False,
)
