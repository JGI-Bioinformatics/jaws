from setuptools import setup
import os

setup(
    name="jaws_run",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Run Service",
    url="https://gitlab.com/jgi-dsi/aa/jaws/run",
    author="The JAWS Team",
    packages=["jaws_run"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-run = jaws_run.cli:jaws", ]},
    zip_safe=False,
)
