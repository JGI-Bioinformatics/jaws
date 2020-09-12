from setuptools import setup
import os

setup(
    name="jaws_worker",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Worker",
    url="https://gitlab.com/jgi-dsi/aa/jaws/worker",
    author="The JAWS Team",
    packages=["jaws_worker"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-worker = jaws_worker.cli:jaws", ]},
    zip_safe=False,
)
