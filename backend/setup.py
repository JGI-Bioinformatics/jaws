from setuptools import setup
import os

setup(
    name="jaws_backend",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Cromwell Backend",
    url="https://gitlab.com/jgi-dsi/aa/jaws/backend",
    author="The JAWS Team",
    packages=["jaws_backend"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-backend = jaws_backend.cli:jaws", ]},
    zip_safe=False,
)
