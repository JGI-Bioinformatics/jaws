from setuptools import setup
import os

setup(
    name="jaws_rpc",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS internal RPC library",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_rpc"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    zip_safe=False,
)
