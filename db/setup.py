from setuptools import setup
import os

setup(
    name="jaws_db",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Database Class",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_db"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    zip_safe=False,
)
