from setuptools import setup
import os

setup(
    name="jaws_catalog",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Catalog",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_catalog"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-catalog = jaws_catalog.cli:jaws", ]},
    zip_safe=False,
)
