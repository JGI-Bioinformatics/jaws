from setuptools import setup, find_packages
import os

setup(
    name="jaws_central",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service Central REST Server",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    package_dir={'': 'jaws_central'},
    packages=find_packages('jaws_central'),
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-central = jaws_central.cli:jaws", ]},
    zip_safe=False,
)
