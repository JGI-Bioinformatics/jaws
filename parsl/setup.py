from setuptools import setup
import os

setup(
    name="jaws_parsl",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service Parsl Backend",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_parsl"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-parsl = jaws_parsl.cli:jaws", ]},
    zip_safe=False,
)
