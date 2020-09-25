from setuptools import setup
import os

setup(
    name="jaws_auth",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS OAuth2 Server",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_auth"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-auth = jaws_auth.cli:jaws", ]},
    zip_safe=False,
)
