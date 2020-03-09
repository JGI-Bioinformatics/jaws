from setuptools import setup
import os

setup(
    name="jaws_auth",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service OAuth2 Server",
    url="https://gitlab.com/jgi-dsi/aa/jaws/auth",
    author="The JAWS Team",
    packages=["jaws_auth"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={"console_scripts": ["jaws-auth = jaws_auth.server:auth", ]},
    zip_safe=False
)
