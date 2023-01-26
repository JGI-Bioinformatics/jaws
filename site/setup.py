from setuptools import setup
import os

setup(
    name="jaws_site",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service Site Server",
    url="https://gitlab.com/jgi-dsi/aa/jaws/site",
    author="The JAWS Team",
    packages=["jaws_site"],
    package_data={'': ['datatransfer_plugins/*']},
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "jaws-site = jaws_site.cli:jaws",
            "pagurus = jaws_site.pagurus:main",
        ]
    },
    zip_safe=False,
)
