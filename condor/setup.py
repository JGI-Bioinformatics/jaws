from setuptools import setup
import os

setup(
    name="jaws_condor",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service Condor Backend Pool Manager",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_condor"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "condor_pool_add=jaws_condor.add:cli",
            "condor_pool_remove=jaws_condor.remove:cli",
            "jaws_condor=jaws_condor.cli:jaws"
        ]
    },
    zip_safe=False,
)
