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
    packages=["jaws_parsl",
              "jaws_parsl.parsl_configs"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "jaws-parsl-recv=jaws_parsl.recv:cli",
            "jaws-parsl-send=jaws_parsl.send:send"
        ]
    },
    zip_safe=False,
)
