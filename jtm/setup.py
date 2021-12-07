from setuptools import setup
import os

setup(
    name="jaws_jtm",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service - JGI Task Manager",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_jtm", "jaws_jtm.lib"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "jtm=jaws_jtm.jtm:jtm",
        ]
    },
    scripts=['jaws_jtm/performance_metrics.py'],
    zip_safe=False,
)
