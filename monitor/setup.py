from setuptools import setup
import os

setup(
    name="jaws_monitor",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JGI Analysis Workflow Service Prometheus Monitor Server",
    url="https://code.jgi.doe.gov/advanced-analysis/jaws",
    author="The JAWS Team",
    packages=["jaws_prometheus"],
    include_package_data=True,
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "jaws-prometheus = jaws_prometheus.jaws_prometheus:main",
            "site-monitor = jaws_prometheus.site_monitor:main",
        ]
    },
    zip_safe=False,
)
