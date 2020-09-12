from setuptools import setup
import os

setup(
    name="jaws_grid",
    version=os.popen("git describe --dirty=-dev --always --tags --abbrev=6")
    .read()
    .strip()
    .replace("-", "+", 1),
    description="JAWS Grid",
    url="https://gitlab.com/jgi-dsi/aa/jaws/grid",
    author="The JAWS Team",
    packages=["jaws_grid"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    # entry_points={"console_scripts": ["jaws-grid = jaws_grid.cli:jaws", ]},  # none yet
    zip_safe=False,
)
