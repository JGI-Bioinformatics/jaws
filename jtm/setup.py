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
    packages=["jaws_jtm"],
    install_requires=[line.strip() for line in open("requirements.txt")],
    entry_points={
        "console_scripts": [
            "jtm=jaws_jtm.jtm:jtm",
            "jtm-check-manager=jaws_jtm.jtm_check_manager:check_manager",
            "jtm-check-worker=jaws_jtm.jtm_check_worker:check_worker",
            "jtm-isalive=jaws_jtm.jtm_isalive:isalive",
            "jtm-kill=jaws_jtm.jtm_kill:kill",
            "jtm-manager=jaws_jtm.jtm_manager:manager",
            "jtm-resource-log=jaws_jtm.jtm_resource_log:resource_log",
            "jtm-remove-pool=jaws_jtm.jtm_remove_pool:remove_pool",
            "jtm-status=jaws_jtm.jtm_status:status",
            "jtm-submit=jaws_jtm.jtm_submit:submit",
            "jtm-worker=jaws_jtm.jtm_worker:worker",
        ]
    },
    zip_safe=False,
)
