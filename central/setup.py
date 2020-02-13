from setuptools import setup
import os

setup(name='jaws_central',
      version=os.popen('git describe --dirty=-dev --always --tags --abbrev=6').read().strip().replace("-", "+", 1),
      description='JGI Analysis Workflow Service Central Server',
      url='https://gitlab.com/jgi-dsi/aa/jaws/site',
      author='The JAWS Team',
      packages=['jaws_central'],
      install_requires=[ line.strip() for line in open("requirements.txt") ], #TODO: this needs jtm requirements
      entry_points={
        'console_scripts': [
            'jaws-central = jaws_central.cli:jaws',
        ]
      },
      zip_safe=False)
