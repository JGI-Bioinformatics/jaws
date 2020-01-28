from setuptools import setup
import os

setup(name='jaws_site',
      version=os.popen('git describe --dirty=-dev --always --tags --abbrev=6').read().strip().replace("-", "+", 1),
      description='JGI Analysis Workflow Service Site RPC Server',
      url='https://gitlab.com/jgi-dsi/aa/jaws/site',
      author='The JAWS Team',
      packages=['jaws_site'],
      install_requires=[ line.strip() for line in open("requirements.txt") ], #TODO: this needs jtm requirements
      entry_points={
        'console_scripts': [
            'jaws-site = jaws_site.cli:jaws',
            'jaws-site-jtm = jaws_site.jtm.cli:jtm', #TODO: this should not be a separate entrypoint, but a subcommand of jaws-site
        ]
      },
      zip_safe=False)
