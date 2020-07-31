from setuptools import setup
import os

setup(name='jaws_client',
      version=os.popen('git describe --dirty=-dev --always --tags --abbrev=6').read().strip().replace("-", "+", 1),
      description='JGI Analysis Workflow Service Client',
      url='https://gitlab.com/jgi-dsi/aa/jaws/client',
      author='Edward Kirton',
      packages=['jaws_client'],
      install_requires=[line.strip() for line in open("requirements.txt")],
      entry_points={
        'console_scripts': [
            'jaws = jaws_client.cli:jaws',
        ]
      },
      zip_safe=False)
