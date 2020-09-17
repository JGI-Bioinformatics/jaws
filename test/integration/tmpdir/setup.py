from setuptools import setup

setup(name='test',
      version='0.1.0',
      packages=['test_it'],
      entry_points={
          'console_scripts': [
              'tt = test_it.tt:main'
          ]
      },
      )
