#!/usr/bin/env python

from setuptools import setup

setup(name='hy454',
      version='0.1.0',
      description='454 UDS analysis tools',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/veg/hy454',
      license='GNU GPL version 3',
      packages=['hy454'],
      package_dir={'hy454': 'src/hy454'},
      package_data={'hy454': ['hyphy/*.bf']},
      data_files=[('/usr/local/bin', ['src/codonaligner', 'src/graphcoverage'])],
      requires=['Bio', 'fakemp', 'matplotlib']
     )
