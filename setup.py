#!/usr/bin/env python

from __future__ import division, print_function

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from hy454 import __version__ as _hy454_version

setup(name='hy454',
      version=_hy454_version,
      description='454 UDS analysis tools',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/veg/hy454',
      license='GNU GPL version 3',
      packages=['hy454'],
      package_dir={'hy454': 'lib/hy454'},
      package_data={'hy454': [
            'hyphy/*.bf',
            'data/fonts/ttf/*.ttf',
            'data/scores/*.dat'
    ]},
      data_files=[
            'scripts/codonaligner',
            'scripts/coveragegrapher',
            'scripts/seqlogo'
      ],
      requires=['Bio', 'BioExt', 'fakemp', 'freetype', 'hypy', 'matplotlib', 'numpy']
     )
