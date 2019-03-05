#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) Patrick Meyers (2019)
#
# This file is part of code for the seismic_radiometer package
#
# seismic_radiometer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# seismic_radiometer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SeisMon.  If not, see <http://www.gnu.org/licenses/>

"""Setup script for seismic_radiometer
"""

import glob
import os.path
from setuptools import (find_packages, setup)
from setuptools.command.install import install
from distutils.core import Extension
from distutils.command.build import build
import sys
import numpy


PACKAGENAME = 'seispy'


# set version information

DESCRIPTION = 'Seismic Radiometer'
LONG_DESCRIPTION = ''
AUTHOR = 'Patrick Meyers'
AUTHOR_EMAIL = 'patrick.meyers@ligo.org'
LICENSE = 'GPLv3'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages(exclude=['utils'])

# find all scripts
scripts = glob.glob('bin/*')

setup(name=PACKAGENAME,
      description=DESCRIPTION,
      scripts=scripts,
      packages=packagenames,
      ext_modules=[],
      requires=[],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      long_description=LONG_DESCRIPTION,
      zip_safe=False,
      test_suite='seispy.tests',
      use_2to3=True
      )
