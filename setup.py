#!/usr/bin/env python 

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='lbcred',
      version='0.1',
      author='Kirsten Casey',
      packages=['lbcred'],
      url='https://github.com/kirstencasey/LBCreduce')
