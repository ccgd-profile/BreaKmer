#! /usr/bin/python

from setuptools import setup

setup(name='BreaKmer',
      version='0.1',
      description=''
      url='https://github.com/a-bioinformatician/BreaKmer'
      author='Ryan Abo'
      author_email='ryan.abo@gmail.com'
      license='MIT'
      py_modules=['BreaKmer'],
      install_requires=[
        'pysam',
        'Bio'
      ]  
      )
