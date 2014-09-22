#! /usr/bin/python

from setuptools import setup

def readme() :
  with open('README.md') as f:
    return f.read()

setup(name='BreaKmer',
      version='0.0.7',
      description='Structural variation detection tool, designed for targeted sequencing data.',
      long_description=readme(),
      url='https://github.com/a-bioinformatician/BreaKmer',
      author='Ryan Abo',
      author_email='ryan.abo@gmail.com',
      license='LICENSE',
      py_modules=['BreaKmer'],
      install_requires=[
        'pysam >= 0.6',
        'biopython >= 1.62'
      ]  
      )
