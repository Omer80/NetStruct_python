# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:38:43 2016

@author: omertz@post.bgu.ac.il
"""
__version__= 0.1
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
from setuptools import setup

setup(name='netstruct_python',
      version='0.1',
      description='Genetic population structure analysis using networks',
      url='http://github.com/storborg/funniest',
      author='Gili Greenbaum and Omer Tzuk',
      author_email='omertz@post.bgu.ac.il',
      license='MIT',
      packages=['netstruct_python'],
      install_requires=[
          'numpy','scipy','igraph','deepdish',
      ],
      zip_safe=False)