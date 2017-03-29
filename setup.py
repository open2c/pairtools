#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import re

import numpy

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

classifiers = """\
    Development Status :: 4 - Beta
    Operating System :: OS Independent
    Programming Language :: Python
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.3
    Programming Language :: Python :: 3.4
    Programming Language :: Python :: 3.5
"""


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text

def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read('pairsamtools', '__init__.py'),
        re.MULTILINE).group(1)
    return version


def get_long_description():
    return _read('README.md')


install_requires = [
    'numpy>=1.10', 
    'cython>=0.25', 
    'click>=6.6', 
]


extensions = [
    Extension(
        "pairsamtools._dedup", ["pairsamtools/_dedup.pyx"],
    ),
]

packages = find_packages()
setup(
    name='pairsamtools',
    author='Mirny Lab',
    author_email='espresso@mit.edu',
    version=get_version(),
    license='BSD3',
    description='CLI tools to process mapped Hi-C data',
    long_description=get_long_description(),
    keywords=['genomics', 'bioinformatics', 'Hi-C', 'contact'],
    url='https://github.com/mirnylab/pairsamtools',
    packages=find_packages(),
    ext_modules = cythonize(extensions),
    zip_safe=False,
    classifiers=[s.strip() for s in classifiers.split('\n') if s],
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
             'pairsamtools = pairsamtools:cli',
        ]
    },
    include_dirs=[numpy.get_include()]
)
