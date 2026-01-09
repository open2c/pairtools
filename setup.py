#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import re
import glob


from setuptools import find_packages, setup
from setuptools.extension import Extension

try:
    from Cython.Distutils import build_ext as _build_ext
    from Cython.Build import cythonize
except ImportError:
    raise ImportError('Cython is now required to build the extension modules.')


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("pairtools", "__init__.py"),
        re.MULTILINE,
    ).group(1)
    return version


def get_ext_modules():
    ext = ".pyx"
    src_files = glob.glob(
        #os.path.join(os.path.dirname(__file__), "pairtools", "lib", "*" + ext)
        os.path.join("pairtools", "lib", "*" + ext)
    )

    ext_modules = []
    for src_file in src_files:
        name = "pairtools.lib." + os.path.splitext(os.path.basename(src_file))[0]
  
        if 'pysam' in name:
            import pysam
            ext_modules.append(
                Extension(
                    name,
                    [src_file],
                    extra_link_args=pysam.get_libraries(),
                    include_dirs=pysam.get_include(),
                    define_macros=pysam.get_defines(),
                )
            )
        elif "regions" in name:
            ext_modules.append(
                Extension(
                    name,
                    [src_file],
                    language="c++",
                )
            )

        else:
            ext_modules.append(Extension(name, [src_file]))

    ext_modules = cythonize(ext_modules)  # , annotate=True

    return ext_modules


class build_ext(_build_ext):
    # Extension module build configuration
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Fix to work with bootstrapped numpy installation
        # http://stackoverflow.com/a/21621689/579416
        # Prevent numpy from thinking it is still in its setup process:
        #__builtins__.__NUMPY_SETUP__ = False
        import numpy

        self.include_dirs.append(numpy.get_include())

    def run(self):

        # Import numpy here, only when headers are needed
        import numpy

        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())

        # Call original build_ext command
        _build_ext.run(self)


setup(
    version=get_version(),
    ext_modules=get_ext_modules(),
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    # entry_points={
    #     "console_scripts": [
    #         "pairtools = pairtools.cli:cli",
    #     ]
    # },
    packages=find_packages(),
)
