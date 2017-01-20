#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_spectrum',
                           sources=['spectrum_wrap.c', 'spectrum.c','likelihood.c'],
                           )

setup (name = 'spectrum',
       version = '??',
       author      = "SWIG Docs",
       description = """Spectrum from LL distribution""",
       ext_modules = [example_module],
       py_modules = ["spectrum"],
       )
