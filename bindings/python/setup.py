# Author(s): Maurice Laveaux
# Copyright: see the accompanying file COPYING or copy at
# https://github.com/mCRL2org/mCRL2/blob/master/COPYING
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#

import os
import sys
import glob

from setuptools import setup

MCRL2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

sys.path += [os.path.join(MCRL2_ROOT, '3rd-party', 'pybind11')]
from pybind11.setup_helpers import Pybind11Extension, ParallelCompile

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension(
        "mcrl2py",
        sorted([
            "mcrl2py.cpp"
            ] + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "atermpp", "source", "*.cpp"))
              + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "process", "source", "*.cpp"))
              + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "data", "source", "*.cpp"))
              + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "utilities", "source", "*.cpp"))
              + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "lps", "source", "*.cpp"))
              + glob.glob(os.path.join(MCRL2_ROOT, "libraries", "core", "source", "*.cpp"))
        ),
        include_dirs=[
            os.path.join(MCRL2_ROOT, "3rd-party", "pybind11", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "atermpp", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "process", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "data", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "utilities", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "lps", "include"),
            os.path.join(MCRL2_ROOT, "libraries", "core", "include"),
        ],
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
        include_pybind11=False,
        cxx_std=17,
    ),
]

# Optional multithreaded build
ParallelCompile("NPY_NUM_BUILD_JOBS").install()

setup(
    name="mcrl2py",
    version=__version__,
    author="Maurice Laveaux",
    author_email="m.laveaux@tue.nl",
    url="https://github.com/mCRL2org/mCRL2",
    description="Python bindings for mCRL2 using pybind",
    license_files = (os.path.join(MCRL2_ROOT, "COPYING")),
    long_description="",

    ext_modules=ext_modules,
    package_data={"mcrl2py": ["mcrl2py.pyi"]},
    extras_require={},
    zip_safe=True,
    python_requires=">=3.7"
)