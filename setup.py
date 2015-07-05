#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

setup(name="libmihmodule",
    ext_modules=[
    #modules=[
        Extension("libmihmodule", ["src/array32.cpp", "src/mihasher.cpp", "src/sparse_hashtable.cpp",  "src/bucket_group.cpp",  "src/reorder.cpp", "src/mihmodule.cpp",], library_dirs=['/'], libraries = ["boost_python", "hdf5"], include_dirs=['include', '/opt/local/include'])
    ])
