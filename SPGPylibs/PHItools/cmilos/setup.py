# python setup.py build_ext --inplace
# python setup.py build_ext --inplace --use-cython

import platform
from distutils.core import setup
from distutils.extension import Extension
# from Cython.Build import cythonize
import numpy
import sys

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.c'

compile_extra_args = []
link_extra_args = [] 

if platform.system() == "Windows":
    compile_extra_args = ["/std:c++latest", "/EHsc"]
elif platform.system() == "Darwin":
    compile_extra_args = ["-mmacosx-version-min=10.13"]
    link_extra_args = ["-mmacosx-version-min=10.13"]

extensions = [Extension(
    name="pymilos",
    sources=["pymilos.pyx"],
    language='c',
    libraries=["milos"],
    library_dirs=["lib"],
    include_dirs=["lib",numpy.get_include()],
    extra_compile_args = compile_extra_args,
    extra_link_args = link_extra_args)]
    # extra_compile_args = ['-Wunknown-warning-option '],                               #addon
    # extra_compile_args = ['-fopenmp'],                               #addon
    # define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')] #addon

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name="pymilos",
    ext_modules = extensions
)
