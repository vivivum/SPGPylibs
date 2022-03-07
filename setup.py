# python setup.py build_ext --inplace
# python setup.py build_ext --inplace --use-cython
# python setup.py install build_ext --inplace 

from setuptools import setup, find_packages
#--pymilos
import sys
import platform
# from distutils.core import setup, find_packages
from distutils.extension import Extension
# from Cython.Build import cythonize
import subprocess
import numpy
#--pymilos

from distutils.command.install import install as _install
class install(_install):
    def run(self):
        subprocess.call(['make','-C', 'SPGPylibs/PHItools/cmilos/','default'])
        _install.run(self)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Parse version number from SPGPylibs/__init__.py:
with open('SPGPylibs/__init__.py') as f:
    info = {}
    for line in f:
        if line.startswith('version'):
            exec(line, info)
            break

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
    sources=["SPGPylibs/PHItools/cmilos/pymilos.pyx"],
    language='c',
    libraries=["milos"],
    library_dirs=["SPGPylibs/PHItools/cmilos/lib"],
    include_dirs=["SPGPylibs/PHItools/cmilos/lib",numpy.get_include()],
    extra_compile_args = compile_extra_args,
    extra_link_args = link_extra_args)]
    # extra_compile_args = ['-Wunknown-warning-option '],                               #addon
    # extra_compile_args = ['-fopenmp'],                               #addon
    # define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')] #addon

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup_info = dict(
    name="SPGPylibs",
    version=info['version'],
    author="SPG",
    author_email="orozco@iaa.es",
    description="Solar Physics Group python tools for PHI and TuMAG instruments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vivivum/SPGPylibs",
    project_urls={
        "Bug Tracker": "https://github.com/vivivum/SPGPylibs/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "./"},
    packages=find_packages(where="./"),
    install_requires=['numpy','scipy','photutils','matplotlib','astropy','tqdm','cython','spiceypy','scikit_learn','pyFFTW'],
    python_requires=">=3.6",
    ext_modules = extensions,
    cmdclass = {'install': install}
) 

print(info['version'])
print(find_packages(where="./"))
setup(**setup_info)