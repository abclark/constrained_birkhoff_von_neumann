"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

#: Installation requirements.
requirements = ['numpy', 'networkx>=2.0']

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='generalized_birkhoff_von_neumann',  # Required
    version='0.0.1.dev1',  # Required
    description='A generalized Birkoff von Neumann decomposition',  # Required
    long_description=__doc__,  # Optional
    url='https://github.com/abclark/generalized_birkhoff_von_neumann',  # Optional
    author='Aubrey Clark',  # Optional
    author_email='aubs.bc@gmail.com',  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required
)
