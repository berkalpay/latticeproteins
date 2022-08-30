from setuptools import setup
from distutils.extension import Extension

contactlooper = Extension('latticeproteins.contactlooper', sources = ['latticeproteins/contactlooper.c'])

setup(
    name = 'latticeproteins',
    version = '0.0',
    description = 'Code for lattice protein simulations.',
    url='https://github.com/berkalpay/latticeproteins',
    packages = ['latticeproteins'],
    ext_modules = [contactlooper],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
