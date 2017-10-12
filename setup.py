#!/usr/bin/env python

from setuptools import setup, find_packages

exec(open('resonance/version.py').read())

setup(
    name='resonance',
    version=__version__,
    author='Jason K. Moore',
    author_email='moorepants@gmail.com',
    url="https://github.com/moorepants/resonance/",
    description='Learning mechanical vibrations through computation.',
    long_description=open('README.rst').read(),
    keywords="engineering vibrations mechanical",
    license='CC-BY 4.0',
    packages=find_packages(),
    install_requires=['numpy>=1.13',
                      'matplotlib>=2.1',
                      'scipy>=0.19',
                      'pandas>=0.20'],
    extras_require={'notebooks': ['notebook', 'ipywidgets']},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics',
        ],
)
