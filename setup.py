#!/usr/bin/env python

# from distutils.core import setup, Extension
from setuptools import setup, Extension
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(os.path.dirname(__file__), 'genice_bondtwist', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))
    
long_desc = "".join(open("README.md").readlines())


setup(
    name='genice_bondtwist',
    version=metadata['version'],
    description='Bondtwist analysis plugin for GenIce.',
    long_description=long_desc,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
    ],
    author='Masakazu Matsumoto',
    author_email='vitroid@gmail.com',
    url='https://github.com/vitroid/genice-bondtwist/',
    keywords=['genice', 'chirality'],

    packages=['genice_bondtwist',
              'genice_bondtwist.formats',
    ],
    
    entry_points = {
        'genice_format_hook0': [
            'bondtwist  = genice_bondtwist.formats.bondtwist:hook0',
            'bondtwist2 = genice_bondtwist.formats.bondtwist2:hook0',
            'bondtwist3 = genice_bondtwist.formats.bondtwist3:hook0',
        ],
        'genice_format_hook1': [
            'bondtwist3 = genice_bondtwist.formats.bondtwist3:hook1',
        ],
        'genice_format_hook2': [
            'bondtwist  = genice_bondtwist.formats.bondtwist:hook2',
            'bondtwist2 = genice_bondtwist.formats.bondtwist2:hook2',
        ],
    },
    install_requires=['genice>=0.23', 'genice-svg'],

    license='MIT',
)
