import sys
from setuptools import setup, find_packages

from oncodriveclustl import __version__, __author__, __author_email__, __description__, __long_description__


# Check the python compatibility
if sys.hexversion < 0x03040000:
    raise RuntimeError('This package requires Python 3.4 or later.')

install_requires = [
    'intervaltree>=2.1.0',
    'numpy>=1.13.3',
    'scipy>=1.0.0',
    'pandas>=0.22.0',
    'statsmodels>=0.8.0',
    'bgreference>=0.5',
    'click>=6.7',
    'daiquiri>=1.3.0',
    'tqdm>=4.19.4',
    'matplotlib>=2.0.2'
]


setup(
    name="oncodriveclustl",
    python_requires='>3.4',
    version=__version__,
    packages=find_packages(),
    package_data={'oncodriveclustl': ['data/*.tsv']},
    author=__author__,
    author_email=__author_email__,
    description=__description__,
    license="GPLv3",
    keywords="",
    url="https://bitbucket.org/bbglab/oncodriveclustl",
    download_url="https://bitbucket.org/bbglab/oncodriveclustl/get/" + __version__ + ".tar.gz",
    long_description=__long_description__,
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'oncodriveclustl = oncodriveclustl.main:main',
        ]
    }
)
