import sys
from os import path
from setuptools import setup, find_packages

from oncodriveclustl import __version__

DESCRIPTION = "OncodriveCLUSTL is a clustering method to identify cancer drivers"

# Check the python compatibility
if sys.hexversion < 0x03050000:
    raise RuntimeError('This package requires Python 3.5 or later.')


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    install_requires = f.read().splitlines()


# Get the long description from the README file
with open(path.join(directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name="oncodriveclustl",
    python_requires='>=3.5',
    version=__version__,
    packages=find_packages(),
    package_data={'oncodriveclustl': ['data/*.tsv']},
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description=DESCRIPTION,
    license="AGPLv3",
    keywords="",
    url="https://bitbucket.org/bbglab/oncodriveclustl",
    download_url="https://bitbucket.org/bbglab/oncodriveclustl/get/" + __version__ + ".tar.gz",
    long_description=long_description,
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'oncodriveclustl = oncodriveclustl.main:main',
            'parse_vcf = oncodriveclustl.parsers.vcf:vcf_to_tsv',
        ]
    }
)
