from setuptools import setup, find_packages
from oncodriveclustl import __version__, __author__, __author_email__, __description__, __long_description__


install_requires = [
    'intervaltree',
    'numpy',
    'scipy',
    'pandas',
    'statsmodels',
    'bgreference',
    'click',
    'daiquiri',
    'tqdm',
    'matplotlib'
]


setup(
    name="oncodriveclustl",
    version=__version__,
    packages=find_packages(),
    package_data={'oncodriveclustl': ['data/*.tsv']},
    author=__author__,
    author_email=__author_email__,
    description=__description__,
    license="Apache License 2",
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
