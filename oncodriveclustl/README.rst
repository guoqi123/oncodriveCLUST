.. _readme:

OncodriveCLUSTL
============

.. _readme license:

License
-------

.. _readme install:

Installation
------------

OncodriveCLUSTL depends on Python 3.5 and some external libraries. You can get the latest code from the repository and install with ``pip``::

        $ git clone git@bitbucket.org:bbglab/oncodriveclustl.git
        $ cd oncodriveclustl
        $ pip install .

.. note::

   OncodriveCLUSTL has a set up dependency with `Cython <http://cython.org/>`_,
   which is required to compile the ``*.pyx`` files.

The first time that you run OncodriveCLUSTL it will download the genome reference from our servers. For the analysis of coding regions, if you specify it, it will also download our latest run from Variant Effect Prediction Tool (VEP)
version 88 (McLaren et al., 2016) from our servers to build the clusters using non-synonymous mutations. VEP tool is freely available from `<https://www.ensembl.org/info/docs/tools/vep/>`_
By default the downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

	$ oncodriveclustl --help

.. _readme example:

Run the example
---------------

Download and extract the example files (if you cloned the repository skip this step)::

   $ wget https://bitbucket.org/bbglab/oncodrivefml/downloads/oncodrivefml-examples_v2.0.tar.gz
   $ tar xvzf oncodrivefml-examples_v2.0.tar.gz

If you want to speed up the download of the genome reference that is also needed,
run this command::

   $ bg-data datasets genomereference hg19

To run the example, we have included a bash script (``run.sh``)
than will execute OncodriveCLUSTL. The script should be executed in
the folder where the files have been extracted::

   $ ./run.sh

The results will be saved in a folder named ``cds``.


.. _readme docs:

Documentation
-------------

Find OncodriveCLUSTL documentation in `ReadTheDocs <http://oncodrivefml.readthedocs.io/en/latest/>`_.

You can also compile the documentation yourself using `Sphinx <http://www.sphinx-doc.org/en/stable/>`_,
if you have cloned the repository.
To do so, install the optional packages in ``optional-requirements.txt`` and build the
documentation in the docs folder::

    $ cd docs
    $ make html
