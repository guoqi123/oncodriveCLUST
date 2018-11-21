.. _readme:

OncodriveCLUSTL
================

One of the main goals of cancer research is the identification of the genomic elements that drive tumorigenesis. OncodriveCLUSTL is a new nucleotide sequence-based clustering algorithm to detect cancer drivers in genomic regions. OncodriveCLUSTL is based on a local background model derived from the nucleotide context mutational probabilities of the chort under study. Our method is able to identify well-known cancer drivers in coding regions. It can be applied to non-coding regions and non-human data.
.. _readme license:

License
-------

OncodriveCLUSTL is available to the general public subject to certain conditions described in its `LICENSE <LICENSE>`_.

.. _readme install:

Installation
------------

OncodriveCLUSTL depends on Python 3.5 and some external libraries. First you need to download and uncompress the oncodriveclustl-0.1.0.tar.gz file.
Then you can install it with ``pip``::

        $ cd oncodriveclustl-0.1.0/
        $ pip install .

We recommend using `conda <https://www.anaconda.com/download/>`_ to install Python 3.5 and OncodriveCLUSTL.

.. note::

    The first time that you run OncodriveCLUSTL it will download the genome reference from our servers. By default the
    downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the
    system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

        $ oncodriveclustl --help

.. _readme example:

Run the example
---------------

After installing OncodriveCLUSTL, you can run an example of TCGA breast adenocarcinomas (Weinstein et al., 2013) for coding regions (Mularoni et al., 2016) using 1000 simulations.
First you need to download and uncompress the example.tar.xz file. Then you can type::

        $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example -sim region_restricted --concatenate -n 1000

The results will be saved in a folder named ``output_example``.

If you want to run a specific gene, you can type::

        $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example -sim region_restricted --concatenate -n 1000 -e PIK3CA

You can plot the observed clusters for a gene::

        $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example --concatenate -e PIK3CA --plot

