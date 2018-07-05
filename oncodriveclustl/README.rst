.. _readme:

OncodriveCLUSTL (MSc version)
============

One of the main goals of cancer research is the identification of the genomic elements that drive tumorigenesis. This is a big challenge specially for mutations in the non-coding genome, for which new computational tools need to be addressed. We present OncodriveCLUSTL, a sequence-based clustering method for the detection of cancer drivers in any genomic region of interest  using a local simulation of the mutational process affecting them.

.. _readme license:

License
-------

OncodriveCLUSTL will be available to the general public subject to certain conditions described in its `LICENSE <LICENSE>`_.

.. _readme install:

Installation
------------

OncodriveCLUSTL depends on Python 3.6 and some external libraries. First you need to downloadn and uncompress the oncodriveclustl-0.1.0.tar.gz file.
Then you can install it with ``pip``::

        $ cd oncodriveclustl-0.1.0/
        $ pip install .

We recommend using `conda <https://www.anaconda.com/download/>`_ to install Python 3.6 and OncodriveCLUSTL.

.. note::

The first time that you run OncodriveCLUSTL it will download the genome reference from our servers. For the analysis of coding regions, if you specify it, it will also download our latest run from Variant Effect Prediction Tool (VEP)
version 88 (McLaren et al., 2016) from our servers to build the clusters using non-synonymous mutations. VEP tool is freely available from `<https://www.ensembl.org/info/docs/tools/vep/>`_
By default the downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

	$ oncodriveclustl --help

.. _readme example:

Run the example
---------------

After installing OncodriveCLUSTL, you can run an example of TCGA breast adenocarcinomas (Weinstein et al., 2013) for coding regions (Mularoni et al., 2016) using 1000 simulations.
First you need to download and uncompress the example.tar.xz file. Then you can type::

   $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example -sim region_restricted --cds -n 1000

The results will be saved in a folder named ``output_example``.

If you want to run a specific gene, you can type::
   $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example -sim region_restricted --cds -n 1000 -e PIK3CA

You can plot the observed clusters for a gene::
   $ oncodriveclustl -i example/BRCA.txt -r example/cds_regions.gz -o example/output_example --cds -e PIK3CA --plot
