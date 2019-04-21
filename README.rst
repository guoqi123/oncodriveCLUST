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

OncodriveCLUSTL depends on Python 3.5 and some external libraries. We recommend to install it using `conda <https://www.anaconda.com/download/>`_::

        $ conda install -c bbglab oncodriveclustl


.. note::

    The first time that you run OncodriveCLUSTL it will download the genome reference from our servers. By default the
    downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the
    system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

        $ oncodriveclustl --help


.. _readme inputdata:

Input data
---------------
OncodriveCLUSTL only requires two main inputs:

- Mutations file. TSV file containing SNVs from whole exome (WXS) or whole genome (WGS) data mapped to a reference genome (e.g., human hg19 or mouse c3h). This file must contain, at least, the following columns with header in the following order:

    1. CHROMOSOME: 1, 2,..., X, Y
    2. POSITION: Integer indicating the position of the mutation
    3. REF: Reference nucleotide
    4. ALT: Alternate nucleotide
    5. SAMPLE: Identifier of the sample

- Annotations file. TSV GZIP compressed file without headers containing the genomic coordinates of genomic elements (GEs):

    1. Chromosome: 1, 2,..., X, Y
    2. Start: Starting position of the genomic region
    3. End: Final position of the genomic region
    4. Strand: Strand of the genomic region ("+" or "-")
    5. Element ID: Identifier of the genomic element.
    6. Transcript ID: Identifier of the transcript
    7. Symbol: Symbol of the genomic element

OncodriveCLUSTL will analyze genomic elements as Symbol + Element ID.

.. note::
    Coordinates of a given GE cannot overlap.

You can check the input formats in the files provided in the example.

If you have a directory containing mutation files in VCF format, you can run our VCF parser to obtain a tabular file compatible with OncodriveCLUSTL input format::

       $ parse_vcf -i [INPUT_DIRECTORY] -o [OUTPUT_FILE]

Please, check 'parsers/vcf.py' module for more details.

.. _readme outputdata:

Output data
---------------
OncodriveCLUSTL generates three output files:

- Elements results file ('elements_results.txt'). TSV file containing results of the analyzed elements:

    1. SYMBOL: GE symbol
    #. ENSID: GE ID
    #. CGC: True if GE in the CGC list
    #. CHROMOSOME: 1, 2,..., X, Y
    #. STRAND: Strand of the GE ("+" or "-")
    #. LENGTH: length (bp) of the GE
    #. TOTAL_MUT: total mutations observed in the GE
    #. CLUSTERED_MUT: number of mutations in a cluster
    #. CLUSTERS: number of clusters
    #. SIM_CLUSTERS: number of simulated clusters
    #. SCORE: GE score
    #. P_EMPIRICAL: empirical p-value of the GE
    #. Q_EMPIRICAL: empirical q-value of the GE
    #. P_ANALYTICAL: analytical p-value of the GE
    #. Q_ANALYTICAL: analytical q-value of the GE
    #. P_TOPCLUSTER: analytical p-value of the cluster with highest cluster score
    #. Q_TOPCLUSTER: analytical q-value of the cluster with highest cluster score


- Clusters results file ('clusters_results.tsv'). TSV file containing results of the clusters observed in the analyzed elements:

    1. RANK: Position of the GE in the list of
    #. SYMBOL: GE symbol
    #. ENSID: GE ID
    #. CGC: True if GE in the CGC list
    #. CHROMOSOME: 1, 2,..., X, Y
    #. STRAND: Strand of the GE ("+" or "-")
    #. COORDINATES: genomic coordinates of the cluster. It can be 'coord1,coord2' for clusters inside a single region or 'coord1,coord2;coord3,coord4' for those spanning regions (--concatenate flag)
    #. MAX_COORD: genomic position with the highest smoothing score inside the cluster
    #. WIDTH: cluster's width (pb)
    #. N_MUT: number of mutations in the cluster
    #. N_SAMPLES: number of samples with a mutation in the cluster
    #. FRA_UNIQ_SAMPLES: proportion of unique samples mutated in the cluster out of the total of mutations in the cluster
    #. SCORE: cluster score
    #. P: analytical p-value of the cluster

- Log file ('results.log'). TXT file containing OncodriveCLUSTL's run information


.. _readme commandline:

Command line
---------------
- '-i', '--input-file': File containing mutations (required)
- '-r', '--regions-file': GZIP compressed file with the genomic regions to analyze (required)
- '-o', '--output-directory': Output directory to be created (required)
- '-sign', '--input-signature': File containing input context based mutational probabilities
- '-ef', '--elements-file': File with the symbol of a set elements to analyze, one per row
- '-e', '--elements': Symbol of the element to analyze
- '-g', '--genome': Genome to use. Default is hg19.
- '-emut', '--element-mutations': Cutoff of element mutations. Default is 2
- '-cmut', '--cluster-mutations': Cutoff of cluster mutations. Default is 2
- '-sw', '--smooth-window': Smoothing window. Default is 11
- '-cw', '--cluster-window': Cluster window. Default is 11
- '-kmer', '--kmer': Kmer-nucleotide context (3 or 5)
- '-n', '--n-simulations': Number of simulations. Default is 1000
- '-sim', '--simulation-mode': Simulation mode. Default is 'mutation_centered'
- '-simw', '--simulation-window': Simulation window. Default is 31
- '-c', '--cores': Number of cores to use in the computation. By default it uses all the available cores
- '--log-level': Verbosity of the logger. Default is 'info'
- '--concatenate': Calculate clustering on concatenated genomic regions (e.g., exons in coding sequences)
- '--groupby': Analysis carried out by groups (e.g., PanCancer cohort analysis).
- '--clustplot': Generate a needle plot with clusters for an element
- '--qqplot': Generate a quantile-quantile (QQ) plot for a dataset
- '--gzip': Gzip compress files

.. note::
    When using simulation mode 'mutation_centered', simulation windows can be simulated outside the genomic element.

.. note::
    When using '--groupby' flag, input mutation file requires an extra column "GROUP_BY" (group identifier). According to it, one mutational signature will be computed for each group. Please check that the number of mutations per group is sufficient for an accurate signatures calculation.

.. _readme example:

Run the example
---------------

After installing OncodriveCLUSTL, you can run an example of TCGA pancreatic adenocarcinomas (Ellrott et al. 2018) for coding regions (Mularoni et al., 2016) using 1000 simulations.
First you need to download the example folder. Then you run OncodriveCLUSTL with default mode and parameters as::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example

The results will be saved in a folder named ``output_example``.

You can compute a more sophisticated analysis using non-default parameters and generate a quantile-quantile plot by typing::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example -sw 15 -cw 15 -simw 35 -sim region_restricted --concatenate --qqplot

If you want to run a specific GE and generate a plot its observed clusters, you can type::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example -sw 15 -cw 15 -simw 35 -sim region_restricted --concatenate --clustplot -e KRAS


