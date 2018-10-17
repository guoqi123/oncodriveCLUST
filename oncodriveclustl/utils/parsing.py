"""
Contains functions to parse input mutations and genomic regions files
"""

import gzip
import csv
from collections import defaultdict
from collections import namedtuple

import daiquiri
from intervaltree import IntervalTree

from oncodriveclustl.utils import exceptions as excep
from oncodriveclustl.utils import preprocessing as prep

Mutation = namedtuple('Mutation', 'position, region, alt, sample, cancertype')
Cds = namedtuple('Cds', 'start, end')


def read_regions(input_regions, elements):
    """
    Parse input genomic regions

    Args:
        input_regions (str): path to input genomic regions
        elements (set): elements to analyze. If the set is empty all the elements in genomic regions will be analyzed

    Returns:
        trees (dict): dictionary of dictionary of intervaltrees containing intervals of genomic elements by chromosome.
        regions_d (dict): dictionary of IntervalTrees with genomic regions for elements
        chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
        strands_d (dict): dictionary of elements (keys) and strands (values)

    """
    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(IntervalTree)
    chromosomes_d = defaultdict()
    strands_d = defaultdict()
    comp = prep.check_compression(input_regions)

    if comp == 'gz':
        with gzip.open(input_regions, 'rt') as fd:
            for line in fd:
                chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
                if elements and symbol not in elements:
                    continue
                if int(start) != int(end):
                    trees[chromosome][int(start): int(end) + 1] = symbol + '//' + ensid
                    regions_d[symbol + '//' + ensid].addi(int(start), (int(end) + 1))
                    chromosomes_d[symbol + '//' + ensid] = chromosome
                    strands_d[symbol + '//' + ensid] = strand
        if not regions_d.keys():
            raise excep.UserInputError('No elements found in genomic regions. Please, check input data')
    else:
        raise excep.UserInputError('Genomic regions are not compressed, please input GZIP compressed file')

    return regions_d, chromosomes_d, strands_d, trees


def map_regions_concatseq(regions_d):
    """
    Calculate position (index) of every region relative to genomic element start

    Args:
        regions_d (dict): dictionary of IntervalTrees with genomic regions for elements

    Returns:
        concat_regions_d (dict): dictionary of dictionaries with relative index of genomic regions

    """
    global Cds
    concat_regions_d = defaultdict(dict)

    for element, regions in regions_d.items():
        start = 0
        for region in sorted(regions):
            length = region.end - region.begin
            end = start + length - 1
            concat_regions_d[element][region.begin] = Cds(start, end)
            start = end + 1

    return concat_regions_d


def read_mutations(input_mutations, trees, is_pancancer):
    """
    Read mutations file (only substitutions) and map to elements' genomic regions

    Args:
        input_mutations (str): path to input file containing mutations
        trees (dict): dictionary of dictionary of IntervalTrees containing intervals of genomic elements per chromosome
        is_pancancer (bool): True, analysis carried out using cohorts available in the input mutations file

    Returns:
        mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
        samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
        cohorts_d (dict): dictionary of elements (keys) and set of cohorts containing element mutations (values)

    """
    global Mutation
    mutations_d = defaultdict(list)
    samples_d = defaultdict(int)
    cohorts_d = defaultdict(set)
    read_function, mode, delimiter, cancer_type_header = prep.check_tabular_csv(input_mutations)
    file_prefix = input_mutations.split('/')[-1].split('.')[0]

    with read_function(input_mutations, mode) as read_file:
        fd = csv.DictReader(read_file, delimiter=delimiter)
        for line in fd:
            chromosome = line['CHROMOSOME']
            position = int(line['POSITION'])
            ref = line['REF']
            alt = line['ALT']
            sample = line['SAMPLE']
            if cancer_type_header and is_pancancer:
                cancer_type = line['CANCER_TYPE']
            else:
                cancer_type = file_prefix
            samples_d[sample] += 1
            # Read substitutions only
            if len(ref) == 1 and len(alt) == 1:
                if ref != alt:
                    if ref != '-' and alt != '-':
                        if trees[chromosome][int(position)] != set():
                            results = trees[chromosome][int(position)]
                            for res in results:
                                m = Mutation(position, (res.begin, res.end), alt, sample, cancer_type)
                                mutations_d[res.data].append(m)
                                cohorts_d[res.data].add(cancer_type)

    return mutations_d, samples_d, cohorts_d


def parse(input_regions, elements, input_mutations, concatenate, is_pancancer):
    """Parse genomic regions and dataset of cancer type mutations

    Args:
        input_regions (str): path to input genomic regions
        elements (set): elements to analyze. If the set is empty all the elements in genomic regions will be analyzed
        input_mutations (str): path to file containing mutations
        concatenate (bool): True calculates clustering on collapsed genomic regions (e.g., coding regions in a gene)
        is_pancancer (bool): True, analysis carried out using cohorts available in the input mutations file

    Returns:
        regions_d (dict): dictionary of IntervalTrees containing genomic regions from all analyzed elements
        concat_regions_d (dict): dictionary of dictionaries with relative to start position (index) of genomic regions
        chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
        strands_d (dict): dictionary of elements (keys) and strands (values)
        mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
        samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
        cohorts_d (dict): dictionary of elements (keys) and set of cohorts containing element mutations (values)

    """
    global logger
    logger = daiquiri.getLogger()

    regions_d, chromosomes_d, strands_d, trees = read_regions(input_regions, elements)
    if concatenate:
        concat_regions_d = map_regions_concatseq(regions_d)
    else:
        concat_regions_d = {}
    logger.info('Regions parsed')
    mutations_d, samples_d, cohorts_d = read_mutations(input_mutations, trees, is_pancancer)
    logger.info('Mutations parsed')

    return regions_d, concat_regions_d, chromosomes_d, strands_d, mutations_d, samples_d, cohorts_d
