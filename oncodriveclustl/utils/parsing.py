"""
Contains functions to parse input mutations and genomic regions files
"""

from collections import defaultdict
from collections import namedtuple

from bgparsers import readers
import daiquiri
from intervaltree import IntervalTree

from oncodriveclustl.utils import exceptions as excep

Mutation = namedtuple('Mutation', 'position, region, alt, sample, group')
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

    # Read regions
    for row in readers.elements(
            input_regions,
            required=['CHROMOSOME', 'START', 'END', 'ELEMENT', 'SYMBOL'],
            extra=['STRAND']
    ):
        chromosome = row['CHROMOSOME']
        start = row['START']
        end = row['END']
        ensid = row['ELEMENT']
        symbol = row['SYMBOL']
        strand = row.get('STRAND', '+')  # If None, '+' by default

        if elements and symbol not in elements:
            continue
        if int(start) != int(end):
            trees[chromosome][int(start): int(end) + 1] = symbol + '//' + ensid
            regions_d[symbol + '//' + ensid].addi(int(start), (int(end) + 1))
            chromosomes_d[symbol + '//' + ensid] = chromosome
            strands_d[symbol + '//' + ensid] = strand
    # Check
    if not regions_d.keys():
        raise excep.UserInputError('No elements found in genomic regions. '
                                   'Please, check input data ({})'.format(input_regions))

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


def read_mutations(input_mutations, trees, signature_group):
    """
    Read mutations file (only substitutions) and map to elements' genomic regions

    Args:
        input_mutations (str): path to input file containing mutations
        trees (dict): dictionary of dictionary of IntervalTrees containing intervals of genomic elements per chromosome
        signature_group (str): header of the column to group signatures calculation

    Returns:
        mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
        samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
        groups_d (dict): dictionary of elements (keys) and groups (values)

    """
    global Mutation
    mutations_d = defaultdict(list)
    samples_d = defaultdict(int)
    groups_d = defaultdict(set)
    file_prefix = input_mutations.split('/')[-1].split('.')[0]

    # Read mutations
    for row in readers.variants(
            input_mutations,
            required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
            extra=['SIGNATURE', 'CANCER_TYPE']
    ):
        chromosome = row['CHROMOSOME']
        position = int(row['POSITION'])
        ref = row['REF']
        alt = row['ALT']
        sample = row['SAMPLE']
        samples_d[sample] += 1
        if signature_group:
            try:
                group = row[signature_group]
            except KeyError:
                raise excep.UserInputError('Mutation file does not contain column "{}" '.format(signature_group))
        else:
            group = file_prefix
        # Read substitutions only
        if len(ref) == 1 and len(alt) == 1:
            if ref != alt:
                if ref != '-' and alt != '-':
                    if trees[chromosome][int(position)] != set():
                        results = trees[chromosome][int(position)]
                        for res in results:
                            m = Mutation(position, (res.begin, res.end), alt, sample, group)
                            mutations_d[res.data].append(m)
                            groups_d[res.data].add(group)

    return mutations_d, samples_d, groups_d


def parse(input_regions, elements, input_mutations, concatenate, signature_group):
    """Parse genomic regions and dataset of cancer type mutations

    Args:
        input_regions (str): path to input genomic regions
        elements (set): elements to analyze. If the set is empty all the elements in genomic regions will be analyzed
        input_mutations (str): path to file containing mutations
        concatenate (bool): True calculates clustering on collapsed genomic regions (e.g., coding regions in a gene)
        signature_group (str): header of the column to group signatures calculation

    Returns:
        regions_d (dict): dictionary of IntervalTrees containing genomic regions from all analyzed elements
        concat_regions_d (dict): dictionary of dictionaries with relative to start position (index) of genomic regions
        chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
        strands_d (dict): dictionary of elements (keys) and strands (values)
        mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
        samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
        groups_d (dict): dictionary of elements (keys) and groups (values)

    """
    global logger
    logger = daiquiri.getLogger()

    regions_d, chromosomes_d, strands_d, trees = read_regions(input_regions, elements)
    if concatenate:
        concat_regions_d = map_regions_concatseq(regions_d)
    else:
        concat_regions_d = {}
    logger.info('Regions parsed')
    mutations_d, samples_d, groups_d = read_mutations(input_mutations, trees, signature_group)
    logger.info('Mutations parsed')

    return regions_d, concat_regions_d, chromosomes_d, strands_d, mutations_d, samples_d, groups_d
