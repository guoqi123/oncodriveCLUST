"""
Contains functions to parse input mutations and genomic regions files
"""

import gzip
import csv
from collections import defaultdict
from collections import namedtuple

import daiquiri
from intervaltree import IntervalTree
import bgreference as bg
import bgdata as bgd
import pickle

from oncodriveclustl.utils import exceptions as excep
from oncodriveclustl.utils import preprocessing as prep

Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample, cancertype')
Cds = namedtuple('Cds', 'start, end')


def read_regions(input_regions, elements, protein):
    """
    Parse input genomic regions

    Args:
        input_regions (str): path to input genomic regions
        elements (set): elements to analyze. If the set is empty all the elements in genomic regions will be analyzed
        protein (bool): True reads transcripts id instead of elements id

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
                if not protein:
                    chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
                else:
                    chromosome, start, end, strand, _, ensid, symbol = line.strip().split('\t')
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


def map_regions_cds(regions_d):
    """
    Calculate cds position of every region relative to genomic element start

    Args:
        regions_d (dict): dictionary of IntervalTrees with genomic regions for elements

    Returns:
        cds_d (dict): dictionary of dictionaries with relative cds index of genomic regions

    """
    global Cds
    cds_d = defaultdict(dict)

    for element, regions in regions_d.items():
        start = 0
        for region in sorted(regions):
            length = region.end - region.begin
            end = start + length - 1
            cds_d[element][region.begin] = Cds(start, end)
            start = end + 1

    return cds_d


def read_mutations(input_mutations, trees, conseq):
    """
    Read mutations file (only substitutions) and map to elements' genomic regions

    Args:
        input_mutations (str): path to file containing mutations
        trees (dict): dictionary of dictionary of IntervalTrees containing intervals of genomic elements per chromosome
        conseq (bool): True, use AA consequence type

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

    if conseq:
        conseq_path = bgd.get_path('oncodriveclustl', 'vep88', 'hg19_canonical_conseq')
        with read_function(input_mutations, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            for line in fd:
                chromosome = line['CHROMOSOME']
                position = int(line['POSITION'])
                ref = line['REF']
                alt = line['ALT']
                sample = line['SAMPLE']
                if cancer_type_header:
                    cancer_type = line['CANCER_TYPE']
                else:
                    cancer_type = file_prefix
                samples_d[sample] += 1
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if trees[chromosome][int(position)] != set():
                            results = trees[chromosome][int(position)]
                            for res in results:
                                ensid = res.data.split('//')[1]
                                path_to_vep_pickle = conseq_path + '/{}.pickle'.format(ensid)
                                try:
                                    with open(path_to_vep_pickle, 'rb') as fd:
                                        conseq_d = pickle.load(fd)
                                        muttype = 0 if position in conseq_d.get(alt, []) else 1
                                except FileNotFoundError as e:
                                    logger.error(
                                        '{}\nVep file for element {} could not be read. Analysis will be done without '
                                        'considering mutations consequence type\n'.format(e, res.data))
                                    muttype = 1
                                m = Mutation(position, (res.begin, res.end), alt, muttype, sample, cancer_type)
                                mutations_d[res.data].append(m)
                                cohorts_d[res.data].add(cancer_type)
    else:
        with read_function(input_mutations, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            for line in fd:
                chromosome = line['CHROMOSOME']
                position = int(line['POSITION'])
                ref = line['REF']
                alt = line['ALT']
                sample = line['SAMPLE']
                if cancer_type_header:
                    cancer_type = line['CANCER_TYPE']
                else:
                    cancer_type = file_prefix
                samples_d[sample] += 1
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if trees[chromosome][int(position)] != set():
                            results = trees[chromosome][int(position)]
                            for res in results:
                                muttype = 1
                                m = Mutation(position, (res.begin, res.end), alt, muttype, sample, cancer_type)
                                mutations_d[res.data].append(m)
                                cohorts_d[res.data].add(cancer_type)

    return mutations_d, samples_d, cohorts_d


def map_transcripts_protein(regions_d, chromosomes_d, strands_d, genome):
    """
    Map transcript to reference protein sequence. Remove transcripts that do not pass control check.

    Args:
    regions_d (dict): dictionary of IntervalTrees with genomic regions for elements
    chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
    strands_d (dict): dictionary of elements (keys) and strands (values)
    genome (str): genome to use

    Returns:
    regions_d (dict): dictionary of IntervalTrees with genomic regions for elements updated

    """
    reverse_d = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    start_codons = ['ATG']  # TODO consider 'TTG', 'CTG'
    elements_to_skip = set()

    for element, regions in regions_d.items():
        cds = []
        # Get nucleotide sequence
        for interval in sorted(regions, reverse=False):
            start = interval[0]
            size = interval[1] - interval[0]  # no need + 1 because interval end already + 1
            sequence = bg.refseq(genome, chromosomes_d[element], start, size)
            if sequence.count('N') == 0:
                cds.extend(sequence)
            else:
                logger.warning('Found N nucleotide in {}. Element is discarded from analysis\n'.format(element))
                elements_to_skip.add(element)
                break
        # Reverse complementary if needed
        if strands_d[element] == '-':
            cds.reverse()
            cds = ''.join([reverse_d.get(i, i) for i in cds])
        else:
            cds = ''.join(cds)
        # Check
        if len(cds) % 3 == 0 and cds[:3] in start_codons:
            pass
        else:
            logger.warning('{} cannot be translated to protein. Element is discarded from analysis'.format(element))
            elements_to_skip.add(element)

    return elements_to_skip


def parse(input_regions, elements, input_mutations, cds, conseq, protein, genome):
    """Parse genomic regions and dataset of cancer type mutations

    Args:
        input_regions (str): path to input genomic regions
        elements (set): elements to analyze. If the set is empty all the elements in genomic regions will be analyzed
        input_mutations (str): path to file containing mutations
        cds (bool): True calculates clustering on collapsed genomic regions (e.g., coding regions in a gene)
        conseq (bool): True, use AA consequence type
        protein (bool): True analyzes clustering in translated protein sequences
        genome (str): genome to use

    Returns:
        regions_d (dict): dictionary of IntervalTrees containing genomic regions from all analyzed elements
        cds_d (dict): dictionary of dictionaries with relative cds index of genomic regions
        chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
        strands_d (dict): dictionary of elements (keys) and strands (values)
        mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
        samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
        cohorts_d (dict): dictionary of elements (keys) and set of cohorts containing element mutations (values)

    """
    global logger
    logger = daiquiri.getLogger()

    regions_d, chromosomes_d, strands_d, trees = read_regions(input_regions, elements, protein)
    if cds:
        cds_d = map_regions_cds(regions_d)
    else:
        cds_d = {}
    logger.info('Regions parsed')
    mutations_d, samples_d, cohorts_d = read_mutations(input_mutations, trees, conseq)
    logger.info('Mutations parsed')

    if protein:
        elements_to_skip = map_transcripts_protein(regions_d, chromosomes_d, strands_d, genome)
        logger.info('Protein sequences calculated')
        # Remove invalid transcripts
        for element in elements_to_skip:
            del regions_d[element]
            del cds_d[element]
            mutations_d.pop(element, None)

    return regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, cohorts_d
