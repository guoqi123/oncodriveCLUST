# Import modules
import gzip
import csv
import daiquiri
import pickle

from collections import defaultdict
from intervaltree import IntervalTree
from os import path


def read_regions(input_regions, elements):
    """
    Parse input regions
    :param input_regions: path tab file chr \t start \t end \t strand  \t ensembl id \t ensembl id \t symbol
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :return:
        trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome.
        regions_d: dictionary containing a list of tuples with the genomic positions of each element.
        chromosomes_d: dict, keys are elements, values are chromosomes
    """

    BAD_ENSEMBL_IDS = ('ENSGR',)  # pseudoautosomal regions

    global logger
    logger = daiquiri.getLogger()

    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(list)
    chromosomes_d = defaultdict()

    # TODO remove this hardcoded file?
    with open(path.join(path.dirname(__file__), '../data/02_cds.regions.mapping.pickle'), 'rb') as fd:
        dict_ids = pickle.load(fd)

    with gzip.open(input_regions, 'rb') as fd:  # rb for binary <-- gz
        for line in fd:
            line = line.decode()  # binary to readable
            chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
            if ensid == dict_ids[symbol]:
                if elements and symbol not in elements:
                    continue
                if all([i not in ensid for i in BAD_ENSEMBL_IDS]):
                    trees[chromosome][int(start): int(end) + 1] = symbol  # int, +1 end
                    regions_d[symbol].append((int(start), int(end) + 1))
                    chromosomes_d[symbol] = chromosome
            else:
                pass

    return regions_d, chromosomes_d, trees


def read_mutations(input_mutations, trees, gz):
    """
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path tab file
    :param trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome
    :param gz: bool, True if '.gz' found in input_file else False
    :return:
        mutations_d: dictionary, key = element, value = list of mutations per element
    """
    mutations_d = defaultdict(list)

    read_function = open if gz is False else gzip.open
    mode = 'r' if gz is False else 'rt'

    with read_function(input_mutations, mode) as csvfile:
        fd = csv.DictReader(csvfile, delimiter='\t')
        for line in fd:
            chromosome = line['CHROMOSOME']
            position = int(line['POSITION'])
            ref = line['REF']
            alt = line['ALT']
            # Read substitutions only
            if len(ref) == 1 and len(alt) == 1:
                if ref != '-' and alt != '-':
                    if trees[chromosome][int(position)] != set():
                        results = trees[chromosome][int(position)]
                        for res in results:
                            mutations_d[res.data].append(int(position))  # dict: key = element, value = list of mutations
    return mutations_d


def parse(input_regions, elements, input_mutations, gz):
    """Parse genomic regions and dataset of cancer type mutations
    :param input_regions: path tab file chr \t start \t end \t strand  \t ensembl id \t ensembl id \t symbol
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param input_mutations: path tab file
    :param gz: bool, True if '.gz' found in input_file else False
    :return
        regions_d: dictionary containing a list of tuples with the genomic positions of each element.
        chromosomes_d: dict, keys are elements, values are chromosomes
        mutations_d: dictionary, key = element, value = list of mutations per element
    """

    regions_d, chromosomes_d, trees = read_regions(input_regions, elements)
    mutations_d = read_mutations(input_mutations, trees, gz)
    return regions_d, chromosomes_d, mutations_d
