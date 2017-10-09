# Import modules
import gzip
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

def read_mutations(input_mutations, trees):
    """
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path tab file
    :param trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome.    :return:
        mutations_d: dictionary, key = element, value = list of mutations per element
    """
    mutations_d = defaultdict(list)

    with open(input_mutations, 'r') as fd:
        next(fd)  # don't read header
        for line in fd:
            chromosome, position, ref, alt, sample, mut_type, tumor = line.strip().split('\t')[:7]
            if trees[chromosome][int(position)] != set() and mut_type == 'subs':  # intersect only subs
                results = trees[chromosome][int(position)]
                for res in results:
                    mutations_d[res.data].append(int(position))  # dict: key = element, value = list of mutations
    return mutations_d


def parse(input_regions, elements, input_mutations):
    """Parse genomic regions and dataset of cancer type mutations
    :param input_regions: path tab file chr \t start \t end \t strand  \t ensembl id \t ensembl id \t symbol
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param input_mutations: path tab file
    :return
        regions_d: dictionary containing a list of tuples with the genomic positions of each element.
        chromosomes_d: dict, keys are elements, values are chromosomes
        mutations_d: dictionary, key = element, value = list of mutations per element
    """

    regions_d, chromosomes_d, trees = read_regions(input_regions, elements)
    mutations_d = read_mutations(input_mutations, trees)

    return regions_d, chromosomes_d, mutations_d