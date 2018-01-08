# Import modules
import gzip
import csv
import daiquiri

from collections import defaultdict
from intervaltree import IntervalTree


def check_inputs(input_regions, input_mutations):
    """
    Check input format
    :param input_regions: path to input genomic regions
    :param input_mutations: path to file containing mutations
    :return:
            gz: bool, True if '.gz' found in input_file else False
    """

    # Check input file compression
    gz = True if '.gz' in input_mutations else False

    read_function = open if gz is False else gzip.open
    mode = 'r' if gz is False else 'rt'

    # magic_dict = {
    #     "\x1f\x8b\x08": "gz",
    #     "\x42\x5a\x68": "bz2",
    #     "\x50\x4b\x03\x04": "zip"
    #     }
    # max_len = max(len(x) for x in magic_dict)
    #
    # with read_function(input_mutations, mode) as csvfile:
    #     file_start = csvfile.read(max_len)
    #     for magic, filetype in magic_dict.items():
    #         if file_start.startswith(magic):
    #             print(filetype)
    #         else:
    #             print('no match')

    """
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }
    
    max_len = max(len(x) for x in magic_dict)
    
    def file_type(filename):
        with open(filename) as f:
            file_start = f.read(max_len)
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype
        return "no match"
    """

    # TODO: pre-processing
    """
    Regions non-overlapping for the same id
    Regions check double ensembl ids
    Regions check if different symbols have repeated ids
    Mutations, Tabular? Header? Extension? 
    """

    return gz


def read_regions(input_regions, elements):
    """
    Parse input genomic regions
    :param input_regions: path to input genomic regions
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :return:
        trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome.
        regions_d: dictionary containing a list of tuples with the genomic positions of each element.
        chromosomes_d: dict, keys are elements, values are chromosomes
    """

    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(list)
    chromosomes_d = defaultdict()

    with gzip.open(input_regions, 'rb') as fd:
        for line in fd:
            line = line.decode()  # binary to readable
            chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
            if elements and symbol not in elements:
                continue
            trees[chromosome][int(start): int(end) + 1] = symbol + '_' + ensid  # int, +1 end
            regions_d[symbol + '_' + ensid].append((int(start), int(end) + 1))
            chromosomes_d[symbol + '_' + ensid] = chromosome

    if not regions_d.keys():
        logger.critical('No elements found in genomic regions. Please, check input data')
        quit()

    return regions_d, chromosomes_d, trees


def read_mutations(input_mutations, trees, gz):
    """
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path to file containing mutations
    :param trees: dictionary of dictionary of intervaltrees containing intervals of genomic elements per chromosome
    :param gz: bool, True if '.gz' found in input_file else False
    :return:
        mutations_d: dictionary, key = element, value = list of mutations per element
    """
    mutations_d = defaultdict(list)

    read_function = open if gz is False else gzip.open
    mode = 'r' if gz is False else 'rt'
    if '.tab' in input_mutations:
        delimiter = '\t'
    elif '.txt' in input_mutations:
        delimiter = '\t'
    elif '.in' in input_mutations:
        delimiter = '\t'
    elif '.csv' in input_mutations:
        delimiter = ','
    else:
        logger.critical('Mutations file format not recognized. Please, provide tabular or csv formatted data')
        quit()

    with read_function(input_mutations, mode) as csvfile:
        fd = csv.DictReader(csvfile, delimiter=delimiter)
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
                            mutations_d[res.data].append(int(position))
    return mutations_d


def parse(input_regions, elements, input_mutations):
    """Parse genomic regions and dataset of cancer type mutations
    :param input_regions: path to file containing mutational data
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param input_mutations: path tab file
    :return
        regions_d: dictionary containing a list of tuples with the genomic positions of each element.
        chromosomes_d: dict, keys are elements, values are chromosomes
        mutations_d: dictionary, key = element, value = list of mutations per element
        gz: bool, True if '.gz' found in input_file else False
    """

    global logger
    logger = daiquiri.getLogger()

    gz = check_inputs(input_regions, input_mutations)
    logger.debug('Pre-processing done')
    regions_d, chromosomes_d, trees = read_regions(input_regions, elements)
    logger.debug('Regions parsed')
    mutations_d = read_mutations(input_mutations, trees, gz)
    logger.debug('Mutations parsed')

    return regions_d, chromosomes_d, mutations_d, gz
