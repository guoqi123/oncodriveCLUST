# Import modules
import gzip
import csv
import daiquiri

from collections import defaultdict
from intervaltree import IntervalTree

from oncodriveclustl.utils import preprocessing as prep


def read_regions(input_regions, elements):
    """
    Parse input genomic regions
    :param input_regions: path to input genomic regions
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :return:
        trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome.
        regions_d: dictionary of IntervalTrees with genomic regions for elements
        chromosomes_d: dict, keys are elements, values are chromosomes
    """

    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(IntervalTree)
    chromosomes_d = defaultdict()
    strands_d = defaultdict()

    comp = prep.check_compression(input_regions)

    if comp == 'gz':
        with gzip.open(input_regions, 'rb') as fd:
            for line in fd:
                line = line.decode()  # binary to readable
                chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
                if elements and symbol not in elements:
                    continue
                trees[chromosome][int(start): int(end) + 1] = symbol + '_' + ensid  # int, +1 end
                regions_d[symbol + '_' + ensid].addi(int(start), (int(end) + 1))
                chromosomes_d[symbol + '_' + ensid] = chromosome
                strands_d[symbol + '_' + ensid] = strand

        if not regions_d.keys():
            logger.critical('No elements found in genomic regions. Please, check input data')
            quit()
    else:
        logger.critical('Genomic regions are not compressed, please input .gz file')
        quit()

    return regions_d, chromosomes_d, strands_d, trees


def read_mutations(input_mutations, trees):
    """
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path to file containing mutations
    :param trees: dictionary of dictionary of intervaltrees containing intervals of genomic elements per chromosome
    :return:
        mutations_d: dictionary, key = element, value = list of mutations per element
    """
    mutations_d = defaultdict(list)
    read_function, mode, delimiter = prep.check_tabular_csv(input_mutations)

    with read_function(input_mutations, mode) as read_file:
        fd = csv.DictReader(read_file, delimiter='\t')
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
    """

    global logger
    logger = daiquiri.getLogger()
    regions_d, chromosomes_d, strands_d, trees = read_regions(input_regions, elements)
    logger.debug('Regions parsed')
    mutations_d = read_mutations(input_mutations, trees)
    logger.debug('Mutations parsed')

    return regions_d, chromosomes_d, strands_d, mutations_d
