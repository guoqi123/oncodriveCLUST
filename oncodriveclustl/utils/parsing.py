# Import modules
import os
import gzip
import csv
from collections import defaultdict
from collections import namedtuple
import json

import daiquiri
from intervaltree import IntervalTree
import bgreference as bg
import bgdata as bgd
import pickle

from oncodriveclustl.utils import preprocessing as prep

Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample, cancertype')
Cds = namedtuple('Cds', 'start, end')


def read_regions(input_regions, elements, protein):
    """
    Parse input genomic regions
    :param input_regions: path to input genomic regions
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param protein: bool, True reads transcripts ID instead of elements id
    :return:
        trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic elements by chromosome.
        regions_d: dictionary of IntervalTrees with genomic regions for elements
        chromosomes_d: dict, keys are elements, values are chromosomes
        strands_d: dict, keys are elements, values are strands
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
            logger.critical('No elements found in genomic regions. Please, check input data')
            quit(-1)
    else:
        logger.critical('Genomic regions are not compressed, please input GZIP compressed file')
        quit(-1)

    return regions_d, chromosomes_d, strands_d, trees


def map_regions_cds(regions_d):
    """
    Calculate cds position of every region relative to genomic element start
    :param regions_d: dictionary of IntervalTrees with genomic regions for elements
    :return:
            cds_d: dictionary of dictionaries with relative cds index of genomic regions
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
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path to file containing mutations
    :param trees: dictionary of dictionary of IntervalTrees containing intervals of genomic elements per chromosome
    :param conseq: True, use aa consequence type
    :return:
        mutations_d: dictionary, key = element, value = list of mutations formatted as namedtuple
        samples_d: dictionary, key = sample, value = number of mutations
        cohorts_d: dictionary, key = element, value = set of cohorts with element mutations
    """
    global Mutation
    mutations_d = defaultdict(list)
    samples_d = defaultdict(int)
    cohorts_d = defaultdict(set)
    read_function, mode, delimiter, cancer_type_header = prep.check_tabular_csv(input_mutations)
    file_prefix = input_mutations.split('/')[-1].split('.')[0]
    conseq_path = bgd.get_path('oncodriveclustl', 'vep88', 'hg19_canonical_conseq')

    if conseq:
        # Read mutations
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


    cds = []

    # Iterate through genomic regions to get their sequences
    for interval in sorted(regions_d[element], reverse=False):
        start = interval[0]
        size = interval[1] - interval[0]  # no need +1 because interval end already +1
        sequence = bgreference.refseq('hg19', chromosomes_d[element], start, size)
        cds.extend(sequence)

    if strands_d[element] == '-':
        # Reverse
        cds.reverse()
        # Complementary
        return ''.join([reverse_d.get(i, i) for i in cds])
    #         return Seq(''.join(cds)).reverse_complement()
    else:
        return ''.join(cds)


def map_transcripts_protein(regions_d, chromosomes_d, strands_d, genome):
    """
    Map transcript to reference protein sequence. Remove transcripts that do not pass control check.
    :param regions_d: dictionary of IntervalTrees with genomic regions for elements
    :param chromosomes_d: dict, keys are elements, values are chromosomes
    :param strands_d: dict, keys are elements, values are strands
    :param genome: str, genome to use
    :return:
    regions_d: dictionary of IntervalTrees with genomic regions for elements updated
    protein_d: dictionary of reference translated protein sequences per transcript
    """
    # TODO remove hardcoded file
    genetic_code_path = '/home/carnedo/projects/oncodriveclustl/oncodriveclustl/data/genetic_code_ncbi_20180727_v1.json'
    with open(genetic_code_path, 'rt') as fd:
        genetic_code = json.load(fd)
    # protein_d_path = os.path.join(cache, 'protein_sequences.pickle')
    reverse_d = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    start_codons = ['ATG']  #, 'TTG', 'CTG']
    protein_d = {}
    elements_to_skip = set()

    for element, regions in regions_d.items():
        cds = []
        # Get nucleotide sequence
        for interval in sorted(regions, reverse=False):
            start = interval[0]
            size = interval[1] - interval[0]  # no need +1 because interval end already +1
            sequence = bg.refseq(genome, chromosomes_d[element], start, size)
            if sequence.count('N') == 0:
                cds.extend(sequence)
            else:
                logger.warning('Found N nucleotide in {}. Element is discarded from analysis\n'.format(element))
                elements_to_skip.add(element)
                break

        # Reverse if needed
        if strands_d[element] == '-':
            # Reverse
            cds.reverse()
            # Complementary
            cds = ''.join([reverse_d.get(i, i) for i in cds])
        else:
            cds = ''.join(cds)

        # Check and translate to protein
        if len(cds) % 3 == 0 and cds[:3] in start_codons:
            protein_d[element] = ''.join([genetic_code[cds[i:i + 3]][0] for i in range(0, len(cds), 3)])
        else:
            logger.warning('{} cannot be translated to protein. Element is discarded from analysis'.format(element))
            elements_to_skip.add(element)

    # # Save to pickle
    # with open(protein_d_path, 'wb') as fd:
    #     pickle.dump(protein_d, fd, protocol=2)

    return elements_to_skip


def parse(input_regions, elements, input_mutations, cds, conseq, protein, genome):
    """Parse genomic regions and dataset of cancer type mutations
    :param input_regions: path to file containing mutational data
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param input_mutations: path tab file
    :param cds: bool, True calculates clustering on cds
    :param conseq: True, use aa consequence type
    :param protein: bool, True analyzes clustering in translated protein sequences
    :return
        regions_d: dictionary of IntervalTrees with genomic regions for elements
        cds_d: dictionary of dictionaries with relative cds position of genomic regions
        chromosomes_d: dict, keys are elements, values are chromosomes
        strands_d: dict, keys are elements, values are strands
        mutations_d: dictionary, key = element, value = list of mutations formatted as namedtuple
        samples_d: dictionary, key = sample, value = number of mutations
        cohorts_d: dictionary, key = element, value = set of cohorts with element mutations
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
        # Remove transcripts
        for element in elements_to_skip:
            del regions_d[element]
            del cds_d[element]
            mutations_d.pop(element, None)

    # set_of_regions = set()
    # for k, v in regions_d.items():
    #     print(k, v)
    #     for i in v:
    #         set_of_regions.add(i[0])
    #         set_of_regions.add(i[1])

    # print(sorted(list(set_of_regions)))
    #
    # for k, v in mutations_d.items():
    #     print(k)
    #     for m in v:
    #         check = m.region[1] > m.position >= m.region[0]
    #         if not check:
    #             print(v, m)

    return regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, cohorts_d
