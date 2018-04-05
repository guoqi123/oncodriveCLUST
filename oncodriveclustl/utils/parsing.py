# Import modules
import os
import gzip
import csv
import pickle
from collections import defaultdict
from collections import namedtuple

import daiquiri
import tabix
import numpy as np
from intervaltree import IntervalTree

from oncodriveclustl.utils import preprocessing as prep


Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample')
Cds = namedtuple('Cds', 'start, end')

def dict_of_sets():
    return defaultdict(set)


def read_regions(input_regions, elements):
    """
    Parse input genomic regions
    :param input_regions: path to input genomic regions
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
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
        with gzip.open(input_regions, 'rb') as fd:
            for line in fd:
                line = line.decode()
                chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
                if elements and symbol not in elements:
                    continue
                if int(start) != int(end):
                    trees[chromosome][int(start): int(end) + 1] = symbol + '_' + ensid
                    regions_d[symbol + '_' + ensid].addi(int(start), (int(end) + 1))
                    chromosomes_d[symbol + '_' + ensid] = chromosome
                    strands_d[symbol + '_' + ensid] = strand
        if not regions_d.keys():
            logger.critical('No elements found in genomic regions. Please, check input data')
            quit()
    else:
        logger.critical('Genomic regions are not compressed, please input .gz file')
        quit()
    if len(elements) == 1 and len(regions_d) != 1:
        logger.warning('{} has more than one Ensembl id'.format(elements))

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


def read_mutations(input_mutations, trees, vep_file, conseq):
    """
    Read mutations file and intersect with regions. Only substitutions.
    :param input_mutations: path to file containing mutations
    :param trees: dictionary of dictionary of IntervalTrees containing intervals of genomic elements per chromosome
    :param vep_file: path to file containing mutations in vep format
    :param conseq: True, use aa consequence type
    :return:
        mutations_d: dictionary, key = element, value = list of mutations formatted as namedtuple
        samples_d: dictionary, key = sample, value = number of mutations
    """
    global Mutation
    mutations_d = defaultdict(list)
    samples_d = defaultdict(int)
    read_function, mode, delimiter = prep.check_tabular_csv(input_mutations)

    if vep_file:
        ident_vep = set()
        vep_moderate_high = {'protein_altering_variant', 'missense_variant', 'inframe_deletion', 'inframe_insertion',
                             'transcript_amplification', 'start_lost', 'stop_lost', 'frameshift_variant', 'stop_gained',
                             'splice_donor_variant', 'splice_acceptor_variant', 'transcript_ablation'}

        # Get ids from mutations impact moderate to high (> missense)
        with gzip.open(vep_file, 'rb') as fd:
            next(fd)
            for line in fd:
                line = line.decode()
                variation, location, _, _, _, _, conseq = line.strip().split('\t')[:7]
                if conseq in vep_moderate_high:
                    ident = variation.split('__')[0]
                    ident_vep.add(ident)
    if conseq:
        # tb = tabix.open('/workspace/datasets/phd_snp_g/input_files_cds/vep_canonical.tsv.gz')
        # dict_of_sets()
        with open('/home/carnedo/projects/inputs/vep/vep_canonical.pickle', 'rb') as fd:
            conseq_d = pickle.load(fd)

        # Read mutations
        with read_function(input_mutations, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            for line in fd:
                chromosome = line['CHROMOSOME']
                position = int(line['POSITION'])
                ref = line['REF']
                alt = line['ALT']
                sample = line['SAMPLE']
                samples_d[sample] += 1
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if trees[chromosome][int(position)] != set():
                            results = trees[chromosome][int(position)]
                            for res in results:
                                #muttype = np.random.choice([0, 1], size=1, p=[0.2, 0.8])

                                ensid = res.data.split('_')[1]
                                try:
                                    muttype = 0 if position in conseq_d[ensid][alt].values() else 1
                                except:
                                    muttype = 1

                                # consequences = [
                                #     c[8] for c in tb.query(chromosome, position - 1, position) if c[4] == alt
                                # ]
                                # if len(set(consequences)) > 1:
                                #     logger.warning('More than one consequence type in {}'.format(res.data))
                                # muttype = 0 if all([i == 'synonymous_variant' for i in consequences]) else 1

                                m = Mutation(position, (res.begin, res.end), alt, muttype, sample)
                                mutations_d[res.data].append(m)
    else:
        with read_function(input_mutations, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            for line in fd:
                chromosome = line['CHROMOSOME']
                position = int(line['POSITION'])
                ref = line['REF']
                alt = line['ALT']
                sample = line['SAMPLE']
                ident = line['ID']
                samples_d[sample] += 1
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if not vep_file:
                            if trees[chromosome][int(position)] != set():
                                results = trees[chromosome][int(position)]
                                for res in results:
                                    muttype = 1
                                    m = Mutation(position, (res.begin, res.end), alt, muttype, sample)
                                    mutations_d[res.data].append(m)
                        else:
                            if ident in ident_vep:
                                if trees[chromosome][int(position)] != set():
                                    results = trees[chromosome][int(position)]
                                    for res in results:
                                        muttype = 1
                                        m = Mutation(position, (res.begin, res.end), alt, muttype, sample)
                                        mutations_d[res.data].append(m)

    return mutations_d, samples_d


def parse(input_regions, elements, input_mutations, cds, vep_file, conseq):
    """Parse genomic regions and dataset of cancer type mutations
    :param input_regions: path to file containing mutational data
    :param elements: set, list of elements to analyze. If the set is empty all the elements will be analyzed
    :param input_mutations: path tab file
    :param cds: bool, True calculates clustering on cds
    :param vep_file: path to vep file
    :param conseq: True, use aa consequence type
    :return
        regions_d: dictionary of IntervalTrees with genomic regions for elements
        cds_d: dictionary of dictionaries with relative cds position of genomic regions
        chromosomes_d: dict, keys are elements, values are chromosomes
        strands_d: dict, keys are elements, values are strands
        mutations_d: dictionary, key = element, value = list of mutations formatted as namedtuple
        samples_d: dictionary, key = sample, value = number of mutations
    """
    global logger
    logger = daiquiri.getLogger()

    regions_d, chromosomes_d, strands_d, trees = read_regions(input_regions, elements)
    if cds:
        cds_d = map_regions_cds(regions_d)
    else:
        cds_d = {}
    logger.info('Regions parsed')
    mutations_d, samples_d = read_mutations(input_mutations, trees, vep_file, conseq)
    logger.info('Mutations parsed')

    return regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d
