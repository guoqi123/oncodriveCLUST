# Import modules
import os
from collections import defaultdict
import pickle
import csv
import gzip

import click
import daiquiri
from tqdm import tqdm
import bgreference as bg

from oncodriveclustl.utils import preprocessing as prep
logger = daiquiri.getLogger()


class Parser:
    """Class to parse the datasets with somatic mutations"""
    def __init__(self):
        self.CHROMOSOME = 'CHROMOSOME'
        self.POSITION = 'POSITION'
        self.REF = 'REF'
        self.ALT = 'ALT'
        self.SAMPLE = 'SAMPLE'
        self.CANCER_TYPE = 'CANCER_TYPE'
        self.KMER = 'KMER'


class Signature:
    """Class to calculate the mutational signatures of a dataset"""

    def __init__(self, kmer, genome, log_level='info', start_at_0=False, mutation_type='subs', pancancer=False):
        self.kmer = kmer
        self.genome = genome
        self.start_at_0 = start_at_0
        self.mutation_type = mutation_type
        self.pancancer = pancancer
        fx = self.triplets if kmer == 3 else self.pentamers
        self.signatures = defaultdict(lambda: {'counts': fx(), 'probabilities': fx()})

    @staticmethod
    def triplets():
        """Create a dictionary with all the possible triplets
        :return: dict
        """
        nucleotides = {'A', 'C', 'T', 'G'}
        results = set()
        for nuc_1 in nucleotides:
            for nuc_2 in nucleotides:
                for nuc_3 in nucleotides:
                    for alt in nucleotides - {nuc_2, }:
                        results.add((''.join([nuc_1, nuc_2, nuc_3]),
                                     ''.join([nuc_1, alt, nuc_3])))
        return {key: 0 for key in results}

    @staticmethod
    def pentamers():
        """Create a dictionary with all the possible pentamers
        :return: dict
        """
        nucleotides = {'A', 'C', 'T', 'G'}
        results = set()
        for nuc_1 in nucleotides:
            for nuc_2 in nucleotides:
                for nuc_3 in nucleotides:
                    for nuc_4 in nucleotides:
                        for nuc_5 in nucleotides:
                            for alt in nucleotides - {nuc_3, }:
                                results.add((''.join([nuc_1, nuc_2, nuc_3, nuc_4, nuc_5]),
                                             ''.join([nuc_1, nuc_2, alt, nuc_4, nuc_5])))
        return {key: 0 for key in results}

    def save(self, directory, prefix):
        """Save the signature to an output file
        :param directory: path to output directory
        :param prefix: str, prefix for output file
        :return None
        """
        for cohort, values in self.signatures.items():
            if cohort == 'COHORT':
                output_file = os.path.join(directory, '{}_kmer_{}.pickle'.format(prefix, self.kmer))
            else:
                output_file = os.path.join(directory, '{}_kmer_{}.pickle'.format(cohort, self.kmer))
            with open(output_file, 'wb') as fd:
                pickle.dump(values, fd, protocol=2)

    @staticmethod
    def load(signature_file):
        """Load precalculated signatures"""
        with open(signature_file, 'rb') as fd:
            signatures = pickle.load(fd)
        return signatures

    def calculate(self, mutations_file):
        """Calculate the signature of a dataset
        :param mutations_file: path to file containing mutations
        :param pancancer: bool, True if pancancer analysis
        :return: None
        """
        parser = Parser()
        read_function, mode, delimiter, _, _ = prep.check_tabular_csv(mutations_file)

        with read_function(mutations_file, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            count = 0
            for line in tqdm(fd):
                chromosome = line[parser.CHROMOSOME]
                position = int(line[parser.POSITION])
                ref = line[parser.REF]
                alt = line[parser.ALT]
                if self.pancancer:
                    cancer_type = line[parser.CANCER_TYPE]
                    kmer = int(line[parser.KMER])
                else:
                    cancer_type = 'COHORT'
                    kmer = self.kmer

                if kmer == self.kmer:
                    # Read substitutions only
                    if len(ref) == 1 and len(alt) == 1:
                        if ref != '-' and alt != '-':
                            if self.kmer == 3:
                                signature_ref = bg.refseq(self.genome, chromosome,  position - 1, 3).upper()
                                # Check reference nucleotide in mutations file equals reference genome nucleotide
                                if signature_ref[1] != ref:
                                    logger.warning('Input REF nucleotide {} in position {} is not equal to '
                                                        'reference genome {} REF nucleotide {}'.format(
                                        ref, position, self.genome, signature_ref[1]
                                    ))
                                signature_alt = ''.join([signature_ref[0], alt, signature_ref[-1]])
                            else:
                                signature_ref = bg.refseq(self.genome, chromosome,  position - 2, 5).upper()
                                signature_alt = ''.join(
                                    [signature_ref[0], signature_ref[1], alt, signature_ref[-2], signature_ref[-1]]
                                )
                            # Control for N nucleotides
                            N_ref = signature_ref.count('N')
                            N_alt = signature_alt.count('N')
                            if N_ref == 0 and N_alt == 0:
                                try:
                                    self.signatures[cancer_type]['counts'][(signature_ref, signature_alt)] += 1
                                    count += 1
                                except KeyError as e:
                                    logger.error('{} not found in dictionary of mutations. Mutation is not taken '
                                                      'into account for signatures calculation'.format(
                                        e, signature_ref, signature_alt
                                    ))

        # Calculate probabilities
        for cohort, values in self.signatures.items():
            try:
                values['probabilities'] = {k: v / count for k, v in values['counts'].items()}
            except ZeroDivisionError as e:
                logger.error('{}. Impossible to calculate signatures, no substitution mutations found in {}'.format(
                    e, mutations_file
                ))