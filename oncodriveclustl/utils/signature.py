"""
Contains the classes and functions to calculate the context based mutational probabilities of a cohort
"""

import os
from collections import defaultdict
import pickle
import csv
import itertools
import gzip  # TODO check conditional import

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


class Signature:
    """Class to calculate the mutational signatures of a dataset"""

    def __init__(self, kmer, genome, start_at_0=False, mutation_type='subs', pancancer=False):
        self.kmer = kmer
        self.genome = genome
        self.start_at_0 = start_at_0
        self.mutation_type = mutation_type
        self.pancancer = pancancer
        self.signatures = defaultdict(lambda: {'counts': self.generate_kmers(kmer),
                                               'probabilities': self.generate_kmers(kmer)})

    @staticmethod
    def generate_kmers(kmer):
        """Create a dictionary with all the possible kmers (trinucleotides or pentanucleotides) and alternates

        Args:
            kmer (int): kmer nucleotides to calculate the signature (3 or 5)

        Returns:
            dict: keys are 'ref kmer --> alt kmer', values are 0
        """
        half_kmer = kmer // 2
        nucleotides = 'ACTG'
        results = set()

        for permutation in itertools.product(nucleotides, repeat=kmer):
            kmer_nucleotide = ''.join(permutation)
            for alt in set(nucleotides).difference(kmer_nucleotide[half_kmer]):
                results.add(
                    (kmer_nucleotide, ''.join([kmer_nucleotide[0:half_kmer], alt, kmer_nucleotide[half_kmer + 1: ]]))
                )

        return {key: 0 for key in results}

    def save(self, directory, prefix):
        """Save the signature to an output file

        Args:
            directory (str): path to output directory
            prefix (str): prefix for output file

        Returns:
            None
        """
        for cohort, values in self.signatures.items():
            if cohort == 'COHORT':
                output_file = os.path.join(directory, '{}_kmer_{}.pickle'.format(prefix, self.kmer))
            else:
                output_file = os.path.join(directory, '{}_kmer_{}.pickle'.format(cohort, self.kmer))
            with open(output_file, 'wb') as fd:
                pickle.dump(values, fd, protocol=2)

    def calculate(self, mutations_file):
        """Calculate the signature of a dataset

        Args:
            mutations_file (str): path to file containing mutations
        Returns:
            None
        """
        parser = Parser()

        read_function, mode, delimiter, _ = prep.check_tabular_csv(mutations_file)

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
                else:
                    cancer_type = 'COHORT'
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if self.kmer == 3:
                            signature_ref = bg.refseq(self.genome, chromosome,  position - 1, 3).upper()
                            # Check reference nucleotide in mutations file equals reference genome nucleotide
                            if signature_ref[1] != ref:
                                logger.warning('Input REF nucleotide {} in position {} is not equal to '
                                               'reference genome {} REF nucleotide {}'.format(
                                                    ref, position, self.genome, signature_ref[1]))
                            signature_alt = ''.join([signature_ref[0], alt, signature_ref[-1]])
                        else:
                            signature_ref = bg.refseq(self.genome, chromosome,  position - 2, 5).upper()
                            signature_alt = ''.join(
                                [signature_ref[0], signature_ref[1], alt, signature_ref[-2], signature_ref[-1]]
                            )
                        # Control for N nucleotides
                        n_ref = signature_ref.count('N')
                        n_alt = signature_alt.count('N')
                        if n_ref == 0 and n_alt == 0:
                            try:
                                self.signatures[cancer_type]['counts'][(signature_ref, signature_alt)] += 1
                                count += 1
                            except KeyError as e:
                                print(signature_ref, signature_alt)
                                print(self.signatures[cancer_type]['counts'][('TCA', 'TGC')])
                                quit()
                                logger.error('{} not found in dictionary of mutations. Mutation is not taken '
                                             'into account for signatures calculation'.format(
                                                e, signature_ref, signature_alt))
        # Calculate probabilities
        for cohort, values in self.signatures.items():
            try:
                values['probabilities'] = {k: v / count for k, v in values['counts'].items()}
            except ZeroDivisionError as e:
                logger.error('{}. Impossible to calculate signatures, no substitution mutations found in {}'.format(
                    e, mutations_file
                ))
