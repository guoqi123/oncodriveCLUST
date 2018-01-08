# Import modules
import logging
import pickle
import csv
import gzip

import click
import daiquiri
from tqdm import tqdm


class Parser:
    """Class to parse the datasets with somatic mutations"""
    def __init__(self):
        self.CHROMOSOME = 'CHROMOSOME'
        self.POSITION = 'POSITION'
        self.REF = 'REF'
        self.ALT = 'ALT'
        self.SAMPLE = 'SAMPLE'
        self.TYPE = 'TYPE'
        self.CANCER_TYPE = 'CANCER_TYPE'
        self.SIGNATURE = 'SIGNATURE'
        self.TRANSCRIPT = 'TRANSCRIPT'
        self.SYMBOL = 'SYMBOL'


class Signature:
    """Class to calculate the mutational signatures of a dataset"""

    def __init__(self, kmer, genome, log_level, start_at_0=False, mutation_type='subs'):
        self.mutation_type = mutation_type
        self.kmer = kmer
        self.start_at_0 = start_at_0
        self.genome_build = self.load_genome(genome)

        fx = self.triplets if kmer == 3 else self.pentamers
        self.signatures = {'counts': fx(), 'probabilities': fx()}

    @staticmethod
    def load_genome(genome):
        """
        Load the appropriate genome build
        :param genome: str, genome
        :return: instance of the bgreference
        """
        if genome == 'hg19':
            from bgreference import hg19 as genome_build
        elif genome == 'mm10':
            from bgreference import mm10 as genome_build
        else:
            from bgreference import c3h as genome_build

        return genome_build

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

    def save(self, signature_file):
        """Save the signature to an output file"""
        with open(signature_file, 'wb') as fd:
            pickle.dump(self.signatures, fd, protocol=2)

    @staticmethod
    def load(signature_file):
        """Load precalculated signatures"""
        with open(signature_file, 'rb') as fd:
            signatures = pickle.load(fd)
        return signatures

    def calculate(self, mutations_file, gz):
        """Calculate the signature of a dataset
        :param mutations_file: path to file containing mutations
        :param gz: boolean, True if gz compressed mutations file
        :return: None
        """
        parser = Parser()

        read_function = open if gz is False else gzip.open
        mode = 'r' if gz is False else 'rt'
        if '.tab' in mutations_file:
            delimiter = '\t'
        elif '.in' in mutations_file:
            delimiter = '\t'
        elif '.txt' in mutations_file:
            delimiter = '\t'
        elif '.csv' in mutations_file:
            delimiter = ','
        else:
            logger.critical('Mutations file format not recognized. Please, provide tabular or csv formatted data')
            quit()

        with read_function(mutations_file, mode) as csvfile:
            fd = csv.DictReader(csvfile, delimiter=delimiter)
            count = 0
            for line in tqdm(fd):
                chromosome = line[parser.CHROMOSOME]
                position = int(line[parser.POSITION])
                ref = line[parser.REF]
                alt = line[parser.ALT]
                # mutation_type = line[parser.TYPE]

                # if mutation_type != 'subs':
                #     continue

                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if self.kmer == 3:
                            signature_ref = self.genome_build(chromosome, position - 1, 3).upper()
                            # Check reference nucleotide in mutations file equals reference genome nucleotide
                            if signature_ref[1] != ref:
                                logger.warning('Input REF nucleotide {} in position {} is not equal to '
                                                    'reference genome {} REF nucleotide {}'.format(
                                    ref, position, self.genome_build.__name__, signature_ref[1]
                                ))
                            signature_alt = ''.join([signature_ref[0], alt, signature_ref[-1]])
                        else:
                            signature_ref = self.genome_build(chromosome, position - 2, 5).upper()
                            signature_alt = ''.join(
                                [signature_ref[0], signature_ref[1], alt, signature_ref[-2], signature_ref[-1]]
                            )
                        # Control for N nucleotides
                        N_ref = signature_ref.count('N')
                        N_alt = signature_alt.count('N')
                        if N_ref == 0 and N_alt == 0:
                            try:
                                self.signatures['counts'][(signature_ref, signature_alt)] += 1
                                count += 1
                            except KeyError as e:
                                logger.error('{} not found in dictionary of mutations. Mutation is not taken '
                                                  'into account for signatures calculation'.format(
                                    e, signature_ref, signature_alt
                                ))
        # Calculate probabilities
        try:
            self.signatures['probabilities'] = {k: v / count for k, v in self.signatures['counts'].items()}
        except ZeroDivisionError as e:
            logger.error('{}. Impossible to calculate signatures, no substitution mutations found in {}'.format(
                e, mutations_file
            ))


@click.command()
@click.argument('input_file')
@click.argument('output_file')
@click.option('-g', '--genome', default='hg19', type=click.Choice(['hg19', 'mm10', 'c3h']), help="genome to use")
@click.option('-k', '--kmer', default='3', type=click.Choice(['3', '5']), help="number of nucleotides' signature")
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
@click.option('--start_at_0', is_flag=True)
def main(input_file, output_file, genome, kmer, start_at_0, log_level, gz):
    """Calculate the signature of a dataset"""

    global logger
    logger = daiquiri.getLogger()

    # Calculate signatures
    obj = Signature(kmer=int(kmer), genome=genome, log_level=log_level, start_at_0=start_at_0)
    obj.calculate(input_file, gz)
    obj.save(output_file)

if __name__ == '__main__':
    main()
