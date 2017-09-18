# Import modules
import pickle
import csv

import colorlog
import click
from tqdm import tqdm

from bgreference import hg19


# Configure the colorlog module
# logger = colorlog.getLogger()


def set_logger(level):
    global logger
    d = {'info': colorlog.colorlog.logging.INFO,
         'warning': colorlog.colorlog.logging.WARNING,
         'error': colorlog.colorlog.logging.ERROR}
    logger.setLevel(d[level])
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter())
    logger.addHandler(handler)


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

    def __init__(self, start_at_0=False, mutation_type='subs'):
        self.mutation_type = mutation_type
        self.start_at_0 = start_at_0
        self.signatures = {'counts': self.triplets(), 'probabilities': self.triplets()}

    @staticmethod
    def triplets():
        """Create a dictionary with all the possible triplets
        :return: dict
        """
        nucleotides = set(['A', 'C', 'T', 'G'])
        results = set()
        for nuc_1 in nucleotides:
            for nuc_2 in nucleotides:
                for nuc_3 in nucleotides:
                    for alt in nucleotides - {nuc_2, }:
                        results.add((''.join([nuc_1, nuc_2, nuc_3]),
                                     ''.join([nuc_1, alt, nuc_3])))
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

    def calculate(self, mutations_file):
        """Calculate the signature of a dataset"""
        parser = Parser()

        # Take into account if the mutations are 0 based or 1 based
        offset = 1 if self.start_at_0 is True else 2

        with open(mutations_file, 'r') as csvfile:
            fd = csv.DictReader(csvfile, delimiter='\t')
            count = 0
            for line in tqdm(fd):
                chromosome = line[parser.CHROMOSOME]
                position = int(line[parser.POSITION])
                ref = line[parser.REF]
                alt = line[parser.ALT]
                mutation_type = line[parser.TYPE]
                if mutation_type != 'subs':
                    continue
                signature_ref = hg19(chromosome, position - 1, 3).upper()
                signature_alt = ''.join([signature_ref[0], alt, signature_ref[-1]])
                self.signatures['counts'][(signature_ref, signature_alt)] += 1
                count += 1
        self.signatures['probabilities'] = {k: v / count for k, v in self.signatures['counts'].items()}


@click.command()
@click.argument('input_file')
@click.argument('output_file')
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['info', 'warning', 'error']))
@click.option('--start_at_0', is_flag=True)
def main(input_file, output_file, start_at_0, log_level):
    """Calculate the signature of a dataset"""
    global logger

    # Configure the colorlog module
    logger = colorlog.getLogger()

    set_logger(log_level)
    signature = Signature(start_at_0=start_at_0)
    signature.calculate(input_file)
    signature.save(output_file)


if __name__ == '__main__':
    main()

# run as: python signature.py inputs/mutations/pancanatlas/LUSC.txt signatures/LUSC.pickle