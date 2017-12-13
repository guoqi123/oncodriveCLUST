# Import modules
import pickle
import csv
import gzip
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

    def __init__(self, nucleotides, start_at_0=False, mutation_type='subs'):
        self.mutation_type = mutation_type
        self.nucleotides = nucleotides
        self.start_at_0 = start_at_0
        fx = self.triplets if nucleotides == 3 else self.pentamers
        self.signatures = {'counts': fx(), 'probabilities': fx()}
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

    @staticmethod
    def pentamers():
        """Create a dictionary with all the possible pentamers
        :return: dict
        """
        nucleotides = set(['A', 'C', 'T', 'G'])
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
        """Calculate the signature of a dataset"""
        parser = Parser()

        # Take into account if the mutations are 0 based or 1 based
        offset = 1 if self.start_at_0 is True else 2

        read_function = open if gz is False else gzip.open
        mode = 'r' if gz is False else 'rt'

        with read_function(mutations_file, mode) as csvfile:
            fd = csv.DictReader(csvfile, delimiter='\t')
            count = 0
            for line in tqdm(fd):
                chromosome = line[parser.CHROMOSOME]
                position = int(line[parser.POSITION])
                ref = line[parser.REF]
                alt = line[parser.ALT]
                # Read substitutions only
                if len(ref) == 1 and len(alt) == 1:
                    if ref != '-' and alt != '-':
                        if self.nucleotides == 3:
                            signature_ref = hg19(chromosome, position - 1, 3).upper()
                            signature_alt = ''.join([signature_ref[0], alt, signature_ref[-1]])
                        else:
                            signature_ref = hg19(chromosome, position - 2, 5).upper()
                            signature_alt = ''.join(
                                [signature_ref[0], signature_ref[1], alt, signature_ref[-2], signature_ref[-1]])
                        self.signatures['counts'][(signature_ref, signature_alt)] += 1
                        count += 1
                    else:
                        continue
                else:
                    continue

        self.signatures['probabilities'] = {k: v / count for k, v in self.signatures['counts'].items()}


@click.command()
@click.argument('input_file')
@click.argument('output_file')
@click.option('--nucleotides', default='3', type=click.Choice(['3', '5']), help="number of nucleotides of the signatures")
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['info', 'warning', 'error']))
@click.option('--start_at_0', is_flag=True)
def main(input_file, output_file, nucleotides, start_at_0, log_level):
    """Calculate the signature of a dataset"""
    global logger

    # Configure the colorlog module
    logger = colorlog.getLogger()

    set_logger(log_level)
    signature = Signature(nucleotides=int(nucleotides), start_at_0=start_at_0)
    signature.calculate(input_file)
    signature.save(output_file)


if __name__ == '__main__':
    main()

# run as: python signature.py inputs/mutations/pancanatlas/LUSC.txt signatures/LUSC.pickle