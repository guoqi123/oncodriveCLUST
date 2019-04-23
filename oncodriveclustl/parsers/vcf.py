"""
This module contains functions to parse input mutations files into tabular format compatible with OncodriveCLUSTL
"""

# Import modules
import os
import gzip
import logging

import click
import daiquiri

@click.command()
@click.option('-i', '--input_path', required=True, type=click.Path(exists=True),
              help='VCF File or directory containing VCF files of somatic mutations')
@click.option('-o', '--output_file', required=True, help='Output TSV file')
def vcf_to_tsv(input_path, output_file):
    """
    Parse 'vcf' or 'vcf.gz' mutation files into tabular format with the following columns:
        'CHROMOSOME'
        'POSITION'
        'REF'
        'ALT'
        'SAMPLE'
        'ID'
    This function will parse a single file or all vcf files in a directory.
    It is assumed that one file corresponds to one sample. Therefore, file name will be used as sample identifier
    ('SAMPLE' column). No filters are taken into account to parse mutations.

    Args:
        input_path (str): path to file or directory containing vcf files
        output_file (str): path to output file

    Returns:
        None

    """

    daiquiri.setup(level=logging.INFO)
    logger = daiquiri.getLogger()

    # File(s) to reformat
    files = set()
    if os.path.isfile(input_path):
        files.add(input_path)
    elif os.path.isdir(input_path):
        for entry in os.scandir(input_path):
            if 'vcf' in entry.name:
                files.add(entry.path)
    # Reformat
    with open(output_file, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID'])))
        for file in files:
            if file.endswith('.vcf.gz'):
                read_function, mode = gzip.open, 'rt'
            else:
                read_function, mode = open, 'r'

            with read_function(file, mode) as fd:
                sample_ident = file.split('/')[-1].split('.')[0]
                for line in fd:
                    if not line.startswith('#'):
                        chrom, pos, mut_ident, ref, alt = line.strip().split('\t')[:5]
                        ofd.write('{}\n'.format('\t'.join([chrom, pos, ref, alt, sample_ident, mut_ident])))
            logger.info('{} PARSED'.format(file))


"""
Once you have installed OncodriveCLUSTL, you can run this function through the command line as: 
~$ parse_vcf -i [INPUT_DIRECTORY] -o [OUTPUT_FILE]
"""
