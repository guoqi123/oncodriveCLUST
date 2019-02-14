"""
This module contains functions to parse input mutations files into tabular format compatible with OncodriveCLUSTL
"""

# Import modules
import os
import gzip

import click


@click.command()
@click.option('-i', '--input_directory', required=True, type=click.Path(exists=True),
              help='Directory containing vcf files of somatic mutations')
@click.option('-o', '--output_file_tsv', required=True, help='Output tsv file')
def vcf_to_tsv(input_directory, output_file_tsv):
    """
    Parse 'vcf' or 'vcf.gz' mutation files in a directory into tabular format with the following columns:
        'CHROMOSOME'
        'POSITION'
        'REF'
        'ALT'
        'SAMPLE'
        'ID'
    This function will parse all vcf files in the directory. It is assumed that one file corresponds to one sample.
    Therefore, file name will be used as sample identifier ('SAMPLE' column). No filters are taken into account to
    parse mutations.

    Args:
        input_directory (str): path to directory containing vcf files
        output_file_tsv (str): path to output file

    Returns:
        None

    """
    output_file_log = '{}.log'.format(output_file_tsv.split('.')[0])

    with open(output_file_log, 'w') as logfd:
        with open(output_file_tsv, 'w') as ofd:
            ofd.write('{}\n'.format('\t'.join(['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID'])))

            for entry in os.scandir(input_directory):

                # Check compression
                if os.path.isfile(entry.path):
                    read_function = None
                    mode = None
                    if entry.name.endswith('.vcf.gz'):
                        read_function = gzip.open
                        mode = 'rt'
                    elif entry.name.endswith('.vcf'):
                        read_function = open
                        mode = 'r'

                        # Read
                    if read_function and mode:
                        with read_function(entry.path, mode) as fd:
                            sample_ident = entry.name.split('.')[0]
                            for line in fd:
                                if not line.startswith('#'):
                                    chrom, pos, mut_ident, ref, alt = line.strip().split('\t')[:5]
                                    ofd.write('{}\n'.format('\t'.join([chrom, pos, ref, alt, sample_ident, mut_ident])))
                        logfd.write('{}\tPARSED\n'.format(entry.path))
                    else:
                        logfd.write('{}\tPASS\n'.format(entry.path))


"""
Once you have installed oncodriveclustl, you can run this function through the command line as: 
~$ parse_vcf -i [INPUT_DIRECTORY] -o [OUTPUT_FILE]
"""
