"""
Contains functions to check input files' format
"""

import os
import gzip
import csv

import daiquiri
import pickle

from oncodriveclustl.utils import exceptions as excep


def check_compression(file):
    """
    Check input file compression

    Args:
        file (str): path to file

    Returns:
        comp (str): compression type

    """

    global logger
    logger = daiquiri.getLogger()

    comp = None
    if os.path.isfile(file):
        try:
            with gzip.open(file, 'rb') as fd:
                for line in fd:
                    comp = 'gz'
                    break
        except OSError:
            try:
                with open(file, 'r') as fd:
                    for line in fd:
                        comp = None
                        break
            except Exception as e:
                logger.critical('{}. Incorrect file format for {}'.format(e, file))
                quit(-1)
    else:
        raise FileNotFoundError('{} file not found'.format(file))

    return comp


def check_tabular_csv(file):
    """
    Check input file tabular or csv format

    Args:
        file (str): path to file

    Returns:
        read_function: function to read file
        mode (str): reading mode
        dialect.delimiter (str): file column delimiter
        cancer_type (bool): True if `CANCER_TYPE` column exists in input file

    """

    chrom = pos = ref = alt = sample = False
    comp = check_compression(file)

    if comp == 'gz':
        read_function = gzip.open
        mode = 'rt'
    else:
        read_function = open
        mode = 'r'

    with read_function(file, mode) as fd:
        for line in fd:
            dialect = csv.Sniffer().sniff(line, delimiters=None)
            chrom = 'CHROMOSOME' in line
            pos = 'POSITION' in line
            ref = 'REF' in line
            alt = 'ALT' in line
            sample = 'SAMPLE' in line
            cancer_type = 'CANCER_TYPE' in line
            break
        if chrom and pos and ref and alt and sample:
            return read_function, mode, dialect.delimiter, cancer_type
        else:
            raise excep.UserInputError('{} does not contain header and/or header is not in correct format'.format(file))


def check_signature(file, kmer):
    """
    Check input signatures format

    Args:
        file (str): path to signatures file
        kmer (int): integer indicanting trinucleotide (3) or pentanucleotide (5) signatures

    Returns:
        error (bool): True if file is not in correct format

    """
    error = False
    signature = None
    try:
        signature = pickle.load(open(file, "rb"))
    except TypeError:
        error = True
    if signature and len(list(signature['probabilities'].keys())[0][0]) != int(kmer):
        error = True

    return error


