# Import modules
import os
import gzip
import csv

import daiquiri
import pickle


def check_compression(file):
    """
    Check input file compression
    :param file: path to file
    :return:
            comp: str
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
        logger.critical('{} file not found'.format(file))
        quit(-1)

    return comp


def check_tabular_csv(file):
    """
    Check input file tabular or csv format
    :param file: path to file
    :return:
            read function
            read mode
            file delimiter
    """

    comp = check_compression(file)

    if comp == 'gz':
        read_function = gzip.open
        mode = 'rt'
    else:
        read_function = open
        mode = 'r'

    with read_function(file) as fd:
        for line in fd:
            if comp:
                line = line.decode()
            dialect = csv.Sniffer().sniff(line, delimiters=None)
            chr = 'CHROMOSOME' in line
            pos = 'POSITION' in line
            ref = 'REF' in line
            alt = 'ALT' in line
            sample = 'SAMPLE' in line
            cancer_type = 'CANCER_TYPE' in line
            break

        if chr == pos == ref == alt == sample == True:
            return read_function, mode, dialect.delimiter, cancer_type

        else:
            logger.critical('{} does not contain header and/or header is not in correct format'.format(file))
            quit(-1)


def check_signature(file, kmer):
    """
    Check input signatures format
    :param file:
    :param kmer:
    :return: error: True if file is not in correct format
    """
    error = False
    try:
        signature = pickle.load(open(file, "rb"))
        # Input signatures calculated for same kmer
        if len(list(signature['probabilities'].keys())[0][0]) != int(kmer):
            error = True
    except Exception:
        # Input signatures cannot be read
        error = True

    return error


