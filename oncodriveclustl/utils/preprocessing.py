# Import modules
import os
import gzip
import csv

import daiquiri



def check_compression(file):
    """
    Check input file compression
    :param file: path to file
    :return:
            comp: str
    """

    global logger
    logger = daiquiri.getLogger()

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
                        comp = 'None'
                        break
            except Exception as e:
                logger.critical('{}. Incorrect file format for {}'.format(e, file))
                quit()
    else:
        logger.critical('{} file not found'.format(file))
        quit()

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
            if comp != 'None':
                line = line.decode()
            dialect = csv.Sniffer().sniff(line, delimiters=None)
            chr = 'CHROMOSOME' in line
            pos = 'POSITION' in line
            ref = 'REF' in line
            alt = 'ALT' in line
            break

    if chr == pos == ref == alt == True:
        header = True
    else:
        logger.critical('{} does not contain header and/or header not in correct format'.format(file))
        quit()

    return read_function, mode, dialect.delimiter


# TODO: check non overlapping regions