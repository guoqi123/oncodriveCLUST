# Import modules
import gzip
import zipfile
import csv

import daiquiri



def check_compression(file):
    """
    Check input file compression (gz or zip)
    :param file: path to file
    :return:
            comp: str
    """

    global logger
    logger = daiquiri.getLogger()

    # zip_file = zipfile.is_zipfile(file)
    # if zip_file is True:
    #     comp = 'zip'

    try:
        with gzip.open(file, 'rb') as fd:
            for line in fd:
                line = line.decode()
                # print(line)
                comp = 'gz'
                break
    except OSError:
        try:
            with open(file, 'r') as fd:
                for line in fd:
                    # print(line)
                    comp = 'None'
                    break
        except Exception as e:
            logger.critical('{}. Incorrect file format for {}'.format(e, file))
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
    # elif comp == 'zip':
    #     read_function = zipfile.ZipFile
    #     mode = 'rt'
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
        logger.critical('{} does not contain header'.format(file))
        quit()

    return read_function, mode, dialect.delimiter

#
# file_txt = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC.txt'
# file_csv = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC_csv.csv'
# file_xls = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC_xls.xlsx'
# file_gz = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC.txt.gz'
# file_zip = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC.txt.zip'
# file_targz = '/home/carnedo/projects/inputs/mutations/pancanatlas/ACC.txt.tar.gz'
# intogen = '/home/carnedo/projects/intogen/intogen_20171111/oncodriveclustl/PCATLAS_WXS_BRCA.nodups.in.gz'
#
# for f in [file_txt, file_csv, file_gz, file_targz, intogen]:
#     print(f)
#     check_tabular_csv(f)