"""
Contains functions to write OncodriveCLUSTL's results
"""
import os
import logging
import csv
import gzip
import shutil
from intervaltree import IntervalTree
from collections import namedtuple

import numpy as np
import pandas as pd
import daiquiri
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt
import bgdata as bgd
import pickle

from oncodriveclustl.utils import preprocessing as prep

# Global variables
logs = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def mtc(p_value, alpha=0.05):
    """
    Adjust p-valus by multiple test correction using false discovery rate (FDR) Benjamini–Hochberg method

    Args:
        p_value (array): list of p-values
        alpha (float): error rate

    Returns:
        (array): list of p-values corrected by multiple test

    """

    return mlpt(p_value, alpha=alpha, method='fdr_bh')[1]


def write_element_results(genome, results, directory, file, is_gzip):
    """
    Save elements results to the output file

    Args:
        genome (str): reference genome
        results (dict): dictionary of results, keys are element's symbols
        directory (str): path to output directory
        file (str): output file name
        is_gzip (bool): True generates gzip compressed output file

    Returns:
        None

    """

    file = os.path.join(directory, file)

    global logger
    logger = daiquiri.getLogger()

    header = ['SYMBOL',
              'ENSID',
              'CHROMOSOME',
              'STRAND',
              'LENGTH',
              'TOTAL_MUT',
              'CLUSTERED_MUT',
              'CLUSTERS',
              'SIM_CLUSTERS',
              'SCORE',
              'P_EMPIRICAL',
              'P_ANALYTICAL',
              'P_TOPCLUSTER',
              'CGC']

    df = pd.DataFrame(columns=header, index=[i for i in range(len(results))])

    i = 0
    for element, values in results.items():
        sym, identif = element.split('//')
        chrom, strand, length, muts, muts_in_clu, obs_clu, sim_clu, obs_score, epval, apval, topcpval, cgc = values
        if genome != 'hg19':
            cgc = 'Non Available'
        df.loc[i] = pd.Series({
            'SYMBOL': sym,
            'ENSID': identif,
            'CHROMOSOME': chrom,
            'STRAND': strand,
            'LENGTH': length,
            'TOTAL_MUT': muts,
            'CLUSTERED_MUT': muts_in_clu,
            'CLUSTERS': obs_clu,
            'SIM_CLUSTERS': sim_clu,
            'SCORE': obs_score,
            'P_EMPIRICAL': epval,
            'P_ANALYTICAL': apval,
            'P_TOPCLUSTER': topcpval,
            'CGC': cgc})
        i += 1

    try:
        # Makes sure the data are float
        for p_value in ['P_ANALYTICAL', 'P_EMPIRICAL', 'P_TOPCLUSTER']:
            df[p_value] = df[p_value].astype(float)
        # Calculate q-values
        df_nonempty = df[np.isfinite(df['P_ANALYTICAL'])].copy()
        df_empty = df[~np.isfinite(df['P_ANALYTICAL'])].copy()
        # Add to df
        for p_value in ['ANALYTICAL', 'EMPIRICAL', 'TOPCLUSTER']:
            df_nonempty['Q_' + p_value] = mtc(df_nonempty['P_' + p_value])
            df_empty['Q_' + p_value] = np.nan
        df = pd.concat([df_nonempty, df_empty])
        # Sort by analytical q-value
        df.sort_values(by=['Q_ANALYTICAL', 'P_ANALYTICAL', 'SCORE', 'CGC'],
                       ascending=[True, True, False, False], inplace=True)

    except Exception as e:
        logger.error('{} in {}. Impossible to calculate q-values'.format(e, file))
        for q_value in ['Q_ANALYTICAL', 'Q_EMPIRICAL', 'Q_TOPCLUSTER']:
            df[q_value] = np.nan

    # Reorder columns
    df = df[['SYMBOL',
             'ENSID',
             'CGC',
             'CHROMOSOME',
             'STRAND',
             'LENGTH',
             'TOTAL_MUT',
             'CLUSTERED_MUT',
             'CLUSTERS',
             'SIM_CLUSTERS',
             'SCORE',
             'P_EMPIRICAL',
             'Q_EMPIRICAL',
             'P_ANALYTICAL',
             'Q_ANALYTICAL',
             'P_TOPCLUSTER',
             'Q_TOPCLUSTER']]

    # Create a sorted list of elements to order the clusters file
    sorted_list_elements = df['SYMBOL'].tolist()

    if is_gzip is True:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False)

    return sorted_list_elements


def write_cluster_results(genome, results, directory, file, sorter, is_gzip, cds_d, protein):
    """Save clusters results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param directory: str, output directory
    :param file: str, output file, if elements in elements file
    :param sorter: list, element symbols ranked by elements p-value
    :param is_gzip: bool, True generates gzip compressed output file
    :param cds_d: dictionary of dictionaries with relative cds position of genomic regions if cds is True
    :param protein: bool, True reverses clusters positions in protein if gene strand is negative
    :return: None
    """
    reverse_cds_d = IntervalTree()
    sorter_index = dict(zip(sorter, range(len(sorter))))
    file = os.path.join(directory, file)

    header = ['RANK',
              'SYMBOL',
              'ENSID',
              'CGC',
              'CHROMOSOME',
              'STRAND',
              '5_COORD',
              'MAX_COORD',
              '3_COORD',
              'WIDTH',
              'N_MUT',
              'N_SAMPLES',
              'FRA_UNIQ_SAMPLES',
              'SCORE',
              'P']

    with open(file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element, values in results.items():
            if cds_d and not protein:
                for genomic, cds in cds_d[element].items():
                    reverse_cds_d.addi(cds[0], cds[1] + 1, genomic)   # end + 1
            sym, identif = element.split('//')
            clustersinfo, chrom, strand, length, cgc = values
            if genome != 'hg19':
                cgc = 'Non Available'
            if type(clustersinfo) != float:
                rank = sorter_index[sym] + 1
                for interval in clustersinfo:
                    for c, v in interval.data.items():
                        left_m = v['left_m'][1]
                        max_cluster = v['max'][1]
                        right_m = v['right_m'][1]

                        # # FIXME problems here!
                        # if cds_d and not protein:
                        #     for i in reverse_cds_d[v['left_m'][0]]:
                        #         start_l = i.data
                        #         end_l = i.data + (i[1] - i[0])
                        #     for i in reverse_cds_d[v['right_m'][0]]:
                        #         start_r = i.data
                        #         end_r = i.data + (i[1] - i[0])
                        #
                        #     if start_l != start_r:
                        #         region_start = (start_l, end_l)
                        #         region_end = (start_r, end_r)
                        #     else:
                        #         region_start = start_l
                        #         region_end = end_l
                        # else:
                        #     region_start = interval[0]
                        #     region_end = interval[1]

                        if protein and strand == '-':
                            left_m = length - left_m
                            max_cluster = length - max_cluster
                            right_m = length - right_m

                        fd.write('{}\n'.format('\t'.join([
                                    rank,
                                    sym,
                                    identif,
                                    cgc,
                                    chrom,
                                    strand,
                                    left_m,
                                    max_cluster,
                                    right_m,
                                    abs(v['right_m'][0] - v['left_m'][0] + 1),
                                    len(v['mutations']),
                                    len(v['samples']),
                                    v['fra_uniq_samples'],
                                    v['score'],
                                    v['p']]
                        )))
    # Sort
    df = pd.read_csv(file, sep='\t', header=0, compression=None)
    df.sort_values(by=['RANK', 'P', 'SCORE'], ascending=[True, True, False], inplace=True)

    if is_gzip is True:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False)


def write_oncohortdrive_results(mutations_file, directory, file, regions_d, vep):
    """
    Generate compressed file with input BED + mutations mapped to clusters (score, p-value, significance)
    :param mutations_file: input mutations file
    :param directory: str, output directory
    :param file: str, clusters output file to read
    :param regions_d: dict, dictionary of IntervalTrees containing genomic regions from all analyzed elements
    :param vep: bool, True considers only non-synonymous mutations
    :return: None
    """
    clusters_tree = IntervalTree()
    Cluster = namedtuple('Cluster', 'sym, ensid, score, p, sig')
    clusters_file = os.path.join(directory, file)
    output_file = os.path.join(directory, 'oncohortdrive_results.out')  # TODO: if elements file, this will overwrite
    conseq_path = bgd.get_path('oncodriveclustl', 'vep88', 'hg19_canonical_conseq')

    # Read clusters file
    if 'gz' in clusters_file:
        read_function = gzip.open
        mode = 'rt'
    else:
        read_function = open
        mode = 'r'

    with read_function(clusters_file, mode) as cf:
        next(cf)
        for line in cf:
            _, sym, ensid, _, _, _, _, clust_l, _, clust_r, _, _, _, _, score, p = line.strip().split('\t')
            element = sym + '_' + ensid
            left_coord = int(clust_l)
            right_coord = int(clust_r)
            sig = 1 if float(p) < 0.05 else 0
            if set(regions_d[element][left_coord]) != set(regions_d[element][right_coord]):
                for interval in regions_d[element][left_coord]:
                    clusters_tree.addi(left_coord, interval[1] + 1, Cluster(sym, ensid, score, p, sig))
                for interval in regions_d[element][right_coord]:
                    clusters_tree.addi(interval[0], right_coord + 1, Cluster(sym, ensid, score, p, sig))
            else:
                clusters_tree.addi(left_coord, right_coord + 1, Cluster(sym, ensid, score, p, sig))

    # Generate output file
    read_function, mode, delimiter, _ = prep.check_tabular_csv(mutations_file)

    with open(output_file, 'w') as of:
        header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'SYM', 'SYM_ENSID', 'SCORE', 'PVALUE', 'SIG_0.05']
        of.write('{}\n'.format('\t'.join(header)))
        with read_function(mutations_file, mode) as read_file:
            fd = csv.DictReader(read_file, delimiter=delimiter)
            # If not vep, write all mutations inside clusters
            if not vep:
                for line in fd:
                    chro = line['CHROMOSOME']
                    pos = int(line['POSITION'])
                    ref = line['REF']
                    alt = line['ALT']
                    sam = line['SAMPLE']

                    # Check substitutions only
                    if len(ref) == 1 and len(alt) == 1:
                        if ref != '-' and alt != '-':
                            if clusters_tree[pos]:
                                for c in clusters_tree[pos]:
                                    of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                        chro, pos, ref, alt, sam,
                                        c.data.sym, c.data.ensid, c.data.score, c.data.p, c.data.sig))
                            else:
                                of.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{5}\t{5}\t{5}\t{5}\n'.format(
                                    chro, pos, ref, alt, sam, float('nan')))
            # If vep, write only non-synonymous mutations
            else:
                for line in fd:
                    chro = line['CHROMOSOME']
                    pos = int(line['POSITION'])
                    ref = line['REF']
                    alt = line['ALT']
                    sam = line['SAMPLE']

                    # Check substitutions only
                    if len(ref) == 1 and len(alt) == 1:
                        if ref != '-' and alt != '-':
                            if clusters_tree[pos]:
                                for c in clusters_tree[pos]:
                                    path_to_vep_pickle = conseq_path + '/{}.pickle'.format(c.data.ensid)
                                    try:
                                        with open(path_to_vep_pickle, 'rb') as fd:
                                            conseq_d = pickle.load(fd)
                                            muttype = 0 if pos in conseq_d.get(alt, []) else 1
                                    except FileNotFoundError:
                                        muttype = 1
                                    if muttype == 1:
                                        of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                            chro, pos, ref, alt, sam,
                                            c.data.sym, c.data.ensid, c.data.score, c.data.p, c.data.sig))
                                    else:
                                        of.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{5}\t{5}\t{5}\t{5}\n'.format(
                                            chro, pos, ref, alt, sam, float('nan')))
                            else:
                                of.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{5}\t{5}\t{5}\t{5}\n'.format(
                                    chro, pos, ref, alt, sam, float('nan')))

    output_file_gz = output_file + '.gz'
    with open(output_file, 'rb') as of:
        with gzip.open(output_file_gz, 'wb') as ofgz:
            shutil.copyfileobj(of, ofgz)
    # Remove not compressed file
    os.remove(output_file)
