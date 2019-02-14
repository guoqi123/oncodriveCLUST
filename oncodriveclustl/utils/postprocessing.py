"""
Contains functions to write OncodriveCLUSTL's results
"""
import os
import logging

import numpy as np
import pandas as pd
import daiquiri
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt

# Global variables
logs = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def mtc(p_value, alpha=0.01):
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
        results (tuple): tuple containing two dictionaries with results, keys are element's symbols
        directory (str): path to output directory
        file (str): output file name
        is_gzip (bool): True generates gzip compressed output file

    Returns:
        sorted_list_elements (list): sorted elements to rank clusters file

    """

    file = os.path.join(directory, file)
    elements_results, info = results

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

    df = pd.DataFrame(columns=header, index=[i for i in range(len(elements_results))])

    i = 0
    for element in elements_results.keys():
        sym, identif = element.split('//')
        _, chrom, strand, _, length, cgc = info[element]
        muts, muts_in_clu, obs_clu, sim_clu, obs_score, epval, apval, topcpval = elements_results[element]
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

    # Makes sure the data are float
    for p_value in ['P_ANALYTICAL', 'P_EMPIRICAL', 'P_TOPCLUSTER']:
        df[p_value] = df[p_value].astype(float)

    # Calculate q-values
    df_nonempty = df[np.isfinite(df['P_ANALYTICAL'])].copy()
    df_empty = df[~np.isfinite(df['P_ANALYTICAL'])].copy()

    if len(df_nonempty) != 0:
        for p_value in ['ANALYTICAL', 'EMPIRICAL', 'TOPCLUSTER']:
            df_nonempty['Q_' + p_value] = mtc(df_nonempty['P_' + p_value])
            df_empty['Q_' + p_value] = np.nan
    else:
        logger.critical('No clusters found in the analysis. No p-values are retrieved')
        for q_value in ['Q_ANALYTICAL', 'Q_EMPIRICAL', 'Q_TOPCLUSTER']:
            df_nonempty[q_value] = np.nan
            df_empty[q_value] = np.nan

    df = pd.concat([df_nonempty, df_empty])

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
    df.sort_values(by=['Q_ANALYTICAL', 'P_ANALYTICAL', 'SCORE', 'CGC'],
                   ascending=[True, True, False, False],
                   inplace=True
                   )
    sorted_list_elements = df['SYMBOL'].tolist()

    if is_gzip is True:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False)

    return sorted_list_elements


def write_cluster_results(genome, results, directory, file, sorter, is_gzip):
    """Save clusters results to the output file. Order according to elements' ranking

    Args:
        genome (str): reference genome
        results (tuple): tuple containing two dictionaries of results, keys are element's symbols
        directory (str): path to output directory
        file (str): output file name
        sorter (list): element symbols ranked by elements p-value to rank clusters
        is_gzip (bool): True generates gzip compressed output file

    Returns:
        None

    """
    sorter_index = dict(zip(sorter, range(len(sorter))))
    file = os.path.join(directory, file)
    clusters_results, info = results

    header = ['RANK',
              'SYMBOL',
              'ENSID',
              'CGC',
              'CHROMOSOME',
              'STRAND',
              'COORDINATES',
              'MAX_COORD',
              'WIDTH',
              'N_MUT',
              'N_SAMPLES',
              'FRA_UNIQ_SAMPLES',
              'SCORE',
              'P']

    with open(file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element in clusters_results.keys():
            sym, identif = element.split('//')
            regions, chrom, strand, _, length, cgc = info[element]
            clustersinfo, _, _ = clusters_results[element]
            if genome != 'hg19':
                cgc = 'Non Available'
            if type(clustersinfo) != float:
                rank = sorter_index[sym] + 1
                for interval in clustersinfo:
                    for c, v in interval.data.items():
                        left_coord = v['left_m'][1]
                        max_coord = v['max'][1]
                        right_coord = v['right_m'][1]
                        if regions[left_coord] == regions[right_coord]:
                            coordinates = '{},{}'.format(left_coord, right_coord)
                        else:
                            coordinates = '{},{};{},{}'.format(left_coord,
                                                               list(regions[left_coord])[0][1] - 1,  # end tree + 1
                                                               list(regions[right_coord])[0][0],
                                                               right_coord)
                        fd.write('{}\n'.format('\t'.join(map(str, [
                                    rank,
                                    sym,
                                    identif,
                                    cgc,
                                    chrom,
                                    strand,
                                    coordinates,
                                    max_coord,
                                    abs(v['right_m'][0] - v['left_m'][0] + 1),
                                    len(v['mutations']),
                                    len(v['samples']),
                                    v['fra_uniq_samples'],
                                    round(v['score'], 4),
                                    v['p']])
                        )))
    # Sort
    df = pd.read_csv(file, sep='\t', header=0, compression=None)
    df.sort_values(by=['RANK', 'P', 'SCORE'], ascending=[True, True, False], inplace=True)

    if is_gzip is True:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=file, sep='\t', na_rep='', index=False)
