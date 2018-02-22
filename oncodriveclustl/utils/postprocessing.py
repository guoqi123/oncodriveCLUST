# Import modules
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


def mtc(p_value):
    """
    Calculate q-value
    :param p_value: array
    :return: q_value
    """

    return mlpt(p_value, alpha=0.05, method='fdr_bh')[1]


def write_element_results(genome, results, directory, file, gzip):
    """Save results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param directory: str, output directory
    :param file: str, output file, if elements in elements file
    :param gzip: bool, True generates gzip compressed output file
    :return: None
    """

    global logger
    logger = daiquiri.getLogger()
    output_file = directory + '/elements_' + file + '.txt'

    # TODO: length if not cds?
    header = ['SYMBOL', 'ENSID', 'CHROMOSOME', 'STRAND', 'LENGTH', 'N_MUT', 'N_CLU',
              'SCORE',
              'P_EMPIRICAL', 'P_ANALYTICAL', 'P_TOPCLUSTER',
              'CGC']
    df = pd.DataFrame(columns=header, index=[i for i in range(len(results))])

    i = 0
    for element, values in results.items():
        sym, id = element.split('_')
        chr, strand, length, muts, obs_clusters, obs_score, epval, apval, topcpval, cgc = values
        if genome != 'hg19':
            cgc = 'Non Available'
        df.loc[i] = pd.Series({
            'SYMBOL':sym, 'ENSID':id, 'CHROMOSOME':chr, 'STRAND':strand, 'LENGTH':length,
            'N_MUT':muts, 'N_CLU':obs_clusters, 'SCORE':obs_score,
            'P_EMPIRICAL':epval, 'P_ANALYTICAL':apval, 'P_TOPCLUSTER':topcpval,'CGC':cgc})
        i += 1

    try:
        # Makes sure the data are float
        df['P_ANALYTICAL'] = df['P_ANALYTICAL'].astype(float)
        df['P_EMPIRICAL'] = df['P_EMPIRICAL'].astype(float)
        df['P_TOPCLUSTER'] = df['P_TOPCLUSTER'].astype(float)

        # Calculate q-values
        df_nonempty = df[np.isfinite(df['P_ANALYTICAL'])].copy()
        df_empty = df[~np.isfinite(df['P_ANALYTICAL'])].copy()
        df_nonempty['Q_EMPIRICAL'] = mtc(df_nonempty['P_EMPIRICAL'])
        df_nonempty['Q_ANALYTICAL'] = mtc(df_nonempty['P_ANALYTICAL'])
        df_nonempty['Q_TOPCLUSTER'] = mtc(df_nonempty['P_TOPCLUSTER'])
        df_empty['Q_EMPIRICAL'] = np.nan
        df_empty['Q_ANALYTICAL'] = np.nan
        df_empty['Q_TOPCLUSTER'] = np.nan
        df = pd.concat([df_nonempty, df_empty])

        # Reorder columns
        df = df[['SYMBOL', 'ENSID', 'CGC', 'CHROMOSOME', 'STRAND', 'LENGTH', 'N_MUT', 'N_CLU', 'SCORE',
                 'P_EMPIRICAL', 'Q_EMPIRICAL','P_ANALYTICAL', 'Q_ANALYTICAL','P_TOPCLUSTER', 'Q_TOPCLUSTER']]

        # Sort by analytical p-value
        df.sort_values(by=['Q_ANALYTICAL', 'P_ANALYTICAL', 'SCORE', 'CGC'],
                       ascending=[True, True, False, False], inplace=True)

    except Exception as e:
        logger.error('{} in {}. Impossible to calculate q-values'.format(e, file))

    # Create a sorted list of elements to order the clusters file
    sorted_list_elements = df['SYMBOL'].tolist()

    if gzip is True:
        output_file = output_file + '.gz'
        df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)

    return sorted_list_elements


def write_cluster_results(genome, results, directory, file, sorter, gzip):
    """Save results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param directory: str, output directory
    :param file: str, output file, if elements in elements file
    :param sorter: list, element symbols ranked by elements p-value
    :param gzip: bool, True generates gzip compressed output file
    :return: None
    """

    sorterindex = dict(zip(sorter, range(len(sorter))))
    output_file = directory + '/clusters_' + file + '.tsv'
    header = ['RANK', 'SYMBOL', 'ENSID', 'CGC', 'CHROMOSOME', 'STRAND', 'REGION',
              '5_COORD', 'MAX_COORD', '3_COORD', 'WIDTH', 'N_MUT', 'SCORE', 'P']

    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element, values in results.items():
            sym, id = element.split('_')
            clustersinfo, chr, strand, cgc= values
            if genome != 'hg19':
                cgc = 'Non Available'
            if type(clustersinfo) != float:
                rank = sorterindex[sym] + 1
                for interval in clustersinfo:
                    for c, v in interval.data.items():
                        fd.write('{}\t{}\t{}\t{}\t{}\t{}\t[{},{}]\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            rank, sym, id, cgc, chr, strand, interval[0], interval[1],
                            v['left_m'][0]+interval[0], v['max'][0]+interval[0], v['right_m'][0]+interval[0],
                            abs(v['right_m'][0] - v['left_m'][0]),
                            sum(v['mutations'].values()), v['score'], v['p']))
    # Sort
    # TODO change to avoid writing not compressed file
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['RANK', 'P', 'SCORE'], ascending=[True, True, False], inplace=True)

    if gzip is True:
        output_file = output_file + '.gz'
        df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False, compression='gzip')
    else:
        df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)


def write_info(input_file,
               output_directory,
               regions_file,
               genome,
               elements_file,
               elements,
               element_mutations,
               cluster_mutations,
               cds,
               smooth_window,
               cluster_window,
               cluster_score,
               element_score,
               kmer,
               n_simulations,
               simulation_mode,
               simulation_window,
               cores,
               seed,
               log_level,
               gzip):
    """Write info to output
        :param input_file: input file
        :param output_directory: output directory path
        :param regions_file: path input genomic regions, tab file
        :param genome: genome to use
        :param elements_file: file containing one element per row
        :param elements: element symbol or file containing elements
        :param element_mutations: int, cutoff of element mutations
        :param cluster_mutations: int, cutoff of cluster mutations
        :param cds: bool, True calculates clustering on cds
        :param smooth_window: int, smoothing window
        :param cluster_window: int, clustering window
        :param cluster_score: cluster score method
        :param element_score: element score method
        :param kmer: int, number of nucleotides of the signature
        :param n_simulations: int, number of simulations
        :param simulation_mode: str, simulation mode
        :param simulation_window: int, window to simulate mutations in hotspot mode
        :param cores: int, number of CPUs to use
        :param seed: int, seed
        :param log_level: verbosity of the logger
        :param gzip: bool, True generates gzip compressed output file
        :return: None
    """

    info_file = output_directory + '/' + output_directory.split('/')[-1] + '.info'

    if not os.path.isfile(info_file):
        with open(info_file, 'w') as fd:
            fd.write('input_file: {}\noutput_directory: {}\nregions_file: {}\ngenome: {}\nelements_file: {}\n'
                     'elements: {}\nelement_mutations: {}\ncluster_mutations: {}\ncds: {}\nsmooth_window: {}\n'
                     'cluster_window: {}\ncluster_score: {}\nelement_score: {}\n'
                     'kmer: {}\nn_simulations: {}\n'
                     'simulation_mode: {}\nsimulation_window: {}\ncores: {}\nseed: {}\n'
                     'log_level: {}\ngzip: {}\n'.format(input_file, output_directory, regions_file, genome,
                                              elements_file, elements, element_mutations, cluster_mutations,
                                              cds, smooth_window, cluster_window,
                                              cluster_score, element_score, kmer, n_simulations, simulation_mode,
                                              simulation_window,
                                              cores, seed, log_level, gzip))
    else:
        pass
