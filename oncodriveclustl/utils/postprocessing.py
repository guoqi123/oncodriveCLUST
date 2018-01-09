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


def write_element_results(genome, results, directory, file, qvalue, gzip):
    """Save results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param directory: str, output directory
    :param file: str, output file, if elements in elements file
    :param qvalue: bool, True calculates empirical and analytical q-values
    :param gzip: bool, True generates gzip compressed output file
    :return: None
    """

    global logger
    logger = daiquiri.getLogger()

    header = ['SYM', 'LEN', 'N_MUT',
              'CLU', 'MEAN_SIM_CLU', 'MEDIAN_SIM_CLU', 'SD_SIM_CLU',
              'SCORE_OBS', 'MEAN_SIM_SCORE', 'MEDIAN_SIM_SCORE', 'SD_SIM_SCORE',
              'E_PVAL', 'A_PVAL', 'CGC']

    output_file = directory + '/elements_' + file + '.txt'
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for gene_name, values in results.items():
            length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, \
            mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc = values
            if genome != 'hg19':
                cgc = 'Non Available'
            if obs_clusters and type(obs_clusters) != float:
                fd.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\t{}\t{}\n'.format(
                    gene_name, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                    obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))
            else:
                fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gene_name, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                    obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))

    # Sort by empirical p-value
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['E_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)

    # Create a sorted list of elements to order the clusters file
    sorted_list_elements = df['SYM'].tolist()

    if qvalue is True:

        try:
            # Calculate empirical q-value
            df_nonempty = df[np.isfinite(df['E_PVAL'])]
            df['E_QVAL'] = pd.DataFrame(mtc(df_nonempty.loc[:, 'E_PVAL']))

            # Sort by analytical p-value and calculate analytical q-value
            df.sort_values(by=['A_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
            df['A_QVAL'] = pd.DataFrame(mtc(df_nonempty.loc[:, 'A_PVAL']))

            # Reorder columns
            df = df[['SYM', 'LEN', 'N_MUT', 'CLU', 'MEAN_SIM_CLU', 'MEDIAN_SIM_CLU', 'SD_SIM_CLU',
                     'SCORE_OBS', 'MEAN_SIM_SCORE', 'MEDIAN_SIM_SCORE', 'SD_SIM_SCORE',
                     'E_PVAL', 'E_QVAL', 'A_PVAL', 'A_QVAL', 'CGC']]
        except Exception as e:
            logger.error('{} in {}. Impossible to calculate q-values'.format(e, file))

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
    header = ['RANK', 'SYMBOL', 'CGC', 'CLUSTER', 'N', 'LEFT_M', 'MAX', 'RIGHT_M', 'WIDTH', 'MUT', 'SCORE']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element, values in results.items():
            clustersinfo, cgc = values
            if genome != 'hg19':
                cgc = 'Non Available'
            if clustersinfo and type(clustersinfo) != float:
                rank = sorterindex[element] + 1
                for c, v in clustersinfo.items():
                    fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rank, element, cgc, v['mode'], c, v['left_m'][0], v['max'][0], v['right_m'][0],
                        abs(v['right_m'][0] - v['left_m'][0]), v['n_mutations'], v['score']))
    # Sort
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['RANK', 'SCORE'], ascending=[True, False], inplace=True)

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
               qvalue,
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
        :param qvalue: bool, True calculates empirical and analytical q-values
        :param gzip: bool, True generates gzip compressed output file
        :return: None
    """

    info_file = output_directory + '/' + output_directory.split('/')[-1] + '.info'

    if not os.path.isfile(info_file):
        with open(info_file, 'w') as fd:
            fd.write('input_file: {}\noutput_directory: {}\nregions_file: {}\ngenome: {}\nelements_file: {}\n'
                     'elements: {}\nelement_mutations: {}\ncluster_mutations: {}\nsmooth_window: {}\n'
                     'cluster_window: {}\ncluster_score: {}\nelement_score: {}\n'
                     'kmer: {}\nn_simulations: {}\n'
                     'simulation_mode: {}\nsimulation_window: {}\ncores: {}\nseed: {}\n'
                     'log_level: {}\nq-value: {}\ngzip: {}\n'.format(input_file, output_directory, regions_file, genome,
                                              elements_file, elements, element_mutations, cluster_mutations,
                                              smooth_window, cluster_window,
                                              cluster_score, element_score, kmer, n_simulations, simulation_mode,
                                              simulation_window,
                                              cores, seed, log_level, qvalue, gzip))
    else:
        pass
