# Oncodriveclustl functions 2
import os.path
import numpy as np
from collections import defaultdict
import pickle
import logging
import daiquiri
from concurrent.futures import ProcessPoolExecutor as Pool

import oncodcl_funct as odf



# Global variables
import_signatures = pickle.load(open("/home/carnedo/projects/oncodriveclustl/signatures/SKCM.pickle", "rb"))
signatures = import_signatures['probabilities']


def normalize(probs):
    """
    Given an array of probabilities, normalize them to 1
    :param probs: array
    :return: array of normalized probabilities
    """
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]


def analysis(region_lists):
    """
    Calculate smoothing, clustering and gene score
    :param region_lists: dict, keys 'genomic', 'probs', 'mutations', 'windows'
    :return: int, gene score
    """
    region_smooth = odf.smoothing(region_lists, window=region_lists['windows'][0])
    clusters = odf.clustering(regions=region_smooth, window=region_lists['windows'][1])
    n_clusters, gene_score = odf.score_gene(clusters=clusters, cutoff=2)

    return gene_score


def simulate(region_lists, pos_prob, n_mut, n_sim):
    """Simulate mutations considering the signature
    :param region_lists: array genomic positions
    :param pos_prob: array, raw probabilities
    :param n_mut: int, number of mutations to simulate
    :param n_sim: int, number of simulations
    :return: list
    """

    # Normalize probabilities
    norm_prob = normalize(pos_prob)

    # Generate simulated mutations
    if len(region_lists['genomic']) == len(norm_prob):
        for i in range(n_sim):
            region_lists['mutations'] = np.random.choice(region_lists['genomic'], size=n_mut, replace=True, p=norm_prob)
            yield region_lists
    else:
        for i in range(n_sim):
            logger.debug('Pre-smoothing...')
            region_lists['mutations'] = []
            yield region_lists


def get_p(observed, simulations):
    """
    Calculate p value
    :return: float, p-value
    """
    return len([x for x in simulations if x >= observed]) / len(simulations)


def run_region(arguments):
    """
    Given a gene, calculate the observed and simulated scores
    :param arguments: tuple (symbol, regions, chromosome, mutations,
    n_simulations, smooth_window, cluster_window, cores)
        symbol: gene name
        regions: list of intervaltrees
        chromosome: str, chromosome of the genomic element analysed
        mutations: list, list of mutations (genomic positions)
        n_simulations: int, number of simulations
        smooth_window: int, smoothing window
        cluster_window: int, clustering window
        cores: int, number of cpus
    :return:
        symbol: str, gene
        observed: float, observed gene score
        p_value: float, p-value
    """
    global logger
    logger = daiquiri.getLogger()

    symbol, regions, chromosome, mutations, n_simulations, smooth_window, cluster_window, cores = arguments

    # Read signatures
    if os.path.isfile('signature.pickle'):
        signatures = pickle.load(open("signature.pickle", "rb"))
        signatures = signatures['probabilities']
        logger.debug('Signatures read')

        # Calculate lists of genomic positions and probabilities
        logger.debug('Pre-smoothing...')
        region_lists = odf.pre_smoothing(symbol, chromosome, regions, signatures)
        region_lists['mutations'] = mutations
        region_lists['windows'] = [smooth_window, cluster_window]
        logger.debug('Pre-smoothing done')

        # Analysis of the observed mutations
        observed = analysis(region_lists)
        logger.debug('Observed mutations analyzed')


        # Simulate mutations (generator)
        sim_scores = []
        simulations = simulate(region_lists,region_lists['probs'],len(mutations),n_simulations)
        # Multiprocess
        logger.debug('Start simulations...')
        with Pool(max_workers=cores) as executor:
            for score in executor.map(analysis, simulations):
                sim_scores.append(score)

        p_value = get_p(observed, sim_scores)

        return symbol, observed, p_value

    else:
        logger.critical('File \'signatures.pickle\' not found')