# Oncodriveclustl functions 2
import numpy as np
from collections import defaultdict
import pickle

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


def simulate(array_positions, n_mutations, pos_prob):
    """Simulate mutations considering the signature
    :param array_positions: array genomic positions
    :param n_mutations: int, number of mutations to simulate
    :param pos_prob: array, raw probabilities
    :return: list
    """

    # Normalize probabilities
    norm_prob = normalize(pos_prob)

    # Calculate list of simulated mutations
    if len(array_positions) == len(norm_prob):
        sim_mutations = np.random.choice(array_positions, size=n_mutations, replace=True, p=norm_prob)
        return sim_mutations
    else:
        return []


def analysis(region_arrays, mutations, window_s, window_c):
    """
    Calculate smoothing, clustering and gene score
    :param region_arrays: dict, keys 'genomic', 'probs'.
    :param mutations: list, list of mutations, mutations == genomic positions
    :param window_s: int, smoothing window
    :param window_c: int, clustering window
    :return: int, gene score
    """
    region_smooth = odf.smoothing(region_arrays, mutations, window=window_s)

    clusters = odf.clustering(regions=region_smooth, mutations=mutations, window=window_c)
    n_clusters, gene_score = odf.score_gene(clusters=clusters, cutoff=2)

    return gene_score


def run_region(arguments, n_simulations=1000):
    """
    Given a gene, calculate the observed and simulated scores
    :param arguments: tuple (symbol, regions, chromosome, mutations)
        symbol: gene name
        regions: list of intervaltrees
    :param n_simulations: int, number of simulations
    :return:
        symbol: str, gene
        gene_results_d: dict, keys 'obs', 'sim'
    """
    gene_results_d = defaultdict()
    symbol, regions, chromosome, mutations = arguments

    # Calculate arrays
    """
    region_lists['symbol']
    region_lists['genomic']
    region_lists['probs']

    """
    region_arrays = odf.pre_smoothing(symbol, chromosome, regions, signatures)

    # Analysis of the observed mutations
    observed_results = analysis(region_arrays, mutations, window_s=50, window_c=50)
    gene_results_d['obs'] = observed_results

    # Analysis of the simulated mutations
    simulations = []
    for _ in range(n_simulations):
        simulated_mutations = simulate(
            array_positions=region_arrays['genomic'],
            n_mutations=len(mutations),
            pos_prob=region_arrays['probs'])
        simulated_results = analysis(region_arrays, simulated_mutations, window_s=50, window_c=50)
        simulations.append(simulated_results)

    gene_results_d['sim'] = simulations

    return symbol, gene_results_d
