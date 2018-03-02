# Import modules
import os
import logging
from intervaltree import IntervalTree

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import click
import daiquiri

# Global variables
logs = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def array_mutations(regions, mutations, tukey_filter, simulation_window):
    """Generate an array for a list of element's mutations
    :param regions: IntervalTree with genomic positions of an element
    :param mutations: list, list of mutations of an element
    :param tukey_filter: numpy array. Length equals smoothing window. The elements sum to 1
    :param simulation_window: int, simulation window
    :return:
        final_mutation_tree, IntervalTree. Interval are genomic regions, data np.array of number of mutations by
        position.
    """


def plot_element(regions, mutations, tukey_filter, simulation_window):
    """
    Plot element clustering
    :param regions: IntervalTree with genomic positions of an element
    :param mutations: list, list of mutations of an element
    :param tukey_filter: numpy array. Length equals smoothing window. The elements sum to 1
    :param simulation_window: int, simulation window
    :return: None
    """

    # Get smoothing

    title = 'OR6K2 COAD'
    mutations = COAD
    elements = set(['OR6K2'])
    simw = 50
    regions_d, chromosomes_d, strands_d, mutations_d = pars.parse(regions, elements, mutations)
    print(mutations_d)
    smooth_tree, mutation_tree = smooth(regions_d['OR6K2_ENSG00000196171'], mutations_d['OR6K2_ENSG00000196171'],
                                        tukey_filter, simw)

    smoothing_array = np.array([])
    mutation_array = np.array([])

    for interval in sorted(smooth_tree):
        smoothing_array = np.append(smoothing_array, interval.data)
    for interval in sorted(mutation_tree):
        mutation_array = np.append(mutation_array, interval.data)

    plt.figure(figsize=(24, 6))
    ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)

    if strands_d['OR6K2_ENSG00000196171'] == '+':
        pass
    elif strands_d['OR6K2_ENSG00000196171'] == '-':
        smoothing_array = smoothing_array[::-1]
        mutation_array = mutation_array[::-1]

    ax1.plot(smoothing_array, c='darkblue')
    ax1.set_ylabel('smoothing score', color='darkblue', fontsize=22)

    ax2 = ax1.twinx()
    ax2.plot(mutation_array, c='red', alpha=.5)
    ax2.set_ylabel('n mutations', color='red', fontsize=22)

    ax1.tick_params('y', colors='darkblue', labelsize=18)
    ax2.tick_params('y', colors='red', labelsize=18)

    plt.title(title, fontsize=26)

    plt.show()
