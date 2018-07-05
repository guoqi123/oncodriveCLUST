# Import modules
import logging

import numpy as np
import matplotlib.pyplot as plt

# Global variables
logs = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def prepare_clusters_data(strand, length, clusters_tree, false_margin=5):
    """
    Prepare clusters data to plot
    :param strand: element strand
    :param length: element length (total regions)
    :param clusters_tree: IntervalTree, clusters
    :param false_margin: int,
    :return:
            plot_max: list, position of maximum in total length
            plot_clusters: list of tuples, position of cluster in total length
    """
    plot_max = []
    plot_clusters = []

    for interval in clusters_tree:
        for key, value in interval.data.items():
            if strand == '-':
                plot_max.append(length - value['max'][0])
                plot_clusters.append(
                    (length - value['left_m'][0] - false_margin, length - value['right_m'][0] + false_margin))
            else:
                plot_max.append(value['max'][0])
                plot_clusters.append((value['left_m'][0] + false_margin, value['right_m'][0] - false_margin))

    return plot_max, plot_clusters


def plot(element, strand, smoothing, mutations, length, filter_clusters_tree):
    """

    :return:
    """
    plt.figure(figsize=(12, 4))
    ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)

    title = element

    if strand == '+':
        pass
    elif strand == '-':
        smoothing = smoothing[::-1]
        mutations = mutations[::-1]

    ax2 = ax1.twinx()
    ax2.plot(mutations, c='red', alpha=.5, linewidth=2.0)
    ax2.set_ylabel('Number of mutations', color='red', fontsize=20)

    ax1.plot(smoothing, c='royalblue', linewidth=3.0)
    ax1.set_ylabel('Smoothing score', color='royalblue', fontsize=20)

    ax1.tick_params('y', colors='royalblue', labelsize=16)
    ax2.tick_params('y', colors='red', labelsize=16)
    ax1.tick_params('x', labelsize=16)


    plot_max3, plot_clusters3 = prepare_clusters_data(strand=strand, length=length,
                                                  clusters_tree=filter_clusters_tree)
    for maxima in plot_max3:
        ax1.plot(maxima, 0.45, '.', ms=10, c='darkblue')
    for cluster in plot_clusters3:
        ax1.plot((cluster[0], cluster[1]), (0.45, 0.45), c='darkblue')

    ax1.set_axisbelow(True)
    ax1.grid(color='lightgrey', linestyle='--', linewidth=1)

    plt.title(title, fontsize=28)
    # plt.title('TERT (promoter)', fontsize=26)
    # plt.savefig('/home/carnedo/Claudia/Inkscape/20180703/TERT.svg', dpi=300, bbox_inches='tight')

    plt.show()


def run_plot(element, mutations, cds_d, strand, chromosome, smooth_window,
             smooth_tree, raw_clusters_tree, merge_clusters_tree, score_clusters_tree, element_score):
    """
    Prepare data and generate a plot for an element
    :param element:
    :param mutations:
    :param cds_d:
    :param strand:
    :param chromosome:
    :param smooth_window
    :param smooth_tree:
    :param raw_clusters_tree:
    :param merge_clusters_tree:
    :param score_clusters_tree:
    :param element_score:
    :return: None
    """

    smoothing_array = np.array([])
    for interval in sorted(smooth_tree):
        smoothing_array = np.append(smoothing_array, interval.data)

    length = len(smoothing_array)

    mutation_array = np.zeros(length)
    for mutation in mutations:
        index = (mutation.position - mutation.region[0]) + cds_d[mutation.region[0]].start
        mutation_array[index] += 1

    plot(element, strand, smoothing_array, mutation_array, length, score_clusters_tree)