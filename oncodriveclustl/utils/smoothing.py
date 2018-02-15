# Import modules
from intervaltree import IntervalTree

import numpy as np


def smooth(regions, mutations, tukey_filter):
    """Generate a smoothing curve for a list of element's mutations
    :param regions: IntervalTree with genomic positions of an element
    :param mutations: list, list of mutations of an element
    :param tukey_filter: numpy array. The elements sum to 1
    :return: dict, region_lists with new keys 'binary' and 'mut_by_pos'
        smooth: length == genomic, smoothing score curve
    """

    smooth_tree = IntervalTree()
    final_smooth_tree = IntervalTree()

    for interval in regions:
        # Add extra bases for smoothing tukey window
        begin = interval.begin
        end = interval.end
        # Remove 1 to tukey filter to allow symmetric // 2
        smooth_tree.addi(begin, end, np.zeros((end - begin) + len(tukey_filter) - 1))

    # Find mutations in regions
    for mutation in mutations:
        for interval in smooth_tree[mutation]:
            # Get index of mutation in region
            index = mutation - interval.begin
            # Smooth mutations
            interval.data[index: index + len(tukey_filter)] += tukey_filter

    for interval in smooth_tree:
        final_smooth_tree.addi(interval.begin, interval.end,
                               interval.data[(len(tukey_filter) // 2): - (len(tukey_filter) // 2)])

    return final_smooth_tree
