# Import modules
from intervaltree import IntervalTree

import numpy as np


def smooth(regions, mutations, tukey_filter, simulation_window):
    """Generate a smoothing curve for a list of element's mutations
    :param regions: IntervalTree with genomic positions of an element
    :param mutations: list, list of mutations of an element
    :param tukey_filter: numpy array. Length equals smoothing window. The elements sum to 1
    :param simulation_window: int, simulation window
    :return:
        final_smooth_tree, IntervalTree. Interval are genomic regions, data np.array of smoothing score
        per position.
    """

    smooth_tree = IntervalTree()
    final_smooth_tree = IntervalTree()
    half_window = simulation_window // 2

    for interval in regions:
        # Add extra bases for smoothing tukey window
        begin = interval.begin - half_window
        end = interval.end + half_window
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
        begin = interval.begin + half_window
        end = interval.end - half_window
        slicer = (len(tukey_filter) // 2) + half_window
        final_smooth_tree.addi(begin, end, interval.data[slicer: - slicer])

    return final_smooth_tree
