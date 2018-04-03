# Import modules
from intervaltree import IntervalTree
from collections import namedtuple

import numpy as np


def smooth(element, regions, cds_d, mutations, tukey_filter, simulation_window):
    """Generate a smoothing curve for a list of element's mutations
    :param regions: IntervalTree with genomic positions of an element
    :param cds_d: dict, keys are start genomic regions, values are cds positions
    :param mutations: list, list of mutations formatted as namedtuple
    :param tukey_filter: numpy array. Length equals smoothing window. The elements sum to 1
    :param simulation_window: int, simulation window
    :return:
        final_smooth_tree, IntervalTree. Interval are genomic regions or cds, data np.array of smoothing score
        by position.
        mutations_in: list of mutations in regions
    """
    first_smooth_tree = IntervalTree()
    final_smooth_tree = IntervalTree()
    mutations_in = []

    """
    print('---------------------------------------------------')
    # Calculate length
    element_length = 0
    for interval in regions:
        element_length += (interval[1] - interval[0])
    print('element: ', element_length)
    print('tukey: ', len(tukey_filter))
    print('simulation: ', simulation_window)
    """

    # Generate smoothing arrays for regions
    for interval in regions:
       # print('region ', interval)
        # Add extra bases for smoothing of simulated mutations that fall outside regions and tukey_filter
        first_smooth_tree.addi(interval.begin, interval.end,
                         np.zeros((interval.end - interval.begin + len(tukey_filter) + simulation_window - 2)))

    if not cds_d:
        # Smooth
        for mutation in mutations:
            for interval in first_smooth_tree[mutation.region[0]]:
                # Get index of mutation in region
                new_begin = interval.begin - (simulation_window + len(tukey_filter) - 2)//2   # always integer
                index = mutation.position - new_begin
                tukey_begin = index - (len(tukey_filter) - 1)//2
                # Smooth mutations
                interval.data[tukey_begin: tukey_begin + len(tukey_filter)] += tukey_filter
            # Get mutations inside regions
            if regions[mutation.position]:
                mutations_in.append(mutation)

        # Remove extra bp
        for interval in first_smooth_tree:
            begin = interval.begin
            end = interval.end
            slicer = (simulation_window + len(tukey_filter) - 2)//2
            final_smooth_tree.addi(begin, end, interval.data[slicer: - slicer])

    else:
        # Smooth simulated mutations outside regions
        for mutation in mutations:
            # TODO TABIX change
            if mutation.muttype == 1:
                if not first_smooth_tree[mutation.position]:
                    for interval in first_smooth_tree[mutation.region[0]]:
                        new_begin = interval.begin - (simulation_window + len(tukey_filter) - 2) // 2  # always integer
                        index = mutation.position - new_begin
                        tukey_begin = index - (len(tukey_filter) - 1) // 2
                        # Smooth mutations
                        interval.data[tukey_begin: tukey_begin + len(tukey_filter)] += tukey_filter

        # Remove extra bp
        for interval in first_smooth_tree:
            begin = interval.begin
            end = interval.end
            slicer = (simulation_window + len(tukey_filter) - 2) // 2
            final_smooth_tree.addi(begin, end, interval.data[slicer: - slicer])

        # Merge sorted regions (one interval == cds) and add tukey//2 to both ends
        cds_tree = IntervalTree()
        cds_array = np.zeros((len(tukey_filter)-1) // 2)
        for interval in sorted(final_smooth_tree):
            cds_array = np.append(cds_array, interval.data)
        cds_array = np.append(cds_array, np.zeros((len(tukey_filter)-1) // 2))
        cds_tree.addi(final_smooth_tree.begin(), final_smooth_tree.end(), cds_array)
        final_smooth_tree = IntervalTree()

        # Smooth mutations inside regions
        for mutation in mutations:
            # TODO TABIX change
            if mutation.muttype == 1:
                if first_smooth_tree[mutation.position]:
                    for interval in cds_tree[mutation.position]:
                        # Get index of mutation in cds
                        index = (mutation.position - mutation.region[0]) + cds_d[mutation.region[0]].start
                        # Smooth mutations
                        interval.data[index: (index + len(tukey_filter))] += tukey_filter
                    mutations_in.append(mutation)

        # Remove extra bp
        for interval in cds_tree:
            begin = interval.begin
            end = interval.end
            slicer = (len(tukey_filter) -1) // 2
            final_smooth_tree.addi(begin, end, interval.data[slicer: -slicer])

    return final_smooth_tree, mutations_in
