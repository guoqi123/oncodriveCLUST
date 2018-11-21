"""
Contains function to apply KDE smoothing of mutations
"""
from intervaltree import IntervalTree
from collections import namedtuple

import numpy as np

Mutation = namedtuple('Mutation', 'position, region, alt, sample, cancertype')


def smooth_nucleotide(regions, concat_regions_d, mutations, tukey_filter, simulation_window):
    """Generate a smoothing curve for a list of element's mutations in the nucleotide sequence

    Args:
        regions (IntervalTree): IntervalTree with genomic positions of an element
        concat_regions_d (dict): keys are start genomic regions, values are positions (index) relative to the start
        mutations (list): list of mutations formatted as namedtuple
        tukey_filter (numpy.ndarray): kde array, length equals smoothing window.
        simulation_window (int): simulation window

    Returns:
        final_smooth_tree (IntervalTree): interval are genomic regions or indexes (concatenate mode),
            data np.array of smoothing score by position
        mutations_in (list): list of mutations in regions
    """
    first_smooth_tree = IntervalTree()
    final_smooth_tree = IntervalTree()
    mutations_in = []

    # Generate smoothing arrays for regions
    for interval in regions:
        # Add extra bases for smoothing of simulated mutations that fall outside regions and tukey_filter
        first_smooth_tree.addi(interval.begin, interval.end,
                               np.zeros((interval.end - interval.begin + len(tukey_filter) + simulation_window - 2)))

    if not concat_regions_d:
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

        # Merge sorted regions (one interval == concatenated sequence) and add tukey//2 to both ends
        concat_tree = IntervalTree()
        concat_array = np.zeros((len(tukey_filter) - 1) // 2)
        for interval in sorted(final_smooth_tree):
            concat_array = np.append(concat_array, interval.data)
        concat_array = np.append(concat_array, np.zeros((len(tukey_filter) - 1) // 2))
        concat_tree.addi(final_smooth_tree.begin(), final_smooth_tree.end(), concat_array)
        final_smooth_tree = IntervalTree()

        # Smooth mutations inside regions
        for mutation in mutations:
            if first_smooth_tree[mutation.position]:
                for interval in concat_tree[mutation.position]:
                    # Get index of mutation in concatenated sequence
                    index = (mutation.position - mutation.region[0]) + concat_regions_d[mutation.region[0]].start
                    # Smooth mutations
                    interval.data[index: (index + len(tukey_filter))] += tukey_filter
                mutations_in.append(mutation)

        # Remove extra bp
        for interval in concat_tree:
            begin = interval.begin
            end = interval.end
            slicer = (len(tukey_filter) - 1) // 2
            final_smooth_tree.addi(begin, end, interval.data[slicer: -slicer])

    return final_smooth_tree, mutations_in
