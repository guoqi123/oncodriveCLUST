# Import modules
import math as m
from collections import defaultdict

import numpy as np
from scipy.signal import argrelextrema
from intervaltree import IntervalTree

from oncodriveclustl.utils import smoothing as smo


def find_locals(regions, mutations, tukey_filter, cds, simulation_window):
    """
    Find local maximum and minimum of a smoothing curve
    :param regions: IntervalTree with genomic positions of an element
    :param mutations: list, list of mutations of an element
    :param tukey_filter: numpy array. The elements sum to 1
    :param cds: bool, True calculates clustering on cds
    :param simulation_window: int, simulation window
    :return:
    # TODO return
    """
    index_tree = IntervalTree()
    length_tree = IntervalTree()
    smooth_tree = smo.smooth(regions, mutations, tukey_filter, simulation_window)

    if cds:
        # Re-write smooth_tree (one interval == cds)
        cds_tree = IntervalTree()
        cds_array = np.array([])
        start = 0
        for interval in sorted(smooth_tree):
            length = interval.end - interval.begin
            length_tree.addi(start, start + length, interval[0])
            start = start + length + 1
            cds_array = np.append(cds_array, interval.data)
        cds_tree.addi(smooth_tree.begin(), smooth_tree.end(), cds_array)
        smooth_tree = cds_tree

    for interval in smooth_tree:
        indexes_info = []
        smooth = interval.data
        # Greater or equal
        max_eq = argrelextrema(smooth, np.greater_equal, order=1, mode='wrap')[0].tolist()
        # Smaller or equal
        min_eq = argrelextrema(smooth, np.less_equal, order=1, mode='wrap')[0].tolist()
        # Find maximums and minimums
        indexes = list(set(max_eq).symmetric_difference(set(min_eq)))
        # Add information to index
        for index in sorted(indexes):
            score = smooth[index]
            local = 0 if index in min_eq else 1  # 0 if minimum, 1 if maximum
            indexes_info.append((local, index, score))
        # Add to new interval tree
        index_tree.addi(interval[0], interval[1], indexes_info)

    return index_tree, length_tree, smooth_tree


def raw_clusters(index_tree):
    """
    Define a cluster per maximum found in a region
    # TODO change
    :param index_tree: list of tuples of 4 elements (min or max, cds region, smoothing score, genomic position).
    Min == 0, max == 1
    :return: clusters, dict of dict
    """

    clusters_tree = IntervalTree()

    for interval in index_tree:
        clusters = defaultdict(dict)
        j = 0
        indexes = interval.data     # TODO add key to indexes

        # Iterate through all maximum and generate a cluster per maximum
        generator_maxs = (i for i in indexes if i[0] == 1)
        for maximum in generator_maxs:
            i = indexes.index(maximum)
            # Add maximum
            clusters[j]['max'] = (maximum[1], maximum[2])
            # Add margins
            # if maximum not in first nor last position
            if i != 0 and i != len(indexes) - 1:
                clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2])
                clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2])
            # if first position
            elif i == 0:
                clusters[j]['left_m'] = (maximum[1], maximum[2])
                clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2])
            # if last position
            else:
                clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2])
                clusters[j]['right_m'] = (maximum[1], maximum[2])
            j += 1
        clusters_tree.addi(interval[0], interval[1], clusters)

    return clusters_tree


def merge_clusters(clusters_tree, window):
    """
    # TODO
    Given a number of clusters in a region, iterate through them to merge them if their maximums
    are closer than a given length.
    :param clusters: dict of dict
    :param window: clustering window
    :param cds: bool, True calculates clustering on cds
    :return:
        clusters: dict of dict
    """
    merged_clusters_tree = IntervalTree()

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = interval.data
        n_clusters = len(clusters.keys())

        # for k, v in clusters.items():
        #    print(k, v)

        # Get the maximums to iterate
        maxs = []
        for c, v in clusters.items():
            maxs.append(v['max'][0])
        maxs_set = set(maxs)
        # print(maxs)

        # Start iterations
        iterate = 1
        while iterate != 0:  # Iterate until no clusters updates occur
            stop = 0
            for x in range(n_clusters):

                # When x is a key in clusters
                if x in clusters.keys():
                    #print('iterate with cluster', x)
                    maximum = clusters[x]['max']
                    left_margin = clusters[x]['left_m']
                    right_margin = clusters[x]['right_m']
                    if right_margin != maximum:
                        # Define the interval of search
                        search_r = set(range(right_margin[0], right_margin[0] + window + 1))
                        # If there is a maximum in the search region
                        if search_r.intersection(maxs_set):
                            # Analyze only the closest cluster
                            intersect_cluster = maxs.index(sorted(list(search_r.intersection(maxs_set)))[0])
                            stop = 1
                            # When the testing max is greater than the intersected max,
                            # expand the right border and delete intersected cluster from clusters
                            if maximum[1] > clusters[intersect_cluster]['max'][1]:
                                clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                                del clusters[intersect_cluster]
                                maxs_set.remove(maxs[intersect_cluster])
                            # When the testing max is smaller than the intersected max,
                            # expand the left border of intersected max and remove the testing cluster (i)
                            elif maximum[1] < clusters[intersect_cluster]['max'][1]:
                                clusters[intersect_cluster]['left_m'] = left_margin
                                del clusters[x]
                                maxs_set.remove(maxs[x])
                            # When the testing max is equal than the intersected max
                            else:
                                # If contiguous maximum's positions, choose the first maximum
                                if maximum[0] == (clusters[intersect_cluster]['max'][0] - 1):
                                    # Expand the right border and delete intersected cluster from clusters
                                    clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                                    del clusters[intersect_cluster]
                                    maxs_set.remove(maxs[intersect_cluster])
                                # If not contiguous positions
                                else:
                                    # Do not merge
                                    # Delete testing cluster (i)
                                    del clusters[x]
                                    maxs_set.remove(maxs[x])
                    else:  # do not need to merge
                        pass

            if stop == 0:
                iterate = 0
        # print('--------------')
        # for k, v in clusters.items():
        #     print(k, v)
        # print('--------------')
        #

        merged_clusters_tree.addi(interval[0], interval[1], clusters)

    return merged_clusters_tree


def get_genomic_position(index_in_cds, length_tree):
    """
    # TODO
    Get genomic position of an index in cds
    :param index_in_cds:
    :param length_tree:
    :return: genomic_in_cds
    """

    # Check if index in cds_tree, empty set when index == interval.end
    correction = 0
    if not length_tree[index_in_cds]:
        correction = 1

    # Find the genomic region where index lays
    region = list(length_tree[index_in_cds - correction])[0][2]

    # Calculate distance from region to cluster index
    distance = index_in_cds - list(length_tree[index_in_cds - correction])[0][0]

    # Define cluster borders in genomic positions
    genomic_in_cds = region + distance

    return genomic_in_cds


def clusters_mut(clusters_tree, length_tree, mutations, cds):
    """
    # TODO
    Get the number of mutations within a cluster
    :param clusters: dict of dict
    :param regions: list of tuples with genomic positions of an element.
    :param mutations: list, list of mutations
    :param cds: bool, True calculates clustering on cds
    :return:
        clusters: dict of dict
    """

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = dict(interval.data)
        for cluster, values in clusters.items():
            # Define cluster borders in genomic positions
            if cds:
                left = get_genomic_position(values['left_m'][0], length_tree)
                right = get_genomic_position(values['right_m'][0], length_tree)
                maximum = get_genomic_position(values['max'][0], length_tree)
            else:
                left = interval[0] + values['left_m'][0]
                right = interval[0] + values['right_m'][0]
                maximum = values['max'][0] + interval[0]
            # Search mutations
            cluster_muts = [i for i in mutations if left <= i <= right]
            cluster_muts = dict((m, cluster_muts.count(m)) for m in cluster_muts)
            interval.data[cluster]['mutations'] = cluster_muts
            # Update clusters information
            interval.data[cluster]['left_m'] = (values['left_m'][0], left)
            interval.data[cluster]['max'] = (values['max'][0], maximum)
            interval.data[cluster]['right_m'] = (values['right_m'][0], right)

    return clusters_tree


def score_clusters(clusters_tree, length_tree, regions_tree, cutoff, method, cds, total_mutations):
    """
    # TODO
    Score clusters
    :param clusters_tree:
    :param index_tree:
    :param cutoff: int, n cluster mutations cutoff
    :param method: str, clustering scoring method
    :param cds: bool, True calculates clustering on cds
    :param total_mutations: int
    :return:
    """
    # TODO cds and !cds in one
    root = m.sqrt(2)

    if not cds:
        # Iterate through all regions
        for interval in clusters_tree:
            clusters = dict(interval.data)
            for cluster, values in clusters.items():
                # For clusters above cutoff mutations
                if sum(values['mutations'].values()) >= cutoff:
                    score = 0
                    # Calculate score iterating through mutated positions
                    for mutation, count in values['mutations'].items():
                        # Calculate distance of position to smoothing maximum
                        distance_to_max = abs(mutation - values['max'][1])
                        # Calculate fraction of mutations if score is fmutations
                        mutations = (count / total_mutations) * 100
                        # Calculate cluster score
                        score += (mutations / m.pow(root, distance_to_max))
                    # Update
                    interval.data[cluster]['score'] = score
                else:
                    del interval.data[cluster]

    else:
        # TODO change
        reverse_length_tree = defaultdict()
        for interval in length_tree:
            key = interval.data
            reverse_length_tree[key] = (interval[0], interval[1])

        # Analyze cds region
        for interval in clusters_tree:
            clusters = dict(interval.data)
            for cluster, values in clusters.items():
                # For clusters above cutoff mutations
                if sum(values['mutations'].values()) >= cutoff:
                    score = 0
                    # Calculate score iterating through mutated positions
                    for mutation, count in values['mutations'].items():
                        # Get mutations genomic region index of mutation in region
                        for region in regions_tree[mutation]:  # 1 region
                            # Get index of mutation in cds (index of region start + index in region)
                            index = reverse_length_tree[region.begin][0] + (mutation - region.begin)
                            # Calculate distance of position to smoothing maximum
                            distance_to_max = abs(index - values['max'][0])
                            # Calculate fraction of mutations if score is fmutations
                            mutations = (count / total_mutations) * 100
                            # Calculate cluster score
                            score += (mutations / m.pow(root, distance_to_max))
                    # Update
                    interval.data[cluster]['score'] = score
                else:
                    del interval.data[cluster]

    return clusters_tree
