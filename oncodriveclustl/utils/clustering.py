# Import modules
import math as m
from collections import defaultdict
from collections import namedtuple

import numpy as np
from scipy.signal import argrelextrema
from intervaltree import IntervalTree


def find_locals(smooth_tree, cds_d):
    """
    Find local maximum and minimum of a smoothing curve
    :param smooth_tree: IntervalTree. Interval are genomic regions or cds, data np.array of smoothing score by position.
    :param cds_d: dict, keys are start genomic regions, values are cds positions
    :return:
    index_tree, IntervalTree, data are lists of tuples containing information from local maximum and minimum
    """
    index_tree = IntervalTree()
    reverse_cds_d = IntervalTree()

    # Reverse dictionary into IntervalTree
    if cds_d:
        for genomic, cds in cds_d.items():
            reverse_cds_d.addi(cds[0], cds[1] + 1, genomic)  # + 1 because end not included

    for interval in smooth_tree:
        indexes_info = []
        smooth = interval.data
        # Greater or equal
        max_eq = argrelextrema(smooth, np.greater_equal, order=1, mode='clip')[0].tolist()
        # Smaller or equal
        min_eq = argrelextrema(smooth, np.less_equal, order=1, mode='clip')[0].tolist()
        # Find maximums and minimums
        indexes = list(set(max_eq).symmetric_difference(set(min_eq)))

        # Add information to index
        for index in sorted(indexes):
            score = smooth[index]
            local = 0 if index in min_eq else 1  # 0 if minimum, 1 if maximum
            if not cds_d:
                genomic = interval.begin + index
            else:
                for info in reverse_cds_d[index]:
                    genomic = info.data + index - info[0]
            indexes_info.append((local, index, genomic, score))

        # Add to new interval tree
        index_tree.addi(interval[0], interval[1], indexes_info)

    return index_tree


def raw_clusters(index_tree):
    """
    Define a cluster per maximum found
    :param index_tree: IntervalTree, data are lists of tuples of 4 elements (min or max, cds region, genomic position,
                       smoothing score).
    :return: clusters_tree, IntervalTree, data are dict of dict
    """

    clusters_tree = IntervalTree()

    for interval in index_tree:
        clusters = defaultdict(dict)
        j = 0
        indexes = interval.data

        # Iterate through all maximum and generate a cluster per maximum
        generator_maxs = (i for i in indexes if i[0] == 1)
        for maximum in generator_maxs:
            i = indexes.index(maximum)
            # Add maximum
            clusters[j]['max'] = (maximum[1], maximum[2], maximum[3])
            # Add margins
            # if maximum not in first nor last position
            if i != 0 and i != len(indexes) - 1:
                # if no contiguous left max
                if indexes[i - 1][0] != 1:
                    clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
                else:
                    clusters[j]['left_m'] = (maximum[1], maximum[2], maximum[3])
                # if no contiguous right max
                if indexes[i + 1][0] != 1:
                    clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])
                else:
                    clusters[j]['right_m'] = (maximum[1], maximum[2], maximum[3])
            # if first position
            elif i == 0:
                clusters[j]['left_m'] = (maximum[1], maximum[2], maximum[3])
                clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])
            # if last position
            else:
                clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
                clusters[j]['right_m'] = (maximum[1], maximum[2], maximum[3])
            j += 1
        clusters_tree.addi(interval[0], interval[1], clusters)

    return clusters_tree


def merge_clusters(clusters_tree, window):
    """
    Given clusters, iterate through them to merge them if their maximums are closer than a given length.
    :param clusters_tree, IntervalTree, data are dict of dict
    :param window: clustering window
    :return:
        clusters: dict of dict
    """
    merged_clusters_tree = IntervalTree()
    missed_clusters = defaultdict(dict)

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = interval.data.copy()
        n_clusters = len(clusters.keys())
        # Get the maximums to iterate
        maxs = []
        for c, v in clusters.items():
            maxs.append(v['max'][0])
        maxs_set = set(maxs)

        # Start iterations
        iterate = 1
        while iterate != 0:  # Iterate until no clusters updates occur
            stop = 0
            for x in range(n_clusters):
                # When x is a key in clusters
                if x in clusters.keys():
                    maximum = clusters[x]['max']
                    left_margin = clusters[x]['left_m']
                    right_margin = clusters[x]['right_m']

                    # Define the interval of search
                    search_r = set(range(right_margin[0], right_margin[0] + window + 1))
                    # If there is a maximum in the search region
                    if search_r.intersection(maxs_set):
                        # Analyze only the closest cluster
                        intersect_clusters = sorted(list(search_r.intersection(maxs_set)))
                        if maxs.index(intersect_clusters[0]) == x:
                            if len(intersect_clusters) > 1:
                                intersect_cluster = maxs.index(intersect_clusters[1])
                            else:
                                intersect_cluster = maxs.index(intersect_clusters[0])
                        else:
                            intersect_cluster = maxs.index(intersect_clusters[0])
                        stop = 1
                        # When the testing max is greater than the intersected max,
                        # expand the right border and delete intersected cluster from clusters
                        if maximum[2] > clusters[intersect_cluster]['max'][2]:
                            clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                            del clusters[intersect_cluster]
                            maxs_set.remove(maxs[intersect_cluster])
                        # When the testing max is smaller than the intersected max,
                        # expand the left border of intersected max and remove the testing cluster (i)
                        elif maximum[2] < clusters[intersect_cluster]['max'][2]:
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
                            elif maximum[0] < (clusters[intersect_cluster]['max'][0]):
                                # Do not merge
                                missed_clusters[x] = clusters[x]
                                del clusters[x]
                                maxs_set.remove(maxs[x])
                            # If same cluster
                            else:
                                # Don't iterate again
                                maxs_set.remove(maxs[x])
            if stop == 0:
                iterate = 0

        for k, v in missed_clusters.items():
            clusters[k] = v
        merged_clusters_tree.addi(interval[0], interval[1], clusters)
        missed_clusters = defaultdict(dict)

    return merged_clusters_tree


def clusters_mut(clusters_tree, mutations_in, cluster_mutations_cutoff):
    """
    Get the number of mutations within a cluster, remove those below cutoff mutations
    :param clusters_tree, IntervalTree, data are dict of dict
    :param mutations_in: list of mutations in regions
    :param cluster_mutations_cutoff: int, n cluster mutations cutoff
    :return:
        filter_clusters_tree: IntervalTree, data are dict of dict
    """
    filter_clusters_tree = IntervalTree()

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = interval.data.copy()
        for cluster, values in interval.data.items():
            left = values['left_m'][1]
            right = values['right_m'][1]
            # Search mutations
            cluster_muts = [i for i in mutations_in if left <= i.position <= right]
            cluster_samples = set()
            for mut in cluster_muts:
                sample = mut.sample
                cluster_samples.add(sample)
            if len(cluster_muts) >= cluster_mutations_cutoff:
                clusters[cluster]['mutations'] = cluster_muts
                clusters[cluster]['samples'] = cluster_samples
                clusters[cluster]['fra_uniq_samples'] = len(cluster_samples)/len(cluster_muts)
            else:
                del clusters[cluster]
        filter_clusters_tree.addi(interval[0], interval[1], clusters)

    return filter_clusters_tree


def fmutations_score(clusters_tree, regions, mutations_in, protein):
    """
    Score clusters with fraction of mutations formula
    :param clusters_tree: IntervalTree, data are dict of dict
    :param regions: IntervalTree with genomic positions of an element
    :param mutations_in: int, len of list of mutations inside regions
    :param protein: bool, True analyzes clustering in protein sequence
    :return:
            score_clusters_tree: IntervalTree, data are dict of dict
    """
    score_clusters_tree = IntervalTree()
    root = m.sqrt(2)

    # Update regions if protein
    if protein:
        regions = IntervalTree()
        regions.addi(clusters_tree.begin(), clusters_tree.end())

    for interval in clusters_tree:
        clusters = interval.data.copy()
        for cluster, values in clusters.items():
            score = 0
            mutated_positions_d = defaultdict(int)

            # Get number of mutations on each mutated position
            for mutation in values['mutations']:
                mutated_positions_d[mutation.position] += 1

            # Map mutated position and smoothing maximum to region
            for position, count in mutated_positions_d.items():
                map_mut_pos = set()
                map_smo_max = set()

                if regions[position]:
                    for i in regions[position]:
                        map_mut_pos = i
                    for i in regions[values['max'][1]]:
                        map_smo_max = i
                    # Calculate distance of position to smoothing maximum
                    if map_mut_pos[0] == map_smo_max[0]:
                        distance_to_max = abs(position - values['max'][1])
                    elif map_mut_pos[0] < map_smo_max[0]:
                        distance_to_max = (map_mut_pos[1] - position) + (values['max'][1] - map_smo_max[0])
                    else:
                        distance_to_max = (map_smo_max[1] - values['max'][1]) + (position - map_mut_pos[0])
                    # Calculate fraction of mutations
                    mutations = (count / mutations_in) * 100

                    # Calculate cluster score
                    denom = m.pow(root, distance_to_max)
                    score += (mutations / denom)

            # Update
            clusters[cluster]['score'] = score

        score_clusters_tree.addi(interval[0], interval[1], clusters)

    return score_clusters_tree
