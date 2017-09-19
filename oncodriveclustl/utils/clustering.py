# Import modules
import math as m
from collections import defaultdict


def find_locals(binary):
    """
    Find local maximum and minimum
    :param element_lists: dict, dictionary containing information for an element
    :return:
        idexes: list of tuples of 3 elements (min or max, cds region, smoothing score). Min == 0, max == 1
        maxs: list of positions of maximum in cds
    """
    indexes = []
    maxs = []
    length = len(binary)

    # for max or min in first position
    i = 0
    b = binary[i]
    c = binary[i + 1]
    if c <= b != 0:
        indexes.append([1, i, b])
        maxs.append(i)
    elif b <= c != 0:
        indexes.append([0, i, b])

    # for start +1, end -1 of binary array
    for i in range(1, length - 1):
        a = binary[i - 1]
        b = binary[i]
        c = binary[i + 1]
        # Max
        # When max score is equal in contiguous positions, assumes the first as max
        if a < b >= c:
            indexes.append([1, i, b])
            # Get a list of maximums sorted by position
            maxs.append(i)
        # Min
        elif a >= b <= c:
            # get a 'false' minimum where minimum score is equal in contiguous positions
            if a or c != 0:
                indexes.append([0, i, b])

    # for max or min in last position
    i = length - 1
    a = binary[i - 1]
    b = binary[i]
    if b > a:
        indexes.append([1, i, b])
        maxs.append(i)
    elif b < a:
        indexes.append([0, i, b])

    return indexes, maxs


def raw_clusters(indexes):
    """
    Define a cluster per maximum found in a region
    :param indexes: list of tuples of 3 elements (min or max, cds region, smoothing score). Min == 0, max == 1
    :return:
    """

    clusters = defaultdict(dict)

    # Iterate through all maxs in indexes
    j = 0
    generator_maxs = (i for i in indexes if i[0] == 1)

    for maximum in generator_maxs:
        i = indexes.index(maximum)
        clusters[j]['max'] = [maximum[1], maximum[2]]

        # if it's not the first nor the last cluster
        if i != 0 and i != len(indexes) - 1:
            clusters[j]['min_l'] = [indexes[i - 1][1], indexes[i - 1][2]]
            clusters[j]['min_r'] = [indexes[i + 1][1], indexes[i + 1][2]]

        # if it's the first cluster
        elif i == 0:
            clusters[j]['min_l'] = []
            clusters[j]['min_r'] = [indexes[i + 1][1], indexes[i + 1][2]]

        # if it's the last cluster
        elif i == len(indexes) - 1:
            clusters[j]['min_l'] = [indexes[i - 1][1], indexes[i - 1][2]]
            clusters[j]['min_r'] = []

        j += 1

    return clusters


def merge_clusters(maxs, clusters, window):
    """
    Given a number of clusters in a region, iterate through them to merge them if their maximums
    are closer than a given length.
    :param clusters: dict of dict
    :param maxs: list of maximum positions in a region
    :param window: clustering window
    :return:
    """
    maxs_set = set(maxs)

    iterate = 1
    while iterate != 0:  # Iterate until no clusters updates occur
        stop = 0
        for x in range(len(clusters.keys())):

            # When x is a key in clusters and min_r exists (clusters without min_r don't do merging):
            if x in clusters.keys():
                if clusters[x]['min_r']:
                    # Define the interval of search
                    search_r = set(range(clusters[x]['max'][0] + 1, clusters[x]['min_r'][0] + window + 1))

                    # When the intersection between the positions of the search and the maximums is not empty
                    if search_r.intersection(maxs_set) != set():
                        # Analyze only the closest max
                        intersect_max = maxs.index(sorted(list(search_r.intersection(maxs_set)))[0])
                        stop = 1

                        # When the testing max is greater or equal than the intersected max,
                        # expand the right border and delete intersected cluster from clusters
                        if clusters[x]['max'][1] >= clusters[intersect_max]['max'][1]:
                            clusters[x]['min_r'] = clusters[intersect_max]['min_r']
                            del clusters[intersect_max]
                            maxs_set.remove(maxs[intersect_max])

                        # When the testing max is smaller than the intersected max,
                        # expand the left border of intersected max and remove the testing cluster (i)
                        elif clusters[x]['max'][1] < clusters[intersect_max]['max'][1]:
                            clusters[intersect_max]['min_l'] = clusters[x]['min_l']
                            del clusters[x]
                            maxs_set.remove(maxs[x])

        if stop == 0:
            iterate = 0

    return clusters


def clusters_mut(clusters, genomic, mutations):
    """
    Calculates the number of mutations within a cluster
    :param clusters: dict of dict
    :return:
    """
    for cluster, values in clusters.items():

        if values['min_l'] and values['min_r']:
            left = int(genomic[values['min_l'][0]])
            right = int(genomic[values['min_r'][0]])

        elif not values['min_l']:
            left = int(genomic[values['max'][0]])
            right = int(genomic[values['min_r'][0]])

        else:   # elif not values['min_r']:
            left = int(genomic[values['min_l'][0]])
            right = int(genomic[values['max'][0]])

        cluster_muts = [i for i in mutations if left <= i <= right]
        clusters[cluster]['n_mutations'] = len(cluster_muts)

    return clusters


def score_clusters(clusters, mutations, mut_by_pos, method):
    """
    Score clusters
    :param clusters: dict of dict
    :param method: str, clustering scoring method
    :return:
    """

    root = m.sqrt(2)

    for cluster, values in clusters.items():
        score = 0

        # Define cluster borders
        if values['min_l'] and values['min_r']:
            a = values['min_l'][0]
            b = values['min_r'][0]
        elif not values['min_l']:
            a = values['max'][0]
            b = values['min_r'][0]
        else:
            a = values['min_l'][0]
            b = values['max'][0]

        # Calculate score per cluster iterating through positions
        denom = len(mutations)
        for position in range(a, b+1):
            muts = mut_by_pos[position]
            if method != 'nobias':
                muts /= denom
            distance = abs(values['max'][0] - position)
            score += (muts / m.pow(root, distance))

        # for position in range(a, b+1):
        #     mutations = element_lists['mut_by_pos'][position] / len(element_lists['mutations'])
        #     distance = abs(values['max'][0] - position)
        #     score.append(mutations / m.pow(root, distance))

        # Update score
        clusters[cluster]['score'] = score

    return clusters
