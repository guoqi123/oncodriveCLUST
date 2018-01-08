# Import modules
import math as m
from collections import defaultdict

from oncodriveclustl.utils import sequence as seq


def find_locals(binary, regions):
    """
    Find local maximum and minimum
    :param binary: length == genomic, smoothing score curve
    :param regions: list of tuples with genomic positions of an element.
    :return:
        idexes: list of tuples of 4 elements (min or max, cds region, smoothing score, genomic position).
                min == 0, max == 1
        maxs: list of positions of maximum in cds
    """
    indexes = []
    maxs = []
    length = len(binary)
    genomic = seq.get_genomic(regions)

    # For first position
    b = binary[0]
    c = binary[1]
    if c < b:
        indexes.append((1, 0, b, genomic[0]))
        maxs.append(0)
    elif b < c:
        indexes.append((0, 0, b, genomic[0]))
    else:  # when b = c
        pass

    # For start +1, end -1 of binary array
    for i in range(1, length - 1):
        a = binary[i - 1]
        b = binary[i]
        c = binary[i + 1]
        # Max
        if a < b >= c:  # when max score is equal in contiguous positions, assumes the first as max
            indexes.append((1, i, b, genomic[i]))
            # Get a list of maximums sorted by position
            maxs.append(i)
        # Min
        elif a > b:
            if b < c:
                indexes.append((0, i, b, genomic[i]))
            if b == c:
                indexes.append((0, i, b, genomic[i]))
        elif a == b < c:
            indexes.append((0, i, b, genomic[i]))
        else:  # a = b = c
            pass

    # for max or min in last position
    i = length - 1
    a = binary[i - 1]
    b = binary[i]
    if b > a:
        indexes.append((1, i, b, genomic[i]))
        maxs.append(i)
    elif b < a:
        indexes.append((0, i, b, genomic[i]))
    else:
        pass

    return indexes, maxs


def raw_clusters(indexes):
    """
    Define a cluster per maximum found in a region
    :param indexes: list of tuples of 4 elements (min or max, cds region, smoothing score, genomic position).
    Min == 0, max == 1
    :return: clusters, dict of dict
    """

    clusters = defaultdict(dict)

    # Iterate through all maxs in indexes
    j = 0
    generator_maxs = (i for i in indexes if i[0] == 1)

    for maximum in generator_maxs:
        i = indexes.index(maximum)
        clusters[j]['max'] = (maximum[1], maximum[2], maximum[3])

        # if it's not the first nor the last cluster
        if i != 0 and i != len(indexes) - 1:
            clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
            clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])

        # if it's the first cluster
        elif i == 0:
            clusters[j]['left_m'] = (maximum[1], maximum[2], maximum[3])
            clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])

        # if it's the last cluster
        elif i == len(indexes) - 1:
            clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
            clusters[j]['right_m'] = (maximum[1], maximum[2], maximum[3])

        j += 1

    return clusters


def merge_clusters(clusters, maxs, window):
    """
    Given a number of clusters in a region, iterate through them to merge them if their maximums
    are closer than a given length.
    :param clusters: dict of dict
    :param maxs: list of maximum positions in a region
    :param window: clustering window
    :return:
        clusters: dict of dict
    """
    maxs_set = set(maxs)

    iterate = 1
    while iterate != 0:  # Iterate until no clusters updates occur
        stop = 0
        for x in range(len(clusters.keys())):
            # When x is a key in clusters
            if x in clusters.keys():
                maximum = clusters[x]['max']
                left_margin = clusters[x]['left_m']
                right_margin = clusters[x]['right_m']

                if right_margin != maximum:
                    # Define the interval of search
                    search_r = set(range(maximum[0] + 1, right_margin[0] + window + 1))

                    # When the intersection between the positions of the search and the maximums is not empty
                    if search_r.intersection(maxs_set) != set():

                        # Analyze only the closest cluster
                        intersect_cluster = maxs.index(sorted(list(search_r.intersection(maxs_set)))[0])
                        stop = 1

                        # When the testing max is greater or equal than the intersected max,
                        # expand the right border and delete intersected cluster from clusters
                        if maximum[1] >= clusters[intersect_cluster]['max'][1]:
                            clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                            del clusters[intersect_cluster]
                            maxs_set.remove(maxs[intersect_cluster])

                        # When the testing max is smaller than the intersected max,
                        # expand the left border of intersected max and remove the testing cluster (i)
                        elif maximum[1] < clusters[intersect_cluster]['max'][1]:
                            clusters[intersect_cluster]['left_m'] = left_margin
                            del clusters[x]
                            maxs_set.remove(maxs[x])
                        else:  # they are equal
                            pass
                else:  # they are equal
                    pass

        if stop == 0:
            iterate = 0

    return clusters


def clusters_mut(clusters, regions, mutations):
    """
    Calculates the number of mutations within a cluster
    :param clusters: dict of dict
    :param regions: list of tuples with genomic positions of an element.
    :param mutations: list, list of mutations
    :return:
        clusters: dict of dict
    """

    genomic = seq.get_genomic(regions)

    for cluster, values in clusters.items():
        # Define cluster borders
        left = int(genomic[values['left_m'][0]])
        right = int(genomic[values['right_m'][0]])
        # Search mutations
        cluster_muts = [i for i in mutations if left <= i <= right]
        clusters[cluster]['n_mutations'] = len(cluster_muts)

    return clusters


def score_clusters(clusters, mutations, regions, method):
    """
    Score clusters
    :param clusters: dict of dict
    :param mutations: list, list of mutations
    :param regions: list of tuples with genomic positions of an element.
    :param method: str, clustering scoring method
    :return:
    """

    root = m.sqrt(2)
    mut_by_pos = seq.get_mutations(mutations)
    genomic = seq.get_genomic(regions)

    for cluster, values in clusters.items():
        score = 0

        # Calculate score per cluster iterating through positions
        denom = len(mutations)
        for position in range(values['left_m'][0], values['right_m'][0]+1):
            coordinate = genomic[position]
            muts = mut_by_pos.get(coordinate, 0)
            if method != 'nobias':
                muts /= denom
            distance = abs(values['max'][0] - position)
            score += (muts / m.pow(root, distance))

        # Update score
        clusters[cluster]['score'] = score

    return clusters
