# Import modules
from collections import defaultdict


def element_score(clusters, cutoff, mode, method):
    """
    Given the clusters of a region, calculate a global score for it
    :param clusters: dictionary of dictionaries
    :param cutoff: int, n cluster mutations cutoff
    :param mode: str, 'obs' for observed or 'sim' for simulated
    :param method: str, scoring method. Default 'mean'
    :return: int, int; number of clusters above cutoff and element score
    """
    n_clusters = []
    score = []
    cutoff_clusters = defaultdict(dict)

    for cluster, value in clusters.items():
        if value['n_mutations'] >= cutoff:
            n_clusters.append(1)
            score.append(value['score'])
            # get clusters information
            value['mode'] = mode
            cutoff_clusters[cluster] = value

    if method == 'mean':
        if n_clusters:
            element_score = sum(score) / sum(n_clusters)
        else:
            element_score = 0
    elif method == 'sum':
            element_score = sum(score)
    else:
        element_score = 0

    # if n_clusters:
    #     element_score = sum(score) / sum(n_clusters)
    # else:
    #     element_score = 0

    return cutoff_clusters, element_score
