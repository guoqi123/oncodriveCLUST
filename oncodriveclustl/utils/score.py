def element_score(clusters_tree, mode, method):
    """
    Given the clusters of an element, calculate a global score for it
    :param clusters_tree: IntervalTree, data are dict of dict
    :param mode: str, 'obs' for observed or 'sim' for simulated
    :param method: str, scoring method. Default 'mean'
    :return: float, element score
    """
    n_clusters = 0
    score = 0

    for interval in clusters_tree:
        clusters = interval.data
        for cluster, values in clusters.items():
            n_clusters += 1
            score += values['score']  # Add up cluster scores
            interval.data[cluster]['mode'] = mode

    return score
