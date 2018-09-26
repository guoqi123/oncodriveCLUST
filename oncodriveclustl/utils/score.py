"""
Contains function to score genomic elements
"""


def element_score(clusters_tree, mode):
    """
    Calculate the score of a genomic element by adding up the scores of its clusters.

    Args:
        clusters_tree (IntervalTree): IntervalTree, data are dict of dict
        mode (str): 'obs' for observed or 'sim' for simulated clusters

    Returns:
        score (float): element score
    """
    n_clusters = 0
    score = 0

    for interval in clusters_tree:
        clusters = interval.data
        for cluster, values in clusters.items():
            n_clusters += 1
            score += values['score']
            interval.data[cluster]['mode'] = mode

    return score
