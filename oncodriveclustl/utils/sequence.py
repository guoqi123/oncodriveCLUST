# Import modules
from collections import defaultdict


def get_length(regions):
    """Calculate genomic length of an element
    :param regions: list of tuples with genomic positions of an element.
    :return int, genomic length of an element
    """

    length = 0
    for interval in regions:
        length += len(range(interval[0], interval[1] + 1))
    return length


def get_mutations(mutations):
    """Generate a dictionary of mutations
    :param mutations: list of mutations as genomic positions
    :return dict, keys are positions >= 1 mutation, values are number of mutations
    """

    mut_by_pos = defaultdict(int)
    for coordinate in mutations:
        mut_by_pos[coordinate] += 1
    return mut_by_pos
