# Import modules
from collections import defaultdict

def get_genomic(regions):
    """Get a list of genomic positions
    :param regions: list of tuples with genomic positions of an element.
    :return genomic: list og genomic positions
    """

    genomic = []
    for interval in regions:
        positions = range(interval[0], interval[1] + 1)
        for position in positions:
            genomic.append(position)

    return genomic


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
