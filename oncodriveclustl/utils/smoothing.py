# Import modules
import numpy as np

from utils import sequence as seq


def smooth(regions_d, mutations, window, tukey_filter):
    """Generate a smoothing curve for a list of element's mutations
    :param regions_d: list of tuples with genomic positions of an element.
    :param mutations: list, list of mutations of an element
    :param window: int, smoothing window
    :param tukey_filter: numpy array. The elements sum to 1
    :return: dict, region_lists with new keys 'binary' and 'mut_by_pos'
        binary: length == genomic, smoothing score curve
    """

    # Get binary, cds with borders added
    binary = np.zeros(seq.get_length(regions_d) + window - 1)
    indexes = np.searchsorted(np.array(seq.get_genomic(regions_d)), mutations)
    for index in indexes:
        binary[index: index + window] += tukey_filter

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: - window // 2 + 1]

    return binary
