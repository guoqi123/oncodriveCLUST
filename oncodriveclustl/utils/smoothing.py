# Import modules
import numpy as np


def smooth(element_lists, window, tukey_filter):
    """Generate a smoothing curve for a list of element's mutations
    :return: dict, region_lists with new keys 'binary' and 'mut_by_pos'
        binary: length == genomic, smoothing score curve
        mut_by_pos: length == genomic, number of mutations per position
    """
    binary = np.zeros(len(element_lists['genomic']) + window - 1)  # binary cds with borders added
    mutations_a = np.zeros(len(element_lists['genomic']) + window - 1)  # mutations cds with borders added

    indexes = np.searchsorted(np.array(element_lists['genomic']), element_lists['mutations'])

    for index in indexes:
        binary[index: index + window] += tukey_filter
        mutations_a[index + window // 2] += 1

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: - window // 2 + 1]
    mutations_a = mutations_a[window // 2: - window // 2 + 1]

    element_lists['binary'] = binary
    element_lists['mut_by_pos'] = mutations_a

    return element_lists
