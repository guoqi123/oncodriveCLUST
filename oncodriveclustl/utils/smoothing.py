# Import modules
import numpy as np


def smooth(genomic, mutations, window, tukey_filter):
    """Generate a smoothing curve for a list of element's mutations
    :return: dict, region_lists with new keys 'binary' and 'mut_by_pos'
        binary: length == genomic, smoothing score curve
        mut_by_pos: length == genomic, number of mutations per position
    """
    binary = np.zeros(len(genomic) + window - 1)  # binary cds with borders added
    mut_by_pos = np.zeros(len(genomic) + window - 1)  # mutations cds with borders added

    indexes = np.searchsorted(np.array(genomic), mutations)

    for index in indexes:
        binary[index: index + window] += tukey_filter
        mut_by_pos[index + window // 2] += 1

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: - window // 2 + 1]
    mut_by_pos = mut_by_pos[window // 2: - window // 2 + 1]

    return binary, mut_by_pos
