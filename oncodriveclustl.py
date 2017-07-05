# Import modules
import os
import gzip
import pickle
from functools import partial
from collections import defaultdict, namedtuple
from concurrent.futures import ProcessPoolExecutor as Pool

import click
import colorlog
import numpy as np
from tqdm import tqdm
from scipy.signal import argrelmax
from intervaltree import IntervalTree


# Configure the colorlog module
logger = colorlog.getLogger()


def set_logger(level):
    global logger
    d = {'info': colorlog.colorlog.logging.INFO,
         'warning': colorlog.colorlog.logging.WARNING,
         'error': colorlog.colorlog.logging.ERROR}
    logger.setLevel(d[level])
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter())
    logger.addHandler(handler)


def set_seed(seed):
    """Set the numpy seed generator
    :param seed: int, seed to use
    :return: None
    """
    np.random.seed(seed)


# Global variables
tukey_filter = None
Results = namedtuple('Results', ' '.join(['symbol', 'geneid', 'pvalues', 'pos_pvalues',
                                          'mutations', 'mutations_index', 'obs_max',
                                          'obs_peaks_values', 'obs_peaks_index', 'borders',
                                          'fragments', 'smoothed', 'num_mutations',
                                          'random_maxs']))
Mutation = namedtuple('Mutation', ' '.join(['chromosome', 'position', 'ref', 'alt',
                                            'sample', 'type', 'cancer_type']))


def tukey(window):
    """Tukey smoothing function
    :param window: int, length of the smoothing window
    :return: numpy array. The elements sum to 1
    """
    half_window = window // 2
    tukey = lambda x: np.maximum(1 - x ** 2, 0) ** 2
    filter_ = tukey(np.arange(-half_window, half_window + 1) / (half_window + 1))
    filter_ = filter_ / sum(filter_)
    return filter_


def build_regions(regions_file, description, elements=None):
    """Create a IntervalTree regions dictionary
    :param regions_file: path, path of the regions file
    :param description: str, description of the file
    :param elements: set, names of the elements to consider
    :return: dict, IntervalTree regions dictionary; dict, set of coordinates per region name
    """
    logger.info("Loading %s regions...", description)
    regions_tree = defaultdict(IntervalTree)
    regions_dict = defaultdict(set)
    with gzip.open(regions_file) as fd:
        for line in fd:
            line = line.decode().strip().split('\t')
            chromosome, start, end = line[:3]
            name = None if len(line) <= 3 else line[-1]
            if elements is not None and len(elements) > 0 and name not in elements:
                continue
            if chromosome.startswith('chr'):
                chromosome = chromosome[3:]
            regions_tree[chromosome][int(start): int(end) + 1] = name
            regions_dict[name].add((chromosome, int(start), int(end) + 1))
    return regions_tree, regions_dict


def map_mutations_to_regions(mutations_file, regions_tree, mutation_type, header=True):
    """Read a mutations file and return those mutations that fall inside the regions_tree
    :param mutations_file: path, file with mutations
    :param regions_tree: IntervalTree of genomic regions
    :param mutation_type: str, type of mutation ('subs', 'indel', 'all')
    :param header: boolean, if True the first line is skipped
    :return: dictionaries, mutations: keys are the elements name and values are mutations in the elements
    """
    logger.info('Mapping mutations ..')
    mutations = defaultdict(set)
    with open(mutations_file, 'r') as fd:
        if header is True:
            next(fd)
        for line in fd:
            mut = Mutation(*line.strip().split('\t')[:7])
            if mutation_type != 'all' and mut.type != mutation_type:
                continue
            for region in regions_tree[mut.chromosome][int(mut.position)]:
                mutations[region.data].add(mut)
    return mutations


def clusters(positions, fragments, window, geneid=None, simulate=False, signature=None):
    """Clustering function
    :param positions: list, position of the mutations
    :param fragments: list, fragments of a genomic region
    :param window: int, size of the smoothing window 
    :param geneid: str, name of the element
    :param simulate: boolean, if True simulates the mutations, else uses the observed mutations
    :param signature: dict or None, signatures to use in the simulations
    :return: binary (numpy.array), index (numpy.array), and borders (dictionary)
    """
    # TODO: this might be computed outside
    # These are the spaces where the mutations can occur
    whole_feature = np.array([])
    borders = {}
    for chromosome, start, end in sorted(fragments, key=lambda x: x[1]):
        begin = len(whole_feature)
        whole_feature = np.append(whole_feature, np.arange(start, end))
        borders[(begin, len(whole_feature))] = (start, end)

    # Turns it into 0s and 1s
    binary = np.zeros(len(whole_feature))
    # Deprecated: the following code is only good for unique positions
    # index = np.where(pd.Index(pd.unique(positions)).get_indexer(whole_feature) >= 0)[0]
    # This works with recurrent positions

    if simulate:
        positions = np.random.choice(
            a=whole_feature, size=len(positions), replace=True, p=signature
        )
        # if signature is not None:
        #     # Example of possible_positions if distance = 3: [0,0,0,1,1,1,2,2,2]
        #     distance = int(len(signature) / 3)
        #     diff = abs(distance - len(whole_feature))
        #     if diff > 0:
        #         print("{} distance difference {}".format(geneid, diff))
        #     tri_whole_feature = np.array([p for tripos in [[i, i, i] for i in range(distance)] for p in tripos])
        # else:
        #     tri_whole_feature = np.arange(len(binary))
        # index = np.random.choice(
        #     a=np.arange(len(binary)), size=len(positions), replace=True, p=signature
        # )

    index = np.searchsorted(whole_feature, positions)

    # Expands 'binary' at the begin and at the end to allow the window to cover the whole array
    binary = np.append(np.append(np.zeros(window // 2), binary), np.zeros(window // 2))

    # Apply the smoothing
    for i in index:
        try:
            binary[i: i + window] += tukey_filter
        except ValueError:
            print("\nError in gene {} with index {}".format(geneid, i))

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: -window // 2 + 1]
    return binary, index, borders


def cluster_element(elementid_mutations, regions, num_simulations, window_size, signature=None, absolute_max_peak=True):
    """Analyze a genomic element for the presence of clusters of mutations
    :param elementid_mutations: tuple, name of the element and its mutations
    :param regions: dictionary, coordinates of the regions
    :param num_simulations: int, number of simulations
    :param window_size: size of the window used in the smoothing
    :param signature: dictionary of trinucleotide signatures or None
    :param absolute_max_peak: boolean, if True only the highest cluster is considered
    :return: 
    """
    name, muts = elementid_mutations
    positions = [int(mut.position) for mut in muts]
    num_mutations = len(positions)
    num_samples = len(set([mut.sample for mut in muts]))

    # Observed clusters: these are the clusters found with the real mutations
    binary, index, borders = clusters(
        positions=positions,
        fragments=regions[name],
        window=window_size,
        geneid=name,
        simulate=False,
    )

    # Simulations
    if len(index) == 1:
        random_max_values = np.repeat(np.max(binary), num_simulations)
        pos_pvalues = np.ones(len(binary))
        pvalues = [np.nan]
        peaks_values = []
        peaks_index = []
    else:
        # Split mutations by tumor
        positions_by_tumors = defaultdict(list)
        for m in muts:
            positions_by_tumors[m.cancer_type].append(int(m.position))

        random_max_values = np.zeros(num_simulations)  # Max value of each randomization
        random_pos_values = np.zeros(len(binary))      # Max value of each position

        # TODO: load signatures
        for i in range(num_simulations):
            random_binary = None
            for tumor, tumor_specific_positions in positions_by_tumors.items():
                tumor_specific_binary, _, _ = clusters(
                    positions=tumor_specific_positions,
                    fragments=regions[name],
                    window=window_size,
                    geneid=name,
                    simulate=True,
                    signature=signature
                )
                if random_binary is None:
                    random_binary = tumor_specific_binary
                else:
                    random_binary += tumor_specific_binary
            # TODO: store the complete set of random_binary
            random_max_values[i] = np.max(random_binary)
            for x, value in enumerate(binary):
                random_pos_values[x] += int(random_binary[x] >= binary[x])

        if absolute_max_peak:
            peaks_index = np.where(binary == np.max(binary))
        else:
            # Get the cluster's peaks
            peaks_index = argrelmax(binary)[0]
            if np.max(binary) not in binary[peaks_index]:
                max_index = np.where(binary == np.max(binary))
                peaks_index = np.append(peaks_index, max_index)

            # If two peaks are closer than the length of the window, they will be joined.
            # The peaks to be joined first are the closest ones
            peaks_index.sort()
            distances = [peaks_index[i] - peaks_index[i - 1] for i in range(1, len(peaks_index))]
            while any([dist <= window_size for dist in distances]):
                min_distance = min(distances)
                min_index = distances.index(min_distance)
                if min_distance <= window_size:
                    peaks_index[min_index] = np.mean(peaks_index[min_index: min_index + 2])
                    peaks_index = np.delete(peaks_index, min_index + 1)
                # Update the distances
                distances = [peaks_index[i] - peaks_index[i - 1] for i in range(1, len(peaks_index))]

        # Peaks values
        peaks_values = binary[peaks_index]
        pvalues = [sum(random_max_values >= peak_value) / num_simulations
                   for peak_value in peaks_values]
        pos_pvalues = random_pos_values / num_simulations

    # Store the results
    result = dict(
        symbol=None,
        geneid=name,
        pvalues=pvalues,
        pos_pvalues=pos_pvalues,
        mutations=positions,
        mutations_index=index,
        obs_max=np.max(binary),
        obs_peaks_values=peaks_values,
        obs_peaks_index=peaks_index,
        borders=borders,
        fragments=regions[name],
        smoothed=binary,
        num_mutations=num_mutations,
        random_maxs=random_max_values
    )
    return name, result


def save(results, output_file):
    """Save the results as a text file and a pickle file
    :param results: dict, dictionary of namedtuple
    :param output_file: path, path of the files to create (a text file and a pickle file)
    :return: None
    """
    with open(output_file, 'w') as tfd, open(output_file + '.pickle', 'wb') as pfd:
        pickle.dump(results, pfd)
        tfd.write('\t'.join(['SYMBOL', 'ELEMENT_ID', 'PVALUE', 'MAX_SCORE']) + '\n')
        for gene, result in results.items():
            tfd.write('{}\t{}\t{}\t{}\n'.format(
                result.geneid, result.symbol, result.pvalues[0], result.obs_max
            ))


# Main function with CLI
@click.command()
@click.argument('mutations_file', type=click.Path(exists=True))  # , help='file with the mutations to simulate')
@click.argument('output_file')  # , help='path of the file to create with the simulated results')
@click.option('-r', '--regions-file', default=None, type=click.Path(exists=True), required=True,
              help='File with the genomic regions to analyze')
@click.option('-e', '--elements', multiple=True, default=None,
              help='Element(s) to analyze. Other elements in the "mutations_file" will be ignored. ' +
                   'By default all the elements in the --regions-list are analyzed')
@click.option('-w', '--window-size', type=click.INT, default=51,
              help='window used by the smoothing algorithm. Default is 51.')
@click.option('-m', '--mutation-type', default='subs', type=click.Choice(['subs', 'indel', 'all']),
              help='type of mutation to simulate: all, subs, or indel. Default is subs.')
@click.option('-S', '--signature', default=None, type=click.Path(exists=True),
              help='signatures to use in the simulations. Default is None.')
@click.option('--absolute-max-peak', is_flag=True,
              help='only consider the maximum peak of each element (genomic region)')
@click.option('--remove-unmappable', is_flag=True,
              help='do not simulate mutations in not mappable regions')
@click.option('--remove-blacklisted', is_flag=True,
              help='do not simulate mutations in UCSC blacklisted regions')
@click.option('--remove-low-complexity', is_flag=True,
              help='do not simulate mutations in low-complexity regions')
@click.option('--start-at-0', is_flag=True,
              help='if True the first position of a chromosome is 0, otherwise is 1')
@click.option('-n', '--num-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('-s', '--seed', type=click.INT, default=-1, help='seed of the random generator')
@click.option('--log-level', default='info', type=click.Choice(['info', 'warning', 'error']),
              help='verbosity of the logger. By default is info.')
def main(mutations_file, output_file, regions_file, window_size, mutation_type, signature,
         start_at_0, absolute_max_peak, remove_unmappable, remove_blacklisted,
         remove_low_complexity, num_simulations, elements, cores, seed, log_level):
    """Run the clustering algorithm"""
    set_logger(log_level)

    # Set the seed
    if seed >= 0:
        logger.info("Set seed to {}".format(seed))
        set_seed(seed)

    # Set the window size
    window_size += 1 - window_size % 2

    # Read the regions file as a dictionary of IntervalTree
    elements = set([]) if elements is None else set(elements)
    regions_tree, regions_dict = build_regions(regions_file, 'genomic', elements)
    regions_names = set(regions_dict.keys())  # set(sum([[j.data for j in i] for i in regions_tree.values()], []))

    if len(elements) > 0:
        for element in elements - regions_names:
            logger.warning('element %s not found', element)

    # Execute tukey and put the tukey_filter in the global namespace
    global tukey_filter
    tukey_filter = tukey(window_size)

    # Read the mutations file
    mutations = map_mutations_to_regions(
        mutations_file=mutations_file,
        regions_tree=regions_tree,
        mutation_type=mutation_type
    )

    results = {}
    with Pool(max_workers=cores) as executor, tqdm(total=len(mutations)) as pbar:
        fx = partial(
            cluster_element,
            regions=regions_dict,
            num_simulations=num_simulations,
            window_size=window_size,
            signature=signature,
            absolute_max_peak=absolute_max_peak
        )
        for geneid, result in executor.map(fx, mutations.items()):
            pbar.update(1)
            results[geneid] = Results(**result)
    save(results, output_file)


if __name__ == "__main__":
    main()

'''
def load_random_maxs(precalculated_path):
    """Load the randomizad values"""
    return np.fromfile(precalculated_path, dtype='float32').astype('float64')


def dump_random_maxs(random_values, precalculated_path):
    """Save the randomized values"""
    np.array(random_values, dtype='float32').tofile(precalculated_path, format='float32')
'''

# TODO:
#
