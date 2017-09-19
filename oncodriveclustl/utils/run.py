# Import modules
import os.path
import re
from concurrent.futures import ProcessPoolExecutor as Pool

import daiquiri
import pickle
from collections import defaultdict
from bgreference import hg19
from tqdm import tqdm
import numpy as np

from oncodriveclustl.utils import smoothing as smo
from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import analyticalpval as ap


# Logger
logger = daiquiri.getLogger()


class Experiment:
    """Class to analyze elements of a cancer dataset"""

    def __init__(self, regions_d, chromosomes_d, mutations_d,
                 element_mutations_cutoff, cluster_mutations_cutoff, smooth_window, cluster_window,
                 cluster_score, element_score, n_simulations, cores, seed):
        """Initialize the Experiment class
        :param regions_d: dict, dictionary containing genomic regions from all analyzed elements
        :param chromosomes_d: dict, dictionary containing chromosomes from all analyzed elements
        :param mutations_d: dict, dictionary containing mutations lists from all analyzed elements
        :param element_mutations_cutoff: int, cutoff of element mutations
        :param cluster_mutations_cutoff: int, cutoff of cluster mutations
        :param smooth_window: int, smoothing window
        :param cluster_window: int, clustering window
        :param cluster_score: cluster score method
        :param element_score: element score method
        :param n_simulations: int, number of simulations
        :param cores: int, number of CPUs to use
        :param seed: int, seed
        :return: None
        """

        self.regions_d = regions_d
        self.chromosomes_d = chromosomes_d
        self.mutations_d = mutations_d
        self.element_mutations_cutoff = element_mutations_cutoff
        self.cluster_mutations_cutoff = cluster_mutations_cutoff
        self.smooth_window = smooth_window + (1 - smooth_window % 2)
        # Calculate tukey filter
        self.tukey_filter = self.tukey(self.smooth_window)
        self.cluster_window = cluster_window
        self.cluster_score = cluster_score
        self.element_score = element_score
        self.n_simulations = n_simulations
        self.cores = cores
        self.seed = seed

        # Read CGC
        # TODO Remove this hardcoded file
        with open(os.path.join(os.path.dirname(__file__), '../data/CGCMay17_cancer_types_TCGA.tsv'), 'r') as fd:
            self.cgc_genes = set([line.split('\t')[0] for line in fd])

    @staticmethod
    def tukey(window):
        """Tukey smoothing function generates tukey_filter for smoothing
        :return: numpy array. The elements sum to 1
        """
        half_window = window // 2
        tukey = lambda x: np.maximum(1 - x ** 2, 0) ** 2
        filter_ = tukey(np.arange(-half_window, half_window + 1) / (half_window + 1))
        filter_ = filter_ / sum(filter_)
        return filter_

    @staticmethod
    def normalize(probs):
        """
        Given an array of probabilities, normalize them to 1
        :param probs: array
        :return: array of normalized probabilities
        """
        prob_factor = 1 / sum(probs)
        return [prob_factor * p for p in probs]

    def pre_smoothing(self, element):
        """
        Generate genomic, trinucleotides and probabilities lists of an element
        :param element: element to calculate pre-smoothing
        :return: element_lists: dictionary containing information for an element. Keys:
                    genomic: list containing all genomic positions in the element analyzed
                    probabilities: list of length == genomic, contains normalized mutation probabilities per position
                    mutations: list containing genomic positions of mutations within an element
        """

        nucleot = {'A', 'C', 'G', 'T'}
        genomic = []
        probabilities = []

        n1 = re.compile('[N]')
        n2 = re.compile('[N]{2}')

        # Read signatures.pickle

        # TODO remove this hardcoded path
        pickle_path = os.path.join(os.path.dirname(__file__),'../cache/signature.pickle')
        if os.path.isfile(pickle_path):
            signatures = pickle.load(open(pickle_path, "rb"))
            signatures = signatures['probabilities']
            logger.debug('Signatures read')
        else:
            raise RuntimeError('File \'signatures.pickle\' not found')

        # Iterate through tuples of coordinates, get genomic regions and sequence
        for pos in self.regions_d[element]:
            positions = range(pos[0], pos[1] + 1)
            sequence = hg19(self.chromosomes_d[element], pos[0] - 1, len(positions) + 2)  # start -1, end +1

            # Search trinucleotide probabilities
            for n in range(1, len(sequence) - 1):  # start to end
                ref_tri = sequence[n - 1:n + 2]
                two_n = n2.search(ref_tri)
                one_n = n1.search(ref_tri)
                if two_n is None:
                    # No N
                    if one_n is None:
                        prob = 0
                        for nuc in nucleot.difference({ref_tri[1]}):
                            search = (ref_tri, ref_tri[0] + nuc + ref_tri[2])
                            prob += signatures[search]
                        probabilities.append(prob)
                    # One N
                    else:
                        new_tris = set()
                        prob = 0
                        for nuc in nucleot:
                            new_tris.add(re.sub('N', nuc, ref_tri))
                        for reference in new_tris:
                            for change in new_tris.difference(set([reference])):
                                search = (reference, change)
                                prob += signatures[search]
                        probabilities.append(prob)
                # Two or three N
                else:
                    probabilities.append(0)
            # Get genomic regions
            for position in positions:
                genomic.append(position)

        return genomic, self.normalize(probabilities), self.mutations_d[element]

    def analysis(self, genomic, mutations, analysis_mode = 'sim'):
        """
        Calculate smoothing, clustering and element score for either observed or simulated mutations
        :param element_lists: dict, keys 'genomic', 'probs', 'mutations'
        :param analysis_mode: observed or simulated; default simulated
        :return:
            dict of dicts: cutoff clusters
            int, gene score
        """

        binary, mut_by_pos = smo.smooth(genomic, mutations, window=self.smooth_window, tukey_filter=self.tukey_filter)
        logger.debug('Smoothing calculated')

        indexes, maxs = clu.find_locals(binary.tolist())
        r_clusters = clu.raw_clusters(indexes=indexes)
        m_clusters = clu.merge_clusters(maxs=maxs, clusters=r_clusters, window=self.cluster_window)
        mut_clusters = clu.clusters_mut(m_clusters, genomic, mutations)
        scored_clusters = clu.score_clusters(mut_clusters, mutations, mut_by_pos, self.cluster_score)
        logger.debug('Clusters calculated')

        cutoff_clusters, gene_score = score.element_score(clusters=scored_clusters, cutoff=self.cluster_mutations_cutoff, mode=analysis_mode, method=self.element_score)
        logger.debug('Element score calculated')

        return cutoff_clusters, gene_score

    def simulate_and_analysis(self, item):
        element, genomic, probs, n_sim = item

        try:
            sim_scores_chunk = []
            sim_cluster_chunk = []

            mutations = np.random.choice(genomic, size=(n_sim, len(self.mutations_d[element])), p=probs)
            for i in range(n_sim):
                cutoff_clusters, gene_score = self.analysis(genomic, mutations[i])
                sim_scores_chunk.append(gene_score)
                sim_cluster_chunk.append(len(cutoff_clusters))

            return element, sim_scores_chunk, sim_cluster_chunk
        except Exception as e:
            return element, None, e

    @staticmethod
    def empirical_pvalue(observed, simulations):
        """
        Calculate empirical p-value
        :return: float, p-value
        """
        return (len([x for x in simulations if x >= observed]) + 1) / len(simulations)

    def post_process(self, item):
        try:
            element, genomic, mutations, sim_clusters_list, sim_scores_list = item

            # Statistics: number of clusters per simulation
            mean_sim_clusters = np.mean(sim_clusters_list)
            std_sim_clusters = np.std(sim_clusters_list)
            median_sim_clusters = np.median(sim_clusters_list)

            # Statistics: element score per simulation
            mean_sim_score = np.mean(sim_scores_list)
            std_sim_score = np.std(sim_scores_list)
            median_sim_score = np.median(sim_scores_list)

            # Calculate observed mutations results
            obs_clusters, obs_score = self.analysis(genomic, mutations, analysis_mode='obs')
            logger.debug('Observed mutations analyzed')

            # Significance
            empirical_pvalue = self.empirical_pvalue(obs_score, sim_scores_list)
            sim_scores_list_1000 = np.random.choice(sim_scores_list, size=1000, replace=False)
            pvalue_generator = ap.AnalyticalPvalue()
            analytical_pvalue, bandwidth = pvalue_generator.calculate(obs_score, sim_scores_list_1000)
            logger.debug('P-values calculated')

            # Get GCG boolean
            cgc = element in self.cgc_genes

            return element, \
                (len(genomic), len(self.mutations_d[element]), len(obs_clusters.keys()),
                    mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, mean_sim_score,
                    median_sim_score, std_sim_score, empirical_pvalue, analytical_pvalue, cgc), \
                (obs_clusters, cgc)
        except Exception as e:
            return element, None, e

    def run(self):
        """
        Analyze elements
        :return:
        """
        elements_results = defaultdict()
        clusters_results = defaultdict()

        # Filter elements >= cutoff mutations
        all_elements = [elem for elem in self.regions_d.keys() if len(self.mutations_d[elem]) >= self.element_mutations_cutoff]
        logger.info('Calculating results {} element{}...'.format(len(all_elements), 's' if len(all_elements) > 1 else ''))

        # Performance tunning
        pf_num_elements = 100
        pf_num_simulations = 100000

        # Run analysis on each element
        all_elements = list(chunkizator(all_elements, pf_num_elements))
        for i, elements in enumerate(all_elements, start=1):
            simulations = []
            genomic = {}
            mutations = {}
            logger.info("Iteration {} of {}".format(i, len(all_elements)))

            for element in elements:

                # Get element pre-smoothing, common for observed and simulated mutations
                logger.debug('Calculating pre-smoothing...')
                genomic[element], probs, mutations[element] = self.pre_smoothing(element)

                # Calculate simulated mutations results
                logger.debug('Start simulations...')

                element_size = ((len(self.mutations_d[element]) * self.n_simulations) // pf_num_simulations) + 1
                chunk_size = (self.n_simulations // element_size) + 1
                simulations += [(element, genomic[element], probs, n_sim) for n_sim in partitions_list(self.n_simulations, chunk_size)]

            with Pool(max_workers=self.cores) as executor:

                # Process simulations
                sim_scores_list = defaultdict(list)
                sim_clusters_list = defaultdict(list)
                for element, sim_scores_chunk, sim_cluster_chunk in tqdm(
                        executor.map(self.simulate_and_analysis, simulations), total=len(simulations), desc="simulations".rjust(30)
                ):
                    if sim_scores_chunk is not None:
                        sim_scores_list[element] += sim_scores_chunk
                        sim_clusters_list[element] += sim_cluster_chunk
                    else:
                        logger.error("At element {}. {}".format(element, str(sim_cluster_chunk)))

                # Post processing
                post_item = [(e, genomic[e], mutations[e], sim_clusters_list[e], sim_scores_list[e]) for e in elements]
                for e, er, cr in tqdm(executor.map(self.post_process, post_item), total=len(post_item), desc="post processing".rjust(30)):
                    if er is not None:
                        elements_results[e] = er
                        clusters_results[e] = cr
                    else:
                        logger.error("At element {}. {}".format(e, str(cr)))

        return elements_results, clusters_results


def partitions_list(total_size, chunk_size):
    """
    Create a list of values less or equal to chunk_size that sum total_size

    :param total_size: Total size
    :param chunk_size: Chunk size
    :return: list of integers
    """
    partitions = [chunk_size for _ in range(total_size // chunk_size)]

    res = total_size % chunk_size
    if res != 0:
        partitions += [res]

    return partitions


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable:
        size (int): elements in the chunk

    Returns:
        list. Chunk

    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk