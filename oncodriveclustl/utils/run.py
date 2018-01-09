# Import modules
import os.path
from concurrent.futures import ProcessPoolExecutor as Pool

import daiquiri
import pickle
from collections import defaultdict
from tqdm import tqdm
import numpy as np

from oncodriveclustl.utils import smoothing as smo
from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import sequence as seq
from oncodriveclustl.utils import analyticalpval as ap

# Logger
logger = daiquiri.getLogger()


class Experiment:
    """Class to analyze elements of a cancer dataset"""

    def __init__(self, regions_d, chromosomes_d, mutations_d, genome, path_pickle,
                 element_mutations_cutoff, cluster_mutations_cutoff, smooth_window, cluster_window,
                 cluster_score, element_score, kmer, n_simulations, simulation_mode, simulation_window, cores, seed):
        """Initialize the Experiment class
        :param regions_d: dict, dictionary containing genomic regions from all analyzed elements
        :param chromosomes_d: dict, dictionary containing chromosomes from all analyzed elements
        :param mutations_d: dict, dictionary containing mutations lists from all analyzed elements
        :param genome: genome to use
        :param element_mutations_cutoff: int, cutoff of element mutations
        :param cluster_mutations_cutoff: int, cutoff of cluster mutations
        :param smooth_window: int, smoothing window
        :param cluster_window: int, clustering window
        :param cluster_score: cluster score method
        :param element_score: element score method
        :param kmer: int, number of nucleotides of the signature
        :param n_simulations: int, number of simulations
        :param cores: int, number of CPUs to use
        :param seed: int, seed
        :return: None
        """

        self.regions_d = regions_d
        self.chromosomes_d = chromosomes_d
        self.mutations_d = mutations_d
        self.genome_build = self.load_genome(genome)
        self.path_pickle = path_pickle
        self.element_mutations_cutoff = element_mutations_cutoff
        self.cluster_mutations_cutoff = cluster_mutations_cutoff
        self.smooth_window = smooth_window + (1 - smooth_window % 2)
        # Calculate tukey filter
        self.tukey_filter = self.tukey(self.smooth_window)
        self.cluster_window = cluster_window
        self.cluster_score = cluster_score
        self.element_score = element_score
        self.kmer = kmer
        self.n_simulations = n_simulations
        self.simulation_mode = simulation_mode
        self.simulation_window = simulation_window

        self.cores = cores
        self.seed = seed

        # Read CGC
        if self.genome_build.__name__ == 'hg19':
            # TODO Remove this hardcoded file
            with open(os.path.join(os.path.dirname(__file__), '../data/CGCMay17_cancer_types_TCGA.tsv'), 'r') as fd:
                self.cgc_genes = set([line.split('\t')[0] for line in fd])
        else:
            self.cgc_genes = set()

    @staticmethod
    def tukey(window):
        """Tukey smoothing function generates tukey_filter for smoothing
        :param window: smoothing window
        :return: numpy array. The elements sum to 1
        """
        half_window = window // 2
        tukey = lambda x: np.maximum(1 - x ** 2, 0) ** 2
        filter_ = tukey(np.arange(-half_window, half_window + 1) / (half_window + 1))
        filter_ = filter_ / sum(filter_)
        return filter_

    @staticmethod
    def load_genome(genome):
        """
        Load the appropriate genome build
        :param genome: str, genome
        :return: instance of the bgreference
        """
        if genome == 'hg19':
            from bgreference import hg19 as genome_build
        elif genome == 'mm10':
            from bgreference import mm10 as genome_build
        elif genome == 'c3h':
            from bgreference import c3h as genome_build
        elif genome == 'car':
            from bgreference import car as genome_build
        else:
            from bgreference import cast as genome_build

        return genome_build

    @staticmethod
    def normalize(element, probs):
        """
        Given an array of probabilities, normalize them to 1
        :param element: genomic element of analysis
        :param probs: array of probabilities
        :return: array of normalized probabilities
        """

        try:
            prob_factor = 1 / sum(probs)
            return [prob_factor * p for p in probs]
        except ZeroDivisionError as e:
            logger.error('{}. No mutational probabilities dervied from signatures in element {}'.format(e, element))
            return False


    def mut_probabilities(self, element, kmer):
        """
        Generate mutational probabilities per position of an element based on the sequence context observed signature
        :param element: element to calculate pre-smoothing
        :param kmer: int, number of nucleotides of the signature
        :return:
            probabilities: list of length == genomic, contains not normalized mutational probabilities per position
        """
        regions = 0
        nu = 0
        sequ = 0
        delta = 1 if kmer == 3 else 2
        nucleot = {'A', 'C', 'G', 'T'}
        probabilities = []

        # Read signatures pickle
        if os.path.isfile(self.path_pickle):
            signatures = pickle.load(open(self.path_pickle, "rb"))
            signatures = signatures['probabilities']
            logger.debug('Signatures read')
        else:
            logger.critical('Signatures pickle {} not found'.format(self.path_pickle))
            quit()

        # Iterate through tuples of coordinates, get genomic regions' sequence
        for pos in self.regions_d[element]:
            regions += 1
            start = pos[0]-delta
            end = pos[1] - pos[0] + 1 + delta*2
            sequence = self.genome_build(self.chromosomes_d[element], start, end)  # genomic start -d, genomic end +d
            sequ += len(sequence)-2

            # Search kmer probabilities
            for n in range(delta, len(sequence)-delta):  # start to end
                nu += 1
                ref_kmer = sequence[n - delta: n + delta + 1]
                N = ref_kmer.count('N')
                if N == 0:
                    prob = 0
                    for alt in nucleot.difference({ref_kmer[kmer//2]}):  # mutation probability to any other kmer
                        alt_kmer = ref_kmer[: kmer//2] + alt + ref_kmer[kmer//2 + 1:]
                        prob += signatures[(ref_kmer, alt_kmer)]
                # elif N == 1:
                #     prob = 0
                #     n_index = ref_kmer.index('N')
                #     for nuc in nucleot:  # All possible reference kmers
                #         new_ref = ref_kmer[:n_index] + nuc + ref_kmer[n_index + 1:]
                #         for alt in nucleot.difference({new_ref[kmer//2]}):
                #             alt_kmer = new_ref[:kmer//2] + alt + new_ref[kmer//2 + 1:]
                #             prob += signatures[(new_ref, alt_kmer)]
                else:
                    prob = 0
                probabilities.append(prob)

        # Check probabilities
        # if probabilities:
        #     # length == genomic?
        #     # sum() == 0?
        # else:
        #     logger.

        return probabilities

    def analysis(self, element, mutations, analysis_mode='sim'):
        """
        Calculate smoothing, clustering and element score for either observed or simulated mutations
        :param element: element of analysis.
        :param mutations: list, list of mutations of an element
        :param analysis_mode: str, observed or simulated; default simulated
        :return:
            dict of dicts: cutoff clusters
            int, gene score
        """

        binary = smo.smooth(self.regions_d[element], mutations, window=self.smooth_window,
                            tukey_filter=self.tukey_filter)
        logger.debug('Smoothing calculated')
        indexes, maxs = clu.find_locals(binary.tolist(), self.regions_d[element])
        r_clusters = clu.raw_clusters(indexes=indexes)
        m_clusters = clu.merge_clusters(maxs=maxs, clusters=r_clusters, window=self.cluster_window)
        mut_clusters = clu.clusters_mut(m_clusters, self.regions_d[element], mutations)
        scored_clusters = clu.score_clusters(mut_clusters, mutations, self.regions_d[element], self.cluster_score)
        logger.debug('Clusters calculated')
        cutoff_clusters, element_score = score.element_score(
            clusters=scored_clusters,cutoff=self.cluster_mutations_cutoff,mode=analysis_mode,method=self.element_score
        )
        logger.debug('Element score calculated')

        return cutoff_clusters, element_score

    def simulate_and_analysis(self, item):
        """
        Simulate mutations and analyze simulations
        :param item: tuple, element of analysis data
        :return:
            element: str, element of analysis
            sim_scores_chunk: list, simulated element's results
            sim_cluster_chunk: list, simulated cluster's results
        """
        element, probs, n_sim = item
        mode = self.simulation_mode
        sim_scores_chunk = []
        sim_cluster_chunk = []
        # TODO: remove genomic
        genomic = seq.get_genomic(self.regions_d[element])

        if mode == 'hotspot':
            index_n_mutations = defaultdict(int)
            segment_mutations = {}
            regions = []

            # Get mutation indexes in genomic
            for mutation in self.mutations_d[element]:
                index_n_mutations[genomic.index(mutation)] += 1
            # Get window
            half_window = self.simulation_window // 2
             # Get range of genomic coordinates per mutation
            for index, n_mut in index_n_mutations.items():
                a = index - half_window
                b = index + half_window

                if a < 0:
                    a = genomic[0]
                else:
                    a = genomic[a]
                if b >= len(genomic):
                    b = genomic[-1]
                else:
                    b = genomic[b]

                segment_mutations[tuple([a, b])] = n_mut
                regions.append(tuple([a, b]))

            # Simulate mutations inside segments of regions
            initializer = 0
            for segment in regions:
                # Get segment start and end indexes in genomic
                a = genomic.index(segment[0])
                b = genomic.index(segment[1])
                # Get number of mutations in the segment
                n_mut_segment = segment_mutations[segment]
                # Simulate mutations
                if n_mut_segment != 0:   # remove? now there are no segments, so it's not needed

                    # TODO: check normalization, by default replace=True
                    mut_segment = np.random.choice(genomic[a:b], size=(n_sim, n_mut_segment),
                                                   p=self.normalize(element, probs[a:b]))
                    if initializer != 0:
                        mutations = np.hstack([mutations, mut_segment])
                    else:
                        mutations = mut_segment
                    initializer = 1

        else:  # mode == 'element'
            mutations = np.random.choice(genomic, size=(n_sim, len(self.mutations_d[element])), p=self.normalize(probs))

        for i in range(n_sim):
            cutoff_clusters, gene_score = self.analysis(element, mutations[i])
            sim_scores_chunk.append(gene_score)
            sim_cluster_chunk.append(len(cutoff_clusters))

        return element, sim_scores_chunk, sim_cluster_chunk

    @staticmethod
    def empirical_pvalue(observed, simulations):
        """
        Calculate empirical p-value
        :param observed: int, observed score
        :param simulations: list, simulated score
        :return: float, p-value
        """
        # TODO: pseudocount
        return (len([x for x in simulations if x >= observed]) + 1) / len(simulations)

    def post_process(self, item):
        """
        Post processing of results
        :param item: tuple, elements results
        :return:
        """

        element, obs_clusters, obs_score, sim_clusters_list, sim_scores_list = item

        if obs_clusters and type(obs_clusters) != float:
            # Clusters
            total_clusters = len(obs_clusters.keys())

            # Statistics: number of clusters per simulation
            mean_sim_clusters = np.mean(sim_clusters_list)
            std_sim_clusters = np.std(sim_clusters_list)
            median_sim_clusters = np.median(sim_clusters_list)

            # Statistics: element score per simulation
            mean_sim_score = np.mean(sim_scores_list)
            std_sim_score = np.std(sim_scores_list)
            median_sim_score = np.median(sim_scores_list)

            # Significance
            empirical_pvalue = self.empirical_pvalue(obs_score, sim_scores_list)
            sim_scores_list_1000 = np.random.choice(sim_scores_list, size=1000, replace=False)
            pvalue_generator = ap.AnalyticalPvalue()
            analytical_pvalue, bandwidth = pvalue_generator.calculate(obs_score, sim_scores_list_1000)
            logger.debug('P-values calculated')

        else:
            if type(obs_clusters) == float:
                total_clusters = obs_score = float('nan')
            else:
                total_clusters = obs_score = 0

            mean_sim_clusters = median_sim_clusters = std_sim_clusters = mean_sim_score = median_sim_score = \
                std_sim_score = empirical_pvalue = analytical_pvalue = float('nan')

        # Get GCG boolean
        cgc = element.split('_')[0] in self.cgc_genes

        return element, \
            (seq.get_length(self.regions_d[element]), len(self.mutations_d[element]), total_clusters,
                mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, mean_sim_score,
                median_sim_score, std_sim_score, empirical_pvalue, analytical_pvalue, cgc), \
            (obs_clusters, cgc)

    def run(self):
        """
        Analyze elements
        :return:
            elements_results: dict of dict, results of elements analysis
            clusters_results dict of dict, results of clusters analysis
        """
        elements_results = defaultdict()
        clusters_results = defaultdict()

        # Filter elements >= cutoff mutations
        analyzed_elements = []
        belowcutoff_elements = []

        for elem in self.regions_d.keys():
            if len(self.mutations_d[elem]) >= self.element_mutations_cutoff:
                analyzed_elements.append(elem)
            else:
                belowcutoff_elements.append(elem)
        logger.info('Calculating results {} element{}...'.format(len(analyzed_elements),
                                                                 's' if len(analyzed_elements) > 1 else ''))

        # Performance tunning
        pf_num_elements = 100
        pf_num_simulations = 100000

        # Run analysis on each element
        analyzed_elements = list(chunkizator(analyzed_elements, pf_num_elements))
        for i, elements in enumerate(analyzed_elements, start=1):
            observed_clusters_d = {}
            observed_scores_d = {}
            simulations = []
            noclusters_elements = []
            logger.info("Iteration {} of {}".format(i, len(analyzed_elements)))

            for element in elements:
                # Get element pre-smoothing, common for observed and simulated mutations
                logger.debug('Calculating pre-smoothing...')
                probs = self.mut_probabilities(element, kmer=self.kmer)

                # Calculate observed mutations results
                observed_clusters_d[element], observed_scores_d[element] = self.analysis(element,
                                                                                         self.mutations_d[element],
                                                                                         analysis_mode='obs')
                logger.debug('Observed mutations analyzed')

                # Calculate simulated mutations results for elements with observed clusters
                logger.debug('Start simulations...')
                if observed_clusters_d[element]:
                    element_size = ((len(self.mutations_d[element]) * self.n_simulations) // pf_num_simulations) + 1
                    chunk_size = (self.n_simulations // element_size) + 1
                    simulations += [(element, probs, n_sim) for n_sim in partitions_list(self.n_simulations,
                                                                                         chunk_size)]
                else:
                    noclusters_elements.append(element)

            with Pool(max_workers=self.cores) as executor:

                # Process simulations
                sim_scores_list = defaultdict(list)
                sim_clusters_list = defaultdict(list)

                for element, sim_scores_chunk, sim_cluster_chunk in tqdm(
                        executor.map(self.simulate_and_analysis, simulations), total=len(simulations),
                        desc="simulations".rjust(30)
                ):
                    sim_scores_list[element] += sim_scores_chunk
                    sim_clusters_list[element] += sim_cluster_chunk

                # Add information of elements without clusters
                for element in noclusters_elements:
                    sim_scores_list[element] = sim_clusters_list[element] = float('nan')

                # Post process
                post_item_analyzed = [(e, observed_clusters_d[e], observed_scores_d[e],
                                       sim_clusters_list[e], sim_scores_list[e]) for e in elements]
                post_item_nan = [(e, float('nan'), float('nan'), float('nan'), float('nan')) for
                                 e in belowcutoff_elements]
                total_items = post_item_analyzed + post_item_nan

                for e, er, cr in tqdm(map(self.post_process, total_items), total=len(total_items),
                                      desc="post processing".rjust(30)):
                    elements_results[e] = er
                    clusters_results[e] = cr

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
    :param iterable:
    :param size: int, elements in the chunk
    :return: list. Chunk
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
