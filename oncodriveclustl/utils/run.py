# Import modules
import os.path
from concurrent.futures import ProcessPoolExecutor as Pool

import daiquiri
import pickle
from collections import defaultdict
from intervaltree import IntervalTree
from tqdm import tqdm
import numpy as np

from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import sequence as seq
from oncodriveclustl.utils import analyticalpval as ap

# Logger
logger = daiquiri.getLogger()


class Experiment:
    """Class to analyze elements of a cancer dataset"""

    def __init__(self, regions_d, chromosomes_d, strands_d, mutations_d, genome, path_pickle,
                 element_mutations_cutoff, cluster_mutations_cutoff, cds, smooth_window, cluster_window,
                 cluster_score, element_score, kmer, n_simulations, simulation_mode, simulation_window, cores, seed):
        """Initialize the Experiment class
        :param regions_d: dict, dictionary of IntervalTrees containing genomic regions from all analyzed elements
        :param chromosomes_d: dict, dictionary containing chromosomes from all analyzed elements
        :param strands_d: dic, dictionary containing strands from all analyzed elements
        :param mutations_d: dict, dictionary containing mutations lists from all analyzed elements
        :param genome: genome to use
        :param element_mutations_cutoff: int, cutoff of element mutations
        :param cluster_mutations_cutoff: int, cutoff of cluster mutations
        :param cds: bool, True calculates clustering on cds
        :param smooth_window: int, smoothing window
        :param cluster_window: int, clustering window
        :param cluster_score: cluster score method
        :param element_score: element score method
        :param kmer: int, number of nucleotides of the signature
        :param n_simulations: int, number of simulations
        :param simulation_mode: str, simulation mode
        :param simulation_window: int, window to simulate mutations in hotspot mode
        :param cores: int, number of CPUs to use
        :param seed: int, seed
        :return: None
        """

        self.regions_d = regions_d
        self.chromosomes_d = chromosomes_d
        self.strands_d = strands_d
        self.mutations_d = mutations_d
        self.genome_build = self.load_genome(genome)
        self.path_pickle = path_pickle
        self.element_mutations_cutoff = element_mutations_cutoff
        self.cluster_mutations_cutoff = cluster_mutations_cutoff
        self.cds = cds
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

        # TODO check what happens here
        except ZeroDivisionError as e:
            logger.error('{}. No mutational probabilities derived from signatures in element {}'.format(e, element))
            return False

    def mut_probabilities(self, element):
        """
        Generate mutational probabilities per position of an element based on the sequence context observed signature
        :param element: element to calculate pre-smoothing
        :return:
            probs_tree: IntervalTree of genomic regions, interval.data are mutational probabilities per position
            (not normalized). Length == genomic + simulation window
            skip: bool, if True skip further analysis
        """
        skip = False
        nu = 0
        sequ = 0
        delta = 1 if self.kmer == 3 else 2
        nucleot = {'A', 'C', 'G', 'T'}
        probs_tree = IntervalTree()

        if os.path.isfile(self.path_pickle):
            # Read signatures pickle
            signatures = pickle.load(open(self.path_pickle, "rb"))
            signatures = signatures['probabilities']
            logger.debug('Signatures read')

            # Iterate through tuples of coordinates, get genomic regions sequence
            for interval in self.regions_d[element]:
                probabilities = []
                start = interval[0] - (self.simulation_window // 2) - delta
                end = interval[1] - interval[0] + 1 + self.simulation_window + delta*2
                # genomic start -d -sw//2, genomic end +d +sw//2
                sequence = self.genome_build(self.chromosomes_d[element], start, end)
                sequ += len(sequence)-2

                # Search kmer probabilities
                for n in range(delta, len(sequence)-delta):  # start to end
                    nu += 1
                    ref_kmer = sequence[n - delta: n + delta + 1]
                    n_number = ref_kmer.count('N')
                    if n_number == 0:
                        prob = 0
                        for alt in nucleot.difference({ref_kmer[self.kmer//2]}):  # mutational prob to any other kmer
                            alt_kmer = ref_kmer[: self.kmer//2] + alt + ref_kmer[self.kmer//2 + 1:]
                            prob += signatures[(ref_kmer, alt_kmer)]
                    else:
                        prob = 0
                    # Append to probabilities
                    probabilities.append(prob)
                # Add probabilities to probs_tree
                probs_tree.addi(interval[0], interval[1], probabilities)

            # Check probabilities
            check = 0
            for interval in probs_tree:
                probabilities = interval.data
                if probabilities:
                    check = 1
                    if sum(probabilities) == 0:
                        logger.critical('All context based mutational probabilities per position in {} equal to 0\n'
                                        '{} analysis is skipped'.format(element, element))
                        skip = True
                    if len(probabilities) != (interval[1] - interval[0] + 1 + self.simulation_window):
                        logger.warning('{} probabilities list length is different than expected, '
                                       'please review results'.format(element))
                        skip = True
            if check == 0:
                logger.critical('Context based mutational probabilities were not calculated for {}\n'
                                '{} analysis is skipped'.format(element, element))
                skip = True

        return probs_tree, skip

    def analysis(self, element, mutations, analysis_mode='sim'):
        """
        # TODO
        Calculate smoothing, clustering and element score for observed or simulated mutations
        :param element: element of analysis.
        :param mutations: list, list of mutations of an element
        :param analysis_mode: str, observed or simulated; default simulated
        :return:
            dict of dicts: cutoff clusters
            int, gene score
        """

        index_tree, length_tree, smooth_tree = clu.find_locals(self.regions_d[element], mutations,
                                        self.tukey_filter, self.cds, self.simulation_window)
        raw_clusters_tree = clu.raw_clusters(index_tree)
        merge_clusters_tree = clu.merge_clusters(raw_clusters_tree, self.cluster_window)
        clusters_mut_tree = clu.clusters_mut(merge_clusters_tree, length_tree, mutations, self.cds)

        score_clusters_tree = clu.score_clusters(clusters_mut_tree, length_tree, self.regions_d[element], self.cluster_mutations_cutoff,
                                                 self.cluster_score, self.cds, len(mutations))
        logger.debug('Clusters calculated')
        element_score = score.element_score(score_clusters_tree, analysis_mode, self.element_score)
        logger.debug('Element score calculated')

        # for i in score_clusters_tree:
        #     print('Interval', i[0], i[1])
        #     for k, v in i.data.items():
        #         print('Cluster number', k)
        #         for k2, v2 in v.items():
        #             print(k2, v2)

        return score_clusters_tree, element_score

    def simulate_and_analysis(self, item):
        """
        # TODO
        Simulate mutations and analyze simulations
        :param item: tuple, element of analysis data
        :return:
            element: str, element of analysis
            sim_scores_chunk: list, simulated element's results
            sim_cluster_chunk: list, simulated cluster's results
        """
        element, probs_tree, n_sim = item
        sim_scores_chunk = []
        sim_cluster_chunk = []
        simulation_hotspots = {}
        mutations = dict((m, self.mutations_d[element].count(m)) for m in self.mutations_d[element])

        if self.simulation_mode == 'hotspot':
            # Get half window
            half_window = self.simulation_window // 2
            # Get genomic coordinates of simulation hotspot inside a region
            for mutation, count in mutations.items():
                # Find region containing observed mutation
                for interval in self.regions_d[element][mutation]:  # unique
                    # Get margins of hotspot
                    simulation_hotspots[tuple([mutation - half_window, mutation + half_window])] = (count,
                                                                                                    interval.begin)
            # Simulate mutations
            initializer = 0
            for hotspot, values in simulation_hotspots.items():
                region_start = values[1] - half_window
                start_to_0 = hotspot[0] - region_start
                start_to_1 = hotspot[1] - region_start
                for interval in probs_tree[values[1]]:  # unique
                    mut_hotspots = np.random.choice(range(hotspot[0], hotspot[1]), size=(n_sim, values[0]),
                                                p=self.normalize(element, interval.data[start_to_0:start_to_1]))
                if initializer != 0:
                    simulated_mutations = np.hstack([simulated_mutations, mut_hotspots])
                else:
                    simulated_mutations = mut_hotspots
                initializer = 1

        # Start analysis
        for i in range(n_sim):
            cutoff_clusters, element_score = self.analysis(element, simulated_mutations[i])
            sim_scores_chunk.append(element_score)
            sim_cluster_chunk.append(cutoff_clusters)

        return element, sim_scores_chunk, sim_cluster_chunk

    def simulate_and_analysis_old(self, item):
        """
        # TODO
        Simulate mutations and analyze simulations
        :param item: tuple, element of analysis data
        :return:
            element: str, element of analysis
            sim_scores_chunk: list, simulated element's results
            sim_cluster_chunk: list, simulated cluster's results
        """
        element, probs_tree, n_sim = item
        sim_scores_chunk = []
        sim_cluster_chunk = []

        if self.simulation_mode == 'hotspot':
            simulation_hotspots = {}
            # Mutations
            mutations = dict((m, self.mutations_d[element].count(m)) for m in self.mutations_d[element])
            # Get half window
            half_window = self.simulation_window // 2

            if self.cds is False:
                # Get genomic coordinates of simulation hotspot inside a region
                for mutation, count in mutations.items():
                    # Find region containing observed mutation
                    for interval in self.regions_d[element][mutation]:  # unique
                        # Get margins of hotspot
                        a = mutation - half_window
                        b = mutation + half_window
                        if a < interval.begin:
                            a = interval.begin
                        if b > interval.end:
                            b = interval.end
                        simulation_hotspots[tuple([a, b])] = (count, interval.begin)
                # Simulate mutations in hotspots
                initializer = 0
                n = 0
                for hotspot, values in simulation_hotspots.items():
                    region_start = values[1]
                    a_to_start = hotspot[0] - region_start
                    b_to_start = hotspot[1] - region_start
                    for interval in probs_tree[values[1]]:  # unique
                        n += values[0]
                        mut_hotspots = np.random.choice(range(hotspot[0], hotspot[1]), size=(n_sim, values[0]),
                                                        p=self.normalize(element, interval.data[a_to_start:b_to_start]))
                    if initializer != 0:
                        simulated_mutations = np.hstack([simulated_mutations, mut_hotspots])
                    else:
                        simulated_mutations = mut_hotspots
                    initializer = 1
            else:

                # TODO create this tree somewhere else
                # Find regions
                length_tree = IntervalTree()
                length_d = defaultdict()
                start = 0
                for interval in sorted(self.regions_d[element]):
                    length = interval.end - interval.begin
                    length_tree.addi(start, start + length, interval[0])
                    length_d[interval[0]] = (start, start + length)
                    start = start + length + 1

                print(length_d)
                quit()

                for mutation, count in mutations.items():
                    # Find region containing observed mutation
                    for interval in self.regions_d[element][mutation]:  # unique
                        # Get margins of hotspot
                        a = mutation - half_window
                        b = mutation + half_window

                        if a < interval.begin:
                            a = interval.begin

                        if b > interval.end:
                            b = interval.end

                        simulation_hotspots[tuple([a, b])] = (count, interval.begin)
                # Simulate mutations in hotspots
                initializer = 0
                n = 0
                for hotspot, values in simulation_hotspots.items():
                    region_start = values[1]
                    a_to_start = hotspot[0] - region_start
                    b_to_start = hotspot[1] - region_start
                    for interval in probs_tree[values[1]]:  # unique
                        n += values[0]
                        mut_hotspots = np.random.choice(range(hotspot[0], hotspot[1]), size=(n_sim, values[0]),
                                                        p=self.normalize(element, interval.data[a_to_start:b_to_start]))
                    if initializer != 0:
                        simulated_mutations = np.hstack([simulated_mutations, mut_hotspots])
                    else:
                        simulated_mutations = mut_hotspots
                    initializer = 1

        # Start analysis
        for i in range(n_sim):
            cutoff_clusters, element_score = self.analysis(element, simulated_mutations[i])
            sim_scores_chunk.append(element_score)
            sim_cluster_chunk.append(cutoff_clusters)

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

        # If analyzed element and element has clusters
        if type(obs_clusters) != float:
            if type(sim_scores_list) != float:

                # Element score empirical p-value
                empirical_pvalue = self.empirical_pvalue(obs_score, sim_scores_list)

                # Element score analytical p-value
                sim_scores_list_1000 = np.random.choice(sim_scores_list, size=1000, replace=False)
                analytical_pvalue, band = ap.AnalyticalPvalue().calculate(obs_score, sim_scores_list_1000)

                # Cluster analytical p-values
                sim_clusters_score = []
                obs_clusters_score = []
                for simulation in sim_clusters_list:
                    for interval in simulation:
                        for cluster, values in interval.data.items():
                            sim_clusters_score.append(values['score'])
                for interval in obs_clusters:
                    for cluster, values in interval.data.items():
                        obs_clusters_score.append(values['score'])
                        cluster_p_value, band = ap.AnalyticalPvalue().calculate(values['score'], sim_clusters_score)
                        interval.data[cluster]['p'] = cluster_p_value
                n_clusters = len(obs_clusters_score)

                # Element top cluster analytical p-value
                top_cluster_pvalue, band = ap.AnalyticalPvalue().calculate(max(obs_clusters_score), sim_clusters_score)
                logger.debug('P-values calculated')

            else:
                n_clusters = obs_score = 0
                empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')
        else:
            n_clusters = obs_score = float('nan')
            empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')

        # Get GCG boolean
        cgc = element.split('_')[0] in self.cgc_genes

        # Calculate length
        element_length = 0
        for interval in self.regions_d[element]:
            element_length += (interval[1] - interval[0] + 1)

        return element, \
            (self.chromosomes_d[element], self.strands_d[element], element_length, len(self.mutations_d[element]),
             n_clusters, obs_score, empirical_pvalue, analytical_pvalue, top_cluster_pvalue, cgc), \
            (obs_clusters, self.chromosomes_d[element], self.strands_d[element], cgc)

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
            simulated_elements = []
            noprobabilities_elements = []
            nocluster_elements = []
            logger.info("Iteration {} of {}".format(i, len(analyzed_elements)))

            for element in elements:
                # Calculate observed results

                observed_clusters_d[element], observed_scores_d[element] = self.analysis(element,
                                                                                         self.mutations_d[element],
                                                                                         analysis_mode='obs')
                logger.debug('Observed mutations analyzed')
                # Calculate simulated mutations results if element has observed clusters at least in one region
                check = 0
                for interval in observed_clusters_d[element]:
                    if interval.data:
                        check = 1
                    if check == 1:
                        break
                # If there is at least a region with clusters
                if check != 0:
                    # Check if probabilities can be calculated
                    logger.debug('Calculating context based mutational '
                                 'probabilities for region {} {}...'.format(interval[0], interval[1]))
                    probs_tree, skip = self.mut_probabilities(element)
                    if not skip:
                        element_size = ((len(self.mutations_d[element]) * self.n_simulations) // pf_num_simulations) + 1
                        chunk_size = (self.n_simulations // element_size) + 1
                        simulations += [(element, probs_tree, n_sim) for n_sim in partitions_list(self.n_simulations,
                                        chunk_size)]
                        simulated_elements.append(element)
                    # Skip analysis if problems with mutational probabilities
                    else:
                        noprobabilities_elements.append(element)
                else:
                    nocluster_elements.append(element)

            with Pool(max_workers=self.cores) as executor:
                # Process simulations
                sim_scores_list = defaultdict(list)
                sim_clusters_list = defaultdict(list)

                # Analyze elements with observed clusters
                for element, sim_scores_chunk, sim_cluster_chunk in tqdm(
                        executor.map(self.simulate_and_analysis, simulations), total=len(simulations),
                        desc="simulations".rjust(30)):

                    sim_scores_list[element] += sim_scores_chunk
                    sim_clusters_list[element] += sim_cluster_chunk

                # Add information of elements without clusters
                for element in nocluster_elements:
                    sim_scores_list[element] = sim_clusters_list[element] = float('nan')

                # Post process
                post_item_simulated = [(e, observed_clusters_d[e], observed_scores_d[e], sim_clusters_list[e],
                                        sim_scores_list[e]) for e in simulated_elements + nocluster_elements]
                post_item_nan = [(e, float('nan'), float('nan'), float('nan'), float('nan')) for
                                 e in noprobabilities_elements + belowcutoff_elements]
                total_items = post_item_simulated + post_item_nan

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
