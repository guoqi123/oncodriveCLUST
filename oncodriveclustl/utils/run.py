# Import modules
import os.path
from concurrent.futures import ProcessPoolExecutor as Pool
import math
from collections import defaultdict
from collections import namedtuple
import pickle

import daiquiri
from intervaltree import IntervalTree
from tqdm import tqdm
import numpy as np
import tabix

from oncodriveclustl.utils import smoothing as smo
from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import analyticalpval as ap
from oncodriveclustl.utils import plots as plot

# Logger
logger = daiquiri.getLogger()
Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample')



class Experiment:
    """Class to analyze elements of a cancer dataset"""

    def __init__(self, regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, genome, path_pickle,
                 element_mutations_cutoff, cluster_mutations_cutoff, smooth_window, cluster_window,
                 cluster_score, element_score, kmer, n_simulations, simulation_mode, simulation_window, cores, seed,
                 conseq, plot):
        """Initialize the Experiment class
        :param regions_d: dict, dictionary of IntervalTrees containing genomic regions from all analyzed elements
        :param cds_d: dictionary of dictionaries with relative cds position of genomic regions if cds is True
        :param chromosomes_d: dict, dictionary containing chromosomes from all analyzed elements
        :param strands_d: dic, dictionary containing strands from all analyzed elements
        :param mutations_d: dictionary, key = element, value = list of mutations formatted as namedtuple(position, sample)
        :param samples_d: dictionary, key = sample, value = number of mutations
        :param genome: genome to use
        :param path_pickle: path to signatures pickle
        :param element_mutations_cutoff: int, cutoff of element mutations
        :param cluster_mutations_cutoff: int, cutoff of cluster mutations
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
        :param conseq: True, use aa consequence type
        :param plot: bool, True generates a clustering plot for an element
        :return: None
        """

        global Mutation
        global Cds

        self.regions_d = regions_d
        self.cds_d = cds_d
        self.chromosomes_d = chromosomes_d
        self.strands_d = strands_d
        self.mutations_d = mutations_d
        self.samples_d = samples_d
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
        self.simulation_window = simulation_window + (1 - simulation_window % 2)
        self.cores = cores
        self.seed = seed
        self.plot = plot
        self.conseq = conseq

        # Read CGC
        if self.genome_build.__name__ == 'hg19':
            # TODO Remove this hardcoded file
            with open(os.path.join(os.path.dirname(__file__), '../data/CGCMay17_cancer_types_TCGA.tsv'), 'r') as fd:
                self.cgc_genes = set([line.split('\t')[0] for line in fd])
        else:
            self.cgc_genes = set()

        # Read vep pickle
        if self.conseq:
            # TODO Remove this hardcoded file
            with open('/home/carnedo/projects/inputs/vep/vep_canonical.pickle', 'rb') as fd:
                self.coding_consequence = pickle.load(fd)
        else:
            self.coding_consequence = {}

    # def __enter__(self):
    #     if self.conseq:
    #         # TODO Remove this hardcoded file
    #         with open('/home/carnedo/projects/inputs/vep/vep_canonical.pickle', 'rb') as fd:
    #             self.coding_consequence = pickle.load(fd)
    #         # self.tb = tabix.open('/workspace/datasets/phd_snp_g/input_files_cds/vep_canonical.tsv.gz')
    #     else:
    #         self.coding_consequence = {}
    #
    #     return self


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

        # TODO check
        except ZeroDivisionError as e:
            logger.error('{}. No mutational probabilities derived from signatures in element {}'.format(e, element))
            return False

    def mut_probabilities(self, element):
        """
        Generate mutational probabilities per position of an element based on the sequence context observed signature
        :param element: element to calculate pre-smoothing
        :return:
            probs_tree: IntervalTree of genomic regions, interval.data are normalized mutational probabilities for each
            position. Length == genomic + simulation window
            skip: bool, if True skip further analysis
        """
        skip = False
        nu = 0
        delta = 1 if self.kmer == 3 else 2
        nucleot = {'A', 'C', 'G', 'T'}
        probs_tree = IntervalTree()

        simulation_window = self.simulation_window
        correction = 1

        if os.path.isfile(self.path_pickle):
            # Read signatures pickle
            signatures = pickle.load(open(self.path_pickle, "rb"))
            signatures = signatures['probabilities']
            logger.debug('Signatures read')

            # Iterate through genomic regions to get their sequences
            for interval in self.regions_d[element]:
                probabilities = []
                start = interval[0] - (simulation_window // 2) - delta
                size = interval[1] - interval[0] + (simulation_window - correction) + delta*2
                # genomic start -d -sw//2, genomic end +d +sw//2
                sequence = self.genome_build(self.chromosomes_d[element], start, size)

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

                # Normalize and add probabilities to probs_tree
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
                    if len(probabilities) != (interval[1] - interval[0] + simulation_window - correction):
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
        Calculate smoothing, clustering and element score for observed or simulated mutations
        :param element: element of analysis.
        :param mutations: list, list of mutations of an element
        :param analysis_mode: str, observed or simulated; default simulated
        :return:
            dict of dicts: cutoff clusters
            int, gene score
        """

        if self.cds_d:
            cds_d = self.cds_d[element]
        else:
            cds_d = {}

        smooth_tree, mutations_in = smo.smooth(element, self.regions_d[element], cds_d, mutations, self.tukey_filter, self.simulation_window)
        index_tree = clu.find_locals(smooth_tree, cds_d)
        raw_clusters_tree = clu.raw_clusters(index_tree)
        merge_clusters_tree = clu.merge_clusters(raw_clusters_tree, self.cluster_window)
        filter_clusters_tree = clu.clusters_mut(merge_clusters_tree, mutations_in, self.cluster_mutations_cutoff)
        score_clusters_tree = clu.fmutations_score(filter_clusters_tree, self.regions_d[element], len(mutations_in))
        logger.debug('Clusters scores calculated')
        element_score = score.element_score(score_clusters_tree, analysis_mode, self.element_score)
        logger.debug('Element score calculated')

        if self.plot:
            logger.info('Generating plot: {}'.format(element))
            return smooth_tree, raw_clusters_tree, merge_clusters_tree, score_clusters_tree, element_score
        else:
            return score_clusters_tree, element_score

    def simulate_and_analysis(self, item):
        """
        Simulate mutations and analyze simulations
        :param item: tuple, element of analysis data
        :return:
            element: str, element of analysis
            sim_scores_chunk: list, simulated element's results
            sim_cluster_chunk: list, simulated cluster's results
        """
        element, probs_tree, n_sim = item
        id = element.split('_')[1]
        sim_scores_chunk = []
        sim_cluster_chunk = []
        df = []
        half_window = (self.simulation_window - 1) // 2

        delta = 1 if self.kmer == 3 else 2
        nucleot = {'A', 'C', 'G', 'T'}
        signatures = pickle.load(open(self.path_pickle, "rb"))
        signatures = signatures['probabilities']

        # Simulate mutations
        for mutation in self.mutations_d[element]:
            # Get hotspot for simulations
            expected_hotspot_begin = mutation.position - half_window
            expected_hotspot_end = mutation.position + half_window

            if self.simulation_mode == 'exon_restricted':
                # Check if hospot outside region
                check_5 = expected_hotspot_begin < mutation.region[0]
                check_3 = expected_hotspot_end > (mutation.region[1] - 1)

                if check_5 and check_3:
                    hotspot_begin = mutation.region[0]
                    hotspot_end = mutation.region[1] - 1  # regions end +1 in tree
                elif check_5:
                    hotspot_begin = mutation.region[0]
                    hotspot_end = mutation.region[0] + self.simulation_window - 1  # window //2 per side
                elif check_3:
                    hotspot_end = mutation.region[1] - 1  # regions end +1 in tree
                    hotspot_begin = hotspot_end - self.simulation_window + 1  # window //2 per side
                else:
                    hotspot_begin = expected_hotspot_begin
                    hotspot_end = expected_hotspot_end
            else:
                hotspot_begin = expected_hotspot_begin
                hotspot_end = expected_hotspot_end

            hotspot = tuple([hotspot_begin, hotspot_end])
            start_index = hotspot[0] - (mutation.region[0] - half_window)
            end_index = hotspot[1] - (mutation.region[0] - half_window) + 1  # +1 because it is a slice
            for interval in probs_tree[mutation.region[0]]:  # unique iteration
                simulations = np.random.choice(range(hotspot[0], hotspot[1] + 1), size=n_sim,
                                               p=self.normalize(element, interval.data[start_index:end_index]))

                # TODO improve
                # Add simulations
                l = []
                if self.coding_consequence:
                    print('got here')
                    quit()
                    for count, pos in enumerate(simulations):
                        # Get alternate for simulated mutation
                        probs_alternates = []
                        changes = []
                        start = mutation.position - delta
                        size = delta * 2 + 1
                        ref_kmer = self.genome_build(self.chromosomes_d[element], start, size)
                        for alt in nucleot.difference({ref_kmer[delta]}):  # mutational prob to any other kmer
                            alt_kmer = ref_kmer[: self.kmer // 2] + alt + ref_kmer[self.kmer // 2 + 1:]
                            probs_alternates.append(signatures[(ref_kmer, alt_kmer)])
                            changes.append(alt)
                        alternate = np.random.choice(changes, size=1, p=self.normalize(element, probs_alternates))

                        # Get consequence
                        ensid = element.split('_')[1]
                        try:
                            muttype = 0 if pos in self.coding_consequence[ensid][alternate] else 1
                        except:
                            muttype = 1

                        # consequences = [
                        #     c[8] for c in self.tb.query(self.chromosomes_d[element], pos - 1, pos) if c[4] == alternate
                        # ]
                        # muttype = 0 if all([i == 'synonymous_variant' for i in consequences]) else 1
                        l.append(Mutation(pos, mutation.region, alternate, muttype, mutation.sample))
                        df.append(l)

                    print(df)
                    quit()

                else:
                    muttype = 1
                    alt = 'N'
                    for count, pos in enumerate(simulations):
                        l.append(Mutation(pos, mutation.region, alt, muttype, mutation.sample))
                    df.append(l)

        # Start analysis
        for simulated_mutations in zip(*df):
            cutoff_clusters, element_score = self.analysis(element, simulated_mutations)
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

    def length(self, element):
        """Calculate length of an element (sum of input regions)
        :param element: element to analyze
        :return length: int, length of an element
        """
        length_ele = 0
        for interval in self.regions_d[element]:
            length_ele += (interval[1] - interval[0])
        return length_ele

    def post_process(self, items):
        """
        Post processing of results
        :param item: tuple, elements results
        :return:
        """

        results = []
        for item in items:
            element, obs_clusters, obs_score, sim_clusters_list, sim_scores_list = item

            # If analyzed element and element has clusters
            if type(obs_clusters) != float:
                if type(sim_scores_list) != float:
                    if sum(sim_scores_list) == 0:
                        logger.warning('No simulated clusters in {}. '
                                       'Observed cluster p-values calculated with pseudocount'.format(element))
                        n_clusters = 0
                        # Get false p-value
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_scores_list)
                        pseudo_pvalue = obj.calculate(obs_score)
                        for interval in obs_clusters:
                            for cluster, values in interval.data.items():
                                interval.data[cluster]['p'] = pseudo_pvalue
                                n_clusters += 1
                        empirical_pvalue = analytical_pvalue = top_cluster_pvalue = pseudo_pvalue
                        n_clusters_sim = False

                    else:
                        # Element score empirical p-value
                        empirical_pvalue = self.empirical_pvalue(obs_score, sim_scores_list)

                        # Element score analytical p-value
                        sim_scores_array_1000 = np.random.choice(sim_scores_list, size=1000, replace=False)
                        obj = ap.AnalyticalPvalue()

                        obj.calculate_bandwidth(sim_scores_array_1000)
                        analytical_pvalue = obj.calculate(obs_score)

                        # Cluster analytical p-values
                        sim_clusters_score = []
                        obs_clusters_score = []
                        # for simulation in sim_clusters_list:
                        #     for interval in simulation:
                        #         for cluster, values in interval.data.items():
                        #             sim_clusters_score.append(values['score'])
                        #
                        # obj = ap.AnalyticalPvalue()
                        # obj.calculate_bandwidth(sim_clusters_score)

                        for interval in obs_clusters:
                            for cluster, values in interval.data.items():
                                obs_clusters_score.append(values['score'])
                                # cluster_p_value = obj.calculate(values['score'])
                                # interval.data[cluster]['p'] = cluster_p_value
                                interval.data[cluster]['p'] = float('nan')

                        n_clusters = len(obs_clusters_score)
                        n_clusters_sim = True

                        # Element top cluster analytical p-value
                        #top_cluster_pvalue = obj.calculate(max(obs_clusters_score))
                        top_cluster_pvalue = float('nan')
                        logger.debug('P-values calculated')

                else:
                    n_clusters = obs_score = 0
                    n_clusters_sim = float('nan')
                    empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')
            else:
                n_clusters = n_clusters_sim = obs_score = float('nan')

                if len(self.mutations_d[element]) == 1:
                    empirical_pvalue = analytical_pvalue = top_cluster_pvalue = 1
                else:
                    empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')

            # Get GCG boolean
            cgc = element.split('_')[0] in self.cgc_genes

            # Calculate length
            element_length = self.length(element)

            results.append((

                element,

                (self.chromosomes_d[element], self.strands_d[element], element_length, len(self.mutations_d[element]),
                 n_clusters, n_clusters_sim, obs_score, empirical_pvalue, analytical_pvalue, top_cluster_pvalue, cgc),

                (obs_clusters, self.chromosomes_d[element], self.strands_d[element], cgc)
            ))

        return results

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
        if len(analyzed_elements) == 0:
            logger.critical('There are no element above cutoff element mutations to analyze')
            quit()
        else:
            logger.info('Calculating results {} element{}...'.format(len(analyzed_elements),
                        's' if len(analyzed_elements) > 1 else ''))
        # Plot
        if self.plot:
            for element in analyzed_elements:
                smooth_tree, raw_clusters_tree, merge_clusters_tree, score_clusters_tree, element_score = \
                    self.analysis(element, self.mutations_d[element], analysis_mode='obs')
                plot.run_plot(element, self.mutations_d[element], self.cds_d[element],
                              self.strands_d[element], self.chromosomes_d[element], self.smooth_window,
                              smooth_tree, raw_clusters_tree, merge_clusters_tree, score_clusters_tree, element_score)
                logger.info('Plots calculated: {}'.format(element))
            logger.info('Finished')
            quit()

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

                    # TODO increase simulations for elements without simulated clusters

                # Add information of elements without clusters
                for element in nocluster_elements:
                    sim_scores_list[element] = sim_clusters_list[element] = float('nan')

                # Post process
                post_item_simulated = [(e, observed_clusters_d[e], observed_scores_d[e], sim_clusters_list[e],
                                        sim_scores_list[e]) for e in simulated_elements + nocluster_elements]
                post_item_nan = [(e, float('nan'), float('nan'), float('nan'), float('nan')) for
                                 e in noprobabilities_elements + belowcutoff_elements]

                total_items_split = list(chunkizator(post_item_simulated, int(math.ceil(len(post_item_simulated) / (self.cores-1)))))
                total_items_split.append(post_item_nan)
                for results in tqdm(executor.map(self.post_process, total_items_split), total=self.cores,
                                      desc="post processing".rjust(30)):
                    for e, er, cr in results:
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
