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
import bgdata as bgd
import bgreference as bgr

from oncodriveclustl.utils import smoothing as smo
from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import analyticalpval as ap
from oncodriveclustl.utils import plots as plot


# Logger
logger = daiquiri.getLogger()
Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample, cancertype')


class Experiment:
    """Class to analyze elements of a cancer dataset"""

    def __init__(self,
                 regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, genome,
                 path_cache, cohorts_d,
                 element_mutations_cutoff, cluster_mutations_cutoff,
                 smooth_window, cluster_window,
                 cluster_score, element_score,
                 kmer, n_simulations, simulation_mode, simulation_window,
                 cores,
                 conseq, protein,
                 plot):
        """Initialize the Experiment class
        :param regions_d: dict, dictionary of IntervalTrees containing genomic regions from all analyzed elements
        :param cds_d: dictionary of dictionaries with relative cds position of genomic regions if cds is True
        :param chromosomes_d: dict, dictionary containing chromosomes from all analyzed elements
        :param strands_d: dic, dictionary containing strands from all analyzed elements
        :param mutations_d: dictionary, key = element, value = list of mutations as namedtuple(position, sample)
        :param samples_d: dictionary, key = sample, value = number of mutations
        :param genome: str, genome to use
        :param path_cache: str, path to pickles directory (cache)
        :param cohorts_d: dictionary, key = element, value = set of cohorts with element mutations
        :param element_mutations_cutoff: int, cutoff of element mutations
        :param cluster_mutations_cutoff: int, cutoff of cluster mutations
        :param smooth_window: int, smoothing window
        :param cluster_window: int, clustering window
        :param cluster_score: str, cluster score method
        :param element_score: str, element score method
        :param kmer: int, number of nucleotides of the signature
        :param n_simulations: int, number of simulations
        :param simulation_mode: str, simulation mode
        :param simulation_window: int, window to simulate mutations in hotspot mode
        :param cores: int, number of CPUs to use
        :param conseq: True, use aa consequence type
        :param protein: bool, True analyzes clustering in translated protein sequences
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
        self.genome = genome
        self.path_cache = path_cache
        self.cohorts_d = cohorts_d
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
        self.plot = plot
        self.conseq = conseq
        self.protein = protein


        # Read CGC
        if self.genome == 'hg19':
            # TODO Remove this hardcoded file
            with open(os.path.join(os.path.dirname(__file__), '../data/CGCMay17_cancer_types_TCGA.tsv'), 'r') as fd:
                self.cgc_genes = set([line.split('\t')[0] for line in fd])
        else:
            self.cgc_genes = set()

        if self.conseq:
            self.conseq_path = bgd.get_path('oncodriveclustl', 'vep88', 'hg19_canonical_conseq')
        else:
            self.conseq_path = None

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
            logger.error('{}. No mutational probabilities derived from signatures in element {}'.format(e, element))
            return False

    def mut_probabilities(self, element):
        """
        Generate mutational probabilities per position of an element using the sequence context observed mutational
        probabilities calculated from the input cohort/s.
        :param element: element to calculate pre-smoothing
        :return:
            probs_tree: IntervalTree of genomic regions. Length == 3*(genomic + simulation window)
            conseq_tree: IntervalTree of genomic regions, interval.data are boolean arrays representing synonymous or
            non-synonymous changes of each alternate per position. Length == 3*(genomic + simulation window). Tree empty
            if self.conseq is False.
            skip: bool, if True skip further analysis
        """
        skip = False
        nu = 0
        delta = 1 if self.kmer == 3 else 2
        nucleot = {'A', 'C', 'G', 'T'}
        probs_tree = defaultdict(IntervalTree)
        conseq_tree = IntervalTree()
        simulation_window = self.simulation_window
        correction = 1   # (simw // 2 ) * 2 == simw - 1

        if self.conseq:
            ensid = element.split('//')[1]
            path_to_vep_pickle = self.conseq_path + '/{}.pickle'.format(ensid)
            try:
                with open(path_to_vep_pickle, 'rb') as fd:
                    conseq_d = pickle.load(fd)
            except FileNotFoundError:
                skip = True

        # Check signatures pickles exist and read
        # TODO to much memory?
        signatures_d = defaultdict()
        for cohort in self.cohorts_d[element]:
            if os.path.isfile(self.path_cache):
                path_to_signature_pickle = self.path_cache
            else:
                path_to_signature_pickle = os.path.join(self.path_cache, '{}_kmer_{}.pickle'.format(cohort, self.kmer))

            if os.path.isfile(path_to_signature_pickle):
                signature = pickle.load(open(path_to_signature_pickle, "rb"))
                signatures_d[cohort] = signature['probabilities']
            else:
                skip = True

        if not skip:
            # Iterate through genomic regions to get their sequences
            for interval in self.regions_d[element]:
                probabilities = defaultdict(list)
                consequences = []
                start = interval[0] - (simulation_window // 2) - delta
                size = interval[1] - interval[0] + (simulation_window - correction) + delta*2
                # genomic start -d -sw//2, genomic end +d +sw//2
                try:
                    sequence = bgr.refseq(self.genome, self.chromosomes_d[element], start, size)
                except ValueError as e:
                    logger.error(e, element, start, size, interval[0], interval[1])

                # Search kmer probabilities
                for n in range(delta, len(sequence)-delta):  # start to end
                    position = start + n
                    nu += 1
                    ref_kmer = sequence[n - delta: n + delta + 1]
                    prob = defaultdict(list)
                    conseq = []
                    if ref_kmer.count('N') == 0:
                        # calculate mutational prob to any other kmer
                        # sort alternates to keep track of them
                        for alt in sorted(list(nucleot.difference({ref_kmer[self.kmer//2]}))):
                            alt_kmer = ref_kmer[: self.kmer//2] + alt + ref_kmer[self.kmer//2 + 1:]
                            for cohort, signature in signatures_d.items():
                                prob[cohort].append(signature[(ref_kmer, alt_kmer)])
                            if self.conseq:
                                # Get consequence
                                con = 0 if position in conseq_d.get(alt, []) else 1
                                conseq.append(con)
                            else:
                                conseq.append(1)
                    else:
                        for cohort, signature in signatures_d.items():
                            prob[cohort].extend([0,0,0])
                        conseq = [0, 0, 0]  # TODO check, give conseq synonymous to pos where ref_kmer has N

                    # Extend position info
                    for cohort in signatures_d.keys():
                        probabilities[cohort].extend(prob[cohort])
                    consequences.extend(conseq)

                # Add to tree
                for cohort in signatures_d.keys():
                    probs_tree[cohort].addi(interval[0], interval[1], probabilities[cohort])
                conseq_tree.addi(interval[0], interval[1], consequences)

            # Check
            skip = False
            for cohort, tree in probs_tree.items():
                for interval in tree:
                    probabilities = interval.data
                    if probabilities:
                        if sum(probabilities) == 0:
                            logger.critical('All context based mutational probabilities per alternate in {} equal to 0\n'
                                            '{} analysis is skipped'.format(element, element))
                            skip = True
                            break
                        if len(probabilities) != 3 * (interval[1] - interval[0] + simulation_window - correction):
                            logger.warning('{} probabilities list length is different than expected, '
                                           'please review results'.format(element))
                            skip = True
                            break

                        if len(probabilities) != len(list(conseq_tree[interval[0]])[0][-1]):
                            logger.warning('{} probabilities list length is different than expected, '
                                           'please review results'.format(element))
                            skip = True
                            break
                if skip:
                    logger.critical('Context based mutational probabilities were not calculated for {}\n'
                                    '{} analysis is skipped'.format(element, element))
                    break

        return probs_tree, conseq_tree, skip

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

        if not self.protein:
            smooth_tree, mutations_in = smo.smooth_nucleotide(self.regions_d[element], cds_d, mutations,
                                                              self.tukey_filter, self.simulation_window)
        else:
            # protein_d_path = '/home/carnedo/projects/oncodriveclustl/oncodriveclustl/outputs/transcript/cache/protein_sequences.pickle'
            # with open(protein_d_path, 'rb') as fd:
            #     protein_d = pickle.load(fd)
            # if self.strands_d[element] == '-':
            #     protein_seq = protein_d[element][::-1]
            # else:
            #     protein_seq = protein_d[element]
            smooth_tree, mutations_in = smo.smooth_aminoacid(self.regions_d[element], self.chromosomes_d[element],
                                                             self.strands_d[element],
                                                             self.genome, self.tukey_filter, cds_d, mutations)
            cds_d = {}  # Next functions performed with cds False

        index_tree = clu.find_locals(smooth_tree, cds_d)

        # for i in smooth_tree:
        #     print(element, i[0], i[1])
        #     print(len(i.data), i.data[-14:-1])
        #
        # for i in index_tree:
        #     print(element, i[0], i[1])
        #     for index in i.data:
        #         print('\t', index)

        raw_clusters_tree = clu.raw_clusters(index_tree)

        # for i in raw_clusters_tree:
        #     print(element, i[0], i[1])
        #     for c, v in i.data.items():
        #         print('\t', c, v)

        merge_clusters_tree = clu.merge_clusters(raw_clusters_tree, self.cluster_window)
        filter_clusters_tree = clu.clusters_mut(merge_clusters_tree, mutations_in, self.cluster_mutations_cutoff)

        # print('-------------------------------------------------------------------------------------')

        # for i in filter_clusters_tree:
        #     print(element, i[0], i[1])
        #     for c, v in i.data.items():
        #         print('\t', c, v['max'], v['left_m'], v['right_m'], len(v['mutations']))
        #         # print('\t', 'cluster:', c, 'data:', i[1]-v['left_m'][1], i[1]-v['max'][1], i[1]-v['right_m'][1], len(v['mutations']))

        score_clusters_tree = clu.fmutations_score(filter_clusters_tree, self.regions_d[element], len(mutations_in), self.protein)
        logger.debug('Clusters scores calculated')
        element_score = score.element_score(score_clusters_tree, analysis_mode, self.element_score)
        logger.debug('Element score calculated')

        # for i in score_clusters_tree:
        #     print(element, i[0], i[1])
        #     for c, v in i.data.items():
        #         print('\t', c, v['max'], v['left_m'], v['right_m'], len(v['mutations']))
        #         print('\t', 'cluster:', c, 'data:', i[1]-v['left_m'][1], i[1]-v['max'][1], i[1]-v['right_m'][1], len(v['mutations']), 'score:', v['score'])
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
        element, probs_tree, conseq_tree, n_sim = item
        sim_scores_chunk = []
        sim_cluster_chunk = []
        df = []
        half_window = (self.simulation_window - 1) // 2
        nucleot = {'A', 'C', 'G', 'T'}

        # Simulate mutations
        for mutation in self.mutations_d[element]:
            # Get hotspot for simulations
            expected_hotspot_begin = mutation.position - half_window
            expected_hotspot_end = mutation.position + half_window

            if self.simulation_mode == 'region_restricted':
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

            # Map to index
            # 3* accounts for alternates in the array of probabilities
            start_index = 3*(hotspot_begin - (mutation.region[0] - half_window))
            end_index = 3*(hotspot_end - (mutation.region[0] - half_window) + 1)  # +1 because it is a range and slice
            for interval in probs_tree[mutation.cancertype][mutation.region[0]]:  # unique iteration
                simulations = np.random.choice(range(start_index, end_index), size=n_sim,
                                               p=self.normalize(element, interval.data[start_index:end_index]))
                # Add info per simulated mutation
                l = []
                for count, index in enumerate(simulations):
                    # print(count, index, index / 3)
                    muttype = list(conseq_tree[interval[0]])[0][2][index]
                    position = mutation.region[0] - half_window + index // 3
                    ref_nucleotide = bgr.refseq(self.genome, self.chromosomes_d[element], position, 1)
                    # print(position, ref_nucleotide)
                    # Calculate sorted alternates and obtain simulated alternated from index
                    if round(index / 3, 1) == (0.7 + index // 3):
                        alternate_index = 2
                    else:
                        alternate_index = 1 if round(index / 3, 1) == (0.3 + index // 3) else 0

                    alternate = sorted(list(nucleot.difference({ref_nucleotide})))[alternate_index]
                    # Simulated mutation
                    l.append(Mutation(position, mutation.region, alternate, muttype, mutation.sample, mutation.cancertype))
                    # print(Mutation(position, mutation.region, alternate, muttype, mutation.sample, mutation.cancertype))
                df.append(l)
        # for a in df:
        #     print(a)

        # Start analysis
        logger.debug('Start analyzing simulations')
        for simulated_mutations in zip(*df):
            cutoff_clusters, element_score = self.analysis(element, simulated_mutations)
            sim_scores_chunk.append(element_score)
            for interval in cutoff_clusters:
                clusters = interval.data.copy()
                for cluster, values in clusters.items():
                    sim_cluster_chunk.append(values['score'])

        return element, sim_scores_chunk, sim_cluster_chunk

    @staticmethod
    def empirical_pvalue(observed, simulations):
        """
        Calculate empirical p-value using pseudocount (1)
        :param observed: int, observed score
        :param simulations: list, simulated score
        :return: float, p-value
        """
        return (len([x for x in simulations if x >= observed]) + 1) / len(simulations)

    def length(self, element):
        """Calculate length of an element (sum of input regions). If clustering in protein sequence, length in aa
        :param element: element to analyze
        :return length: int, length of an element (bp or aa)
        """
        length_ele = 0

        for interval in self.regions_d[element]:
            length_ele += (interval[1] - interval[0])
        if self.protein:
            return length_ele//3
        else:
            return length_ele

    def post_process(self, items):
        """
        Post processing of results
        :param item: tuple, elements results
        :return:
        """
        pseudo_pvalue = 1.1102230246251566e-19
        results = []
        for item in items:
            element, obs_clusters, obs_score, sim_clusters_scores, sim_element_scores = item
            mut_in_clusters = 0

            # Get GCG boolean
            cgc = element.split('//')[0] in self.cgc_genes

            # Calculate length
            element_length = self.length(element)

            # If analyzed element and element has clusters
            if type(obs_clusters) != float:
                if type(sim_element_scores) != float:

                    n_clusters_sim = len(sim_clusters_scores)

                    if sum(sim_element_scores) == 0:
                        logger.warning('No simulated clusters in {}. '
                                       'Observed cluster p-values calculated with pseudocount'.format(element))

                        # Calculate p-value for element
                        n_clusters = 0
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_element_scores)
                        obs_pvalue = obj.calculate(obs_score)

                        # Add pseudocount p-value for clusters
                        for interval in obs_clusters:
                            for cluster, values in interval.data.items():
                                interval.data[cluster]['p'] = pseudo_pvalue
                                n_clusters += 1
                                mut_in_clusters += len(values['mutations'])

                        empirical_pvalue = self.empirical_pvalue(obs_score, sim_element_scores)
                        analytical_pvalue = top_cluster_pvalue = obs_pvalue

                    else:
                        # Element score empirical p-value
                        empirical_pvalue = self.empirical_pvalue(obs_score, sim_element_scores)

                        # Element score analytical p-value
                        sim_scores_array_1000 = np.random.choice(sim_element_scores, size=1000, replace=False)
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_scores_array_1000)
                        analytical_pvalue = obj.calculate(obs_score)

                        # Cluster analytical p-values
                        obs_clusters_score = []

                        if n_clusters_sim < 3:  #TODO Check
                            for interval in obs_clusters:
                                for cluster, values in interval.data.items():
                                    interval.data[cluster]['p'] = pseudo_pvalue
                                    obs_clusters_score.append((values['score'], pseudo_pvalue))
                                    mut_in_clusters += len(values['mutations'])

                        else:
                            if n_clusters_sim > 1000:
                                # Random choice 1000 simulated cluster scores
                                sim_clusters_scores = np.random.choice(sim_clusters_scores, size=1000, replace=False)

                            # Fit distribution
                            obj = ap.AnalyticalPvalue()
                            obj.calculate_bandwidth(sim_clusters_scores)
                            for interval in obs_clusters:
                                for cluster, values in interval.data.items():
                                    cluster_p_value = obj.calculate(values['score'])
                                    interval.data[cluster]['p'] = cluster_p_value
                                    obs_clusters_score.append((values['score'], cluster_p_value))
                                    mut_in_clusters += len(values['mutations'])

                        n_clusters = len(obs_clusters_score)

                        # Element top cluster analytical p-value
                        top_cluster_pvalue = sorted(obs_clusters_score, key=lambda k: (k[0], -k[1]), reverse=True)[0][1]
                        logger.debug('P-values calculated')

                else:
                    n_clusters = obs_score = mut_in_clusters = 0
                    n_clusters_sim = float('nan')
                    empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')
            else:
                n_clusters = n_clusters_sim = obs_score = mut_in_clusters = float('nan')
                empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')

            # for i in obs_clusters:
            #     print(i[0], i[1])
            #     for c, v in i.data.items():
            #         print('\t', 'cluster:', c, 'data:', v['left_m'], v['max'], v['right_m'], len(v['mutations']), 'score:', v['score'])
            #         # print('\t', 'cluster:', c, 'data:', i[1]-v['left_m'][1], i[1]-v['max'][1], i[1]-v['right_m'][1], len(v['mutations']), 'score:', v['score'])

            results.append((
                element,
                (self.chromosomes_d[element], self.strands_d[element], element_length, len(self.mutations_d[element]),
                 mut_in_clusters, n_clusters, n_clusters_sim, obs_score,
                 empirical_pvalue, analytical_pvalue, top_cluster_pvalue, cgc),
                (obs_clusters, self.chromosomes_d[element], self.strands_d[element], element_length, cgc)
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
                    probs_tree, conseq_tree, skip = self.mut_probabilities(element)
                    if not skip:
                        element_size = ((len(self.mutations_d[element]) * self.n_simulations) // pf_num_simulations) + 1
                        chunk_size = (self.n_simulations // element_size) + 1
                        simulations += [(element, probs_tree, conseq_tree, n_sim) for n_sim in partitions_list(
                                        self.n_simulations,
                                        chunk_size
                        )]
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
