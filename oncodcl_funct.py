# Oncodriveclustl functions
import gzip
import math as m
from collections import defaultdict
import numpy as np
from intervaltree import IntervalTree


def regions(input_regions):
    """
    Parse input regions
    :param input_regions: path tab file chr \t start \t end \t strand  \t ensembl id \t ensembl id \t symbol
    :return:
        trees: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic regions by chromosome.
        regions_d: dictionary containing a list of tuples with the coding sequences of each gene.
    """
    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(list)

    with gzip.open(input_regions, 'rb') as fd:  # rb for binary <-- gz
        for line in fd:
            line = line.decode()  # binary to readable
            chromosome, start, end, strand, ensid, _, symbol = line.strip().split('\t')
            if 'ENSGR' not in ensid:
                trees[chromosome][int(start): int(end) + 1] = symbol  # int, +1 end
                regions_d[symbol].append((int(start), int(end) + 1))

    return trees, regions_d


# Read mutations file and intersect with regions
def read_mutations(input_mutations, trees):
    """
    Parse substitution mutations
    :param input_mutations: path tab file
    chr \t position \t ref \t alt \t sample \t mut_type \t tumor \t signature \ transcript \ symbol
    :param trees: dictionary of dictionary of intervaltrees containing intervals of genomic regions by chromosome.
    :return:
        mutations_d: dictionary, key = gene, value = list of mutations per gene
    """
    mutations_d = defaultdict(list)

    with open(input_mutations, 'r') as fd:
        next(fd)  # don't read header
        for line in fd:
            chromosome, position, ref, alt, sample, mut_type, tumor = line.strip().split('\t')[:7]
            if trees[chromosome][int(position)] != set() and mut_type == 'subs':  # intersect only subs
                results = trees[chromosome][int(position)]
                for res in results:
                    mutations_d[res.data].append(int(position))  # dict: key = gene, value = list of mutations

    return mutations_d


# Smooth
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


def smoothing(symbol, regions, mutations, window):
    """
    Given mutations in a region, generate a smoothing score for position
    :param symbol: gene symbol
    :param regions: regions analysed (ex. cds in a gene)
    :param mutations: gene mutations
    :param window: smoothing window
    :return:
        region_info: dictionary containing information for a gene. Keys:
            symbol: gene identifier
            genomic: np.array containing all postions in the region analyzed
            binary: np.array of length == genomic, contains the sum of smoothing scores of all mutations in the region
            mutations: np.array of length == genomic, contains the number of mutations by position
    """
    # Define smoothing window
    window += 1 - window % 2
    smooth = tukey(window)

    # Generate genomic, binary and mutation arrays and store them in a dictionary
    region_info = defaultdict()
    genomic = np.array([])  # genomic coordinates

    for element in regions:  # iterate through tuples of coordinates
        for position in range(element[0], element[1] + 1):
            genomic = np.append(genomic, [position])

    binary = np.zeros(len(genomic) + window - 1)  # binary cds with borders added
    mutations_a = np.zeros(len(genomic) + window - 1)  # mutations cds with borders added

    indeces = np.searchsorted(np.array(genomic), mutations)

    for index in indeces:
        binary[index: index + window] += smooth
        mutations_a[index + window // 2] += 1

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: - window // 2 + 1]
    mutations_a = mutations_a[window // 2: - window // 2 + 1]

    region_info['symbol'] = symbol
    region_info['genomic'] = genomic
    region_info['binary'] = binary
    region_info['mutations'] = mutations_a

    return region_info


# Clusters
def find_locals(regions):
    """
    For a given region, find local maximum and minimum
    :param regions:  dictionary containing information for a gene (symbol, genomic, binary, mutations)
    :return:
        idexes: list of tuples of 3 elements (min or max, cds region, smoothing score). Min == 0, max == 1
        maxs: list of positions of maximum in cds
    """
    indexes = []
    maxs = []
       
    # Iterate through the binary array
    for i in range(len(regions['binary'])):  # don't add 1, using indexes ---> range(0,len) == {0, len-1}
        
        # for max and min in range(1,len-1)
        if i != 0 and i != len(regions['binary'])-1:
            # Max
            # When max score is equal in contiguous positions, assumes the first as max
            if regions['binary'][i] > regions['binary'][i - 1] and regions['binary'][i] >= regions['binary'][i + 1]:
                indexes.append([1, i, regions['binary'][i]])
                # Get a list of maximums sorted by position
                maxs.append(i)
            # Min
            elif regions['binary'][i] <= regions['binary'][i - 1] and regions['binary'][i] <= regions['binary'][i + 1]:
                # get a 'false' minimum where minimum score is equal in contiguous positions
                if regions['binary'][i - 1] or regions['binary'][i + 1] != 0:
                    indexes.append([0, i, regions['binary'][i]])            
    
        # for max or min in first position
        elif i == 0:
            if regions['binary'][i] >= regions['binary'][i + 1] and regions['binary'][i] != 0: 
                indexes.append([1, i, regions['binary'][i]])
                maxs.append(i)
            elif regions['binary'][i] <= regions['binary'][i + 1] and regions['binary'][i + 1] != 0: 
                indexes.append([0, i, regions['binary'][i]])            
      
        # for max or min in last position
        elif i == len(regions['binary'])-1:
            if regions['binary'][i] > regions['binary'][i - 1]: 
                indexes.append([1, i, regions['binary'][i]])
                maxs.append(i)
            elif regions['binary'][i] < regions['binary'][i - 1]: 
                indexes.append([0, i, regions['binary'][i]])  
                
    return indexes, maxs


def raw_clusters(indexes):
    """
    Define a cluster per maximum found in a region
    :param indexes: list of tuples of 3 elements (min or max, cds region, smoothing score). Min == 0, max == 1
    :return:
        clusters: dictionary of dictionaries of clusters for a gene. Each cluster, named by position from 0 to n of
        clusters, is a dictionary of keys 'min_l', 'max', 'min_r'.
    """
    clusters = defaultdict(dict)
    
    # Iterate through all maxs in indexes
    j = 0 
    generator_maxs = (i for i in indexes if i[0] == 1)
        
    for maximum in generator_maxs:
        i = indexes.index(maximum)
        clusters[j]['max'] = [maximum[1], maximum[2]]     

        # if it's not the first nor the last cluster
        if i != 0 and i != len(indexes)-1: 
            clusters[j]['min_l'] = [indexes[i - 1][1], indexes[i - 1][2]]
            clusters[j]['min_r'] = [indexes[i + 1][1], indexes[i + 1][2]]

        # if it's the first cluster
        elif i == 0: 
            clusters[j]['min_l'] = []
            clusters[j]['min_r'] = [indexes[i + 1][1], indexes[i + 1][2]]

        # if it's the last cluster
        elif i == len(indexes)-1: 
            clusters[j]['min_l'] = [indexes[i - 1][1], indexes[i - 1][2]]
            clusters[j]['min_r'] = []
        
        j += 1

    return clusters


def merge_clusters(clusters, maxs, window):
    """
    Given a number of clusters in a region, iterate through them to merge them if their maximums are closer than a given
    length.
    :param clusters: dictionary of dictionary of clusters for a gene.
    :param maxs: list of maximum positions in a region
    :param window: int, clustering window. It is added to cluster right minimum
    :return:
        clusters: dictionary of dictionaries
    """
    maxs_set = set(maxs)

    iterate = 1 
    while iterate != 0: # Iterate until no clusters updates occur
        stop = 0
        for x in range(len(clusters.keys())):

            # When x is a key in clusters and min_r exists (clusters without min_r don't do merging): 
            if x in clusters.keys() and clusters[x]['min_r'][0] != []: 

                # Define the interval of search
                search_r = set(range(clusters[x]['max'][0] + 1, clusters[x]['min_r'][0] + window + 1))

                # When the intersection between the positions of the search and the maximums is not empty
                if search_r.intersection(maxs_set) != set():
                    # Analyze only the closest max
                    intersect_max = maxs.index(sorted(list(search_r.intersection(maxs_set)))[0])
                    stop = 1

                    # When the testing max is greater or equal than the intersected max,
                    # expand the right border and delete intersected cluster from clusters
                    if clusters[x]['max'][1] >= clusters[intersect_max]['max'][1]:
                        clusters[x]['min_r'] = clusters[intersect_max]['min_r']
                        del clusters[intersect_max]
                        maxs_set.remove(maxs[intersect_max])

                    # When the testing max is smaller than the intersected max,
                    # expand the left border of intersected max and remove the testing cluster (i)
                    elif clusters[x]['max'][1] < clusters[intersect_max]['max'][1]:
                        clusters[intersect_max]['min_l'] = clusters[x]['min_l']
                        del clusters[x]
                        maxs_set.remove(maxs[x])

        if stop == 0:
            iterate = 0

    return clusters


def clusters_mut(clusters, regions, mutations):
    """
    Calculates the number of mutations within a cluster
    :param clusters: dictionary of dictionary containing clusters
    :param regions: dictionary containing information for a gene (symbol, genomic, binary, mutations)
    :param mutations: list of mutations for a region
    :return:
        clusters: dictionary of dictionaries with new key 'n_mutations' added
    """
    for cluster, values in clusters.items():

        if values['min_l'] and values['min_r']:
            cluster_range = range(int(regions['genomic'][values['min_l'][0]]), int(regions['genomic'][values['min_r'][0]] + 1))
            cluster_muts = [i for i in mutations if i in cluster_range]
            clusters[cluster]['n_mutations'] = len(cluster_muts)

        elif not values['min_l']:
            cluster_range = range(int(regions['genomic'][values['max'][0]]), int(regions['genomic'][values['min_r'][0]] + 1))
            cluster_muts = [i for i in mutations if i in cluster_range]
            clusters[cluster]['n_mutations'] = len(cluster_muts)

        elif not values['min_r']:
            cluster_range = range(int(regions['genomic'][values['min_l'][0]]), int(regions['genomic'][values['min_r'][0]] + 1))
            cluster_muts = [i for i in mutations if i in cluster_range]
            clusters[cluster]['n_mutations'] = len(cluster_muts)

    return clusters


def score_clusters(clusters, regions):
    """
    Score clusters
    :param clusters: dictionary of dictionaries of clusters
    :param regions: dictionary containing information for a gene (symbol, genomic, binary, mutations)
    :return: dictionary of dictionaries, with key 'score' added
    """
    root = m.sqrt(2)

    for cluster, values in clusters.items():
        score = []

        if values['min_l'] and values['min_r']:
            for position in range(values['min_l'][0], values['min_r'][0] + 1):  # include mins
                fraction_mut = regions['mutations'][position]
                distance = abs(values['max'][0] - position)
                score.append(fraction_mut / m.pow(root, distance))

        elif not values['min_l']:
            for position in range(values['max'][0] + 1, values['min_r'][0] + 1):
                fraction_mut = regions['mutations'][position]
                distance = abs(values['max'][0] - position)
                score.append(fraction_mut / m.pow(root, distance))

        elif not values['min_r']:
            for position in range(values['min_l'][0], values['max'][0]):
                fraction_mut = regions['mutations'][position]
                distance = abs(values['max'][0] - position)
                score.append(fraction_mut / m.pow(root, distance))

        # Update score
        clusters[cluster]['score'] = sum(score)

    return clusters


def clustering(regions, mutations, window):
    """
    Given a gene, calculate and score its clusters
    :param regions: dictionary containing information for a gene (symbol, genomic, binary, mutations)
    :param mutations: list of mutations for a region
    :param window: int, clustering window for clusters merging
    :return: dictionary of dictionaries of clusters within a region
    """
    indexes, maxs_list = find_locals(regions=regions)
    r_clusters = raw_clusters(indexes)
    m_clusters = merge_clusters(clusters=r_clusters, maxs=maxs_list, window=window)
    f_clusters = clusters_mut(clusters=m_clusters, regions=regions, mutations=mutations)
    s_clusters = score_clusters(clusters=f_clusters, regions=regions)
    
    return s_clusters


# Score genes
def score_gene(clusters, cutoff):
    """
    Given the clusters of a region, calculate a global score for it
    :param clusters: dictionary of dictionaries
    :param cutoff: int, n cluster mutations cutoff
    :return: int, int; number of clusters and gene score
    """
    n_clusters = []
    score = []

    for cluster, value in clusters.items():
        if value['n_mutations'] >= cutoff:
            n_clusters.append(1)
            score.append(value['score'])

    gene_score = 0
    if n_clusters:
        gene_score = sum(score) / sum(n_clusters)

    return sum(n_clusters), gene_score
