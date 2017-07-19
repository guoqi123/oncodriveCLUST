
# Oncodriveclustl functions

import gzip
from collections import defaultdict

import numpy as np
from intervaltree import IntervalTree
from scipy.signal import argrelmax
import math as m

### Parse regions
def regions(input_regions):
    """
    Parse input regions
    :param input_regions: path. File name? chr \t start \t end \t strand  \t ensembl id \t ensembl id \t symbol
    :return: dictionary of dictionary of intervaltrees (trees) containing intervals of genomic regions by chromosome.
    Dictionary containing a list of tuples with the coding sequences of each gene.
    """
    trees = defaultdict(IntervalTree)
    regions_d = defaultdict(list)

    with gzip.open(input_regions, 'rb') as fd:  # rb for binary <-- gz
        for line in fd:
            line = line.decode()  # binary to readable
            chromosome, start, end, strand, _, _, symbol = line.strip().split('\t')
            trees[chromosome][int(start): int(end) + 1] = symbol  # int, +1 end
            regions_d[symbol].append((int(start), int(end) + 1))

    return trees, regions_d


### Read mutations file and intersect with regions

def read_mutations(input_mutations, trees):
    """
    Parse substitution mutations
    :param input_mutations: path. File name? chr \t position \t ref \t alt \t sample \t mut_type \t tumor \t signature \ transcript \ symbol
    :param trees: dictionary of dictionary of intervaltrees containing intervals of genomic regions by chromosome.
    :return: dictionary containing a list of mutations by gene
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

### Smoothing function

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

### Smooth region

def smoothing(symbol, regions, mutations, window):
    """

    :param symbol: gene identifier
    :param mutations: list of mutations (genomic position) of the analyzed gene
    :param window: window (integer) of smoothing function
    :return:
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

    for element in regions[symbol]:  # iterate through tuples of coordinates
        for position in range(element[0], element[1] + 1):
            genomic = np.append(genomic, [position])

    binary = np.zeros(len(genomic) + window - 1)  # binary cds with borders added
    mutations_a = np.zeros(len(genomic) + window - 1)  # mutations cds with borders added

    indeces = np.searchsorted(np.array(genomic), mutations)


    indeces = np.searchsorted(np.array(genomic), mutations)

    try:
        for index in indeces:
            binary[index: index + window] += smooth
            mutations_a[index + window // 2] += 1
    except ValueError as e:
        print("Error, gene %s in pseudoautosomal region" %(symbol))

    # Remove the extra bases at the beginning and at the end
    binary = binary[window // 2: - window // 2 + 1]
    mutations_a = mutations_a[window // 2: - window // 2 + 1]

    region_info['symbol'] = symbol
    region_info['genomic'] = genomic
    region_info['binary'] = binary
    region_info['mutations'] = mutations_a


    return region_info

### Clusters

"""
- Find local maximum and minimum after smoothing for whole cds, according to smoothing score
- Define raw clusters: min(left), max, min(right) ----> Add filter? (ex. clusters wit max > smoothing score of 1 mutation)
- Starting from first maximum in the sequence, search for other relative maximum close to its min(right) border. 
- If there is one, check which of both maximum has a higher smoothing score and merge clusters by updating borders of the highest maximum. The lowest maximum (cluster) is removed from the dictionary of clusters. Iterate through all clusters until no updates are observed. 
- Score clusters. Two scores implemented: using position of smoothing maximum or mutation maximum within the cluster. 
"""

def clustering(regions, mutations, window):
    """
    For a given region, first find local maximum and minimum. Second, define a cluster per maximum. Third, merge
    clusters iteratively. Finally, score them.
    Scoring function: sum of scores of each position within the cluster. The score by position (i) is defined as:
    fraction of mutations (i) / sqrt(2) ** distance(i - position max n mutations or max smoothing)
    Fraction of mutations by position is the percentage of the number of mutations observed in the position out
    of the total number of mutations observed in the analyzed region across samples.
    
    :param regions: dictionary containing keys:
        symbol: gene identifier
        genomic: np.array containing all postions in the region analyzed
        binary: np.array of length == genomic, contains the sum of smoothing scores of all mutations in the region
        mutations: np.array of length == genomic, contains the number of mutations by position
    :param mutations: list of mutations (genomic position) of the analyzed region
    :param window: window (integer) of clustering
    :return:
        clusters: dictionary of dictionaries. Each cluster, named by position, is a dictionary of keys 'min_l', 'max', 
        'min_r'.    
    """
    
    indexes = []
    maxs = []
    clusters = defaultdict(dict)
    root = m.sqrt(2)

    ####  Find maximum and minumum 
       
    # Iterate through the binary array
    for i in range(len(regions['binary'])): # don't add 1, using indexes ---> range(0,len) == {0, len-1}
        
        # for max and min in range(1,len-1)
        if i !=0 and i != len(regions['binary'])-1: 
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
        
    maxs_set = set(maxs)
   
    """
    Generate raw clusters: 
    Each max defines a cluster
    Cluster: min_left(if there is one), max, min:right(if there is one)
    >>>>>Add filter? (ex. clusters with max > smoothing score of 1 mutation)
    """
  
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
                                
    #### Merge the clusters
    iterate = 1 
    while iterate != 0: # Iterate until no clusters updates occur
        stop = 0

        for x in range(len(clusters.keys())): 

            # When x is a key in clusters and min_r exists (clusters without min_r don't do merging): 
            if x in clusters.keys() and clusters[x]['min_r'][0] != []: 

                # Define the interval of search
                search_r = set(range(clusters[x]['max'][0] + 1, clusters[x]['min_r'][0] + window + 1))

                # When the intersection betweem the positions of the search and the maximums is not empty
                if search_r.intersection(maxs_set) != set():
                    intersect_max = maxs.index((list(search_r.intersection(maxs_set))[0])) # Analyze only the closest max
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


    #### Score clusters  
    for cluster, values in clusters.items():
        score = []

        if values['min_l'] != [] and values['min_r'] != []: 
            for position in range(values['min_l'][0], values['min_r'][0] + 1): # include mins
                fraction_mut = (regions['mutations'][position] / len(mutations)) * 100
                distance = abs(position - values['max'][0])
                score.append(fraction_mut / m.pow(root, distance))

        elif values['min_l'] == []:
            for position in range(values['max'][0] + 1, values['min_r'][0] + 1): 
                fraction_mut = (regions['mutations'][position] / len(mutations)) * 100
                distance = abs(position - values['max'][0])
                score.append(fraction_mut / m.pow(root, distance))

        elif values['min_r'] == []:
            for position in range(values['min_l'][0], values['max'][0]): 
                fraction_mut = (regions['mutations'][position] / len(mutations)) * 100
                distance = abs(position - values['max'][0])
                score.append(fraction_mut / m.pow(root, distance))

        # Update score
        clusters[cluster]['score'] = sum(score)


    # ADD SCORING == MUTATIONS?
   
    return clusters





### Individual functions

def find_locals1(binary):
    """
    Find local maximum and minimum
    :param binary: np.array of length == genomic, contains the sum of smoothing scores of all mutations in the region
    :return:
        indexes: list of lists [0 == min or 1 == max, position in cds, smoothing score] containing max and min
        plot_max: np.array for plot, 1 where a max is found, else nan
        plot_min: np.array for plot, 1 where a min is found, else nan
    """

    indexes = []
    plot_max = np.zeros(len(binary))
    plot_min = np.zeros(len(binary))

    # Iterate through the binary array
    for index in range(0, len(binary) - 1):

        # Max
        # when max score is equal in contiguous positions, assumes the first as max
        if binary[index] > binary[index - 1] and binary[index] >= binary[index + 1]:
            indexes.append([1, index, binary[index]])
            plot_max[index] += 1

            # Min
        elif binary[index] <= binary[index - 1] and binary[index] <= binary[index + 1]:
            # get a 'false' minimum where minimum score is equal in contiguous positions
            if binary[index - 1] or binary[index + 1] != 0:
                indexes.append([0, index, binary[index]])
                plot_min[index] += 1

    # Change 0 to nan for matplotlib
    plot_max[plot_max == 0.0] = np.nan
    plot_min[plot_min == 0.0] = np.nan

    return indexes, plot_max, plot_min

def cluster1(indexes, window, mutations, scoring):
    """
    Find clusters and score them.
    Scoring function: sum of scores of each position within the cluster. The score by position (i) is defined as:
    fraction of mutations (i) / sqrt(2) ** distance(i - position max n mutations or max smoothing)
    Fraction of mutations by position is the percentage of the number of mutations observed in the position out
    of the total number of mutations observed in the analyzed region across samples.

    :param indexes: list of lists [0 == min or 1 == max, position in cds, smoothing score]
    :param window: integer. Expands the right minimum (right border) of a cluster to search for nearer clusters
    :param scoring: choose one scoring function, "mutation" or "smoothing".
        mutation: the scoring function calculates the distance from the evaluated nucleotide to the position of maximum
        number of mutations within the cluster.
        smoothing: the scoring function calculates the distance from the evaluated nucleotide to the maximum peak of the
        scoring function within the cluster.

    return:
        clusters_d: dictionary of dictionaries. Contains clusters named by descending order of positions in the
        studied region. Clusters contain: max, min_left, min_right. Each of them has a list with the position in the
        region and the smoothing score.
    """

    total_maxs = []
    clusters_d = defaultdict(dict)

    # Get a list of maximums sorted by position
    for element in indexes:
        if element[0] == 1:
            total_maxs.append(element[1])

    total_maxs_set = set(total_maxs)

    # Generate raw clusters: cluster: min_left, max, min:right.
    # Name clusters by position, from 0 to len(maxs)-1 (integer >= 0 )
    # Add filter? (ex. clusters wit max > smoothing score of 1 mutation)
    j = 0
    for i in range(0, len(indexes) - 2):
        if indexes[i][0] == 0 and indexes[i + 1][0] == 1:
            clusters_d[j]['min_l'] = [indexes[i][1], indexes[i][2]]
            clusters_d[j]['max'] = [indexes[i + 1][1], indexes[i + 1][2]]
            clusters_d[j]['min_r'] = [indexes[i + 2][1], indexes[i + 2][2]]
            j += 1

    # Merge the clusters
    iterations = len(clusters_d.keys())
    iterate = 1

    # Iterate until no clusters_d updates occur
    while iterate != 0:
        stop = set()

        # Iterate through initial number of clusters
        # IMPROVE?
        for i in range(iterations):

            # When i is a key in clusters_d == when i is a cluster
            if i in clusters_d:
                # Define the interval of search
                search_r = set(range(clusters_d[i]['max'][0] + 1, clusters_d[i]['min_r'][0] + window + 1))

                # When the intersection betweem the positions of the search and the maximums is not empty
                if search_r.intersection(total_maxs_set) != set():
                    # Analyze only the closest max
                    intersect_max = total_maxs.index(
                        list(search_r.intersection(total_maxs_set))[0])  # sort / assume its sorted?
                    stop.add(1)

                    # When the testing max is greater or equal than the intersected max,
                    # expand the right border and delete intersected cluster from clusters_d
                    if clusters_d[i]['max'][1] >= clusters_d[intersect_max]['max'][1]:
                        clusters_d[i]['min_r'] = clusters_d[intersect_max]['min_r']
                        del clusters_d[intersect_max]
                        total_maxs_set.remove(total_maxs[intersect_max])


                    # When the testing max is smaller than the intersected max,
                    # expand the left border of intersected max and remove the testing cluster (i)
                    elif clusters_d[i]['max'][1] < clusters_d[intersect_max]['max'][1]:
                        clusters_d[intersect_max]['min_l'] = clusters_d[i]['min_l']
                        del clusters_d[i]
                        total_maxs_set.remove(total_maxs[i])

        if stop == set():
            iterate = 0

    # Score clusters
    root = m.sqrt(2)

    if scoring == "smoothing": 
        for cluster, values in clusters_d.items():
            score = 0

            for position in range(values['min_l'][0], values['min_r'][0] + 1):  # includes mins
                fraction_mut = (mutations[position] / len(mutations)) * 100
                distance = abs(position - values['max'][0])
                score += (fraction_mut / m.pow(root, distance))

            # Update score
            clusters_d[cluster]['score'] = score


    elif scoring == "mutation":
        indices = [i for i, x in enumerate(mutations) if x > 0]

        for cluster, values in clusters_d.items():

            # For each cluster, find position with max number of mutations
            interval_mut = range(values['min_l'][0] + 1, values['min_r'][0])  # don't includes mins
            cluster_mutations_in = list(filter(lambda x: x in interval_mut, indices))

            cluster_mutations_val = []
            
            for element in cluster_mutations_in:
                cluster_mutations_val.append((element, mutations[element]))
                max_mutation_index = max(cluster_mutations_val, key=lambda x: x[1])

            clusters_d[cluster]['max_mut'] = max_mutation_index[
                0]  # What if cluster with positions with same number of mutations???

            # Score cluster
            score = 0
            
            for position in range(values['min_l'][0], values['min_r'][0] + 1):  # includes mins
                fraction_mut = (mutations[position] / len(mutations)) * 100
                distance = abs(position - values['max_mut'])
                score += fraction_mut / m.pow(root, distance)
            # Update score
            clusters_d[cluster]['score'] = score

    return clusters_d