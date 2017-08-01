

# Oncodriveclustl functions
import sys
import argparse
import re
from concurrent.futures import ProcessPoolExecutor as Pool
from collections import defaultdict

import pandas as pd
from tqdm import tqdm
import oncodcl_funct as odf

parser = argparse.ArgumentParser(description="""OncodriveCLUSTL testing""")

parser.add_argument("-GM", "--genemut",
                    dest = "GM",
                    action = "store",
                    default = 2,
                    type = int,
                    help = "Input number of mutations required in a gene")

parser.add_argument("-CM", "--clustermut",
                    dest = "CM",
                    action = "store",
                    default = 2,
                    type = int,
                    help = "Input number of mutations required in a cluster")

parser.add_argument("-SW", "--smoothw",
                    dest = "SW",
                    action = "store",
                    default = None,
                    required = True,
                    type = int,
                    help = "Smoothing window")

parser.add_argument("-CW", "--clusterw",
                    dest = "CW",
                    action = "store",
                    default = None,
                    required = True,
                    type = int,
                    help = "Cluster window")

parser.add_argument("-CS", "--clusterscore",
                    dest = "clusterscore",
                    default = "nobias",
                    action = "store",
                    help = "Cluster score: 'nmutations', 'genelength', 'nobias' ")

parser.add_argument("-GS", "--genescore",
                    dest = "genescore",
                    action = "store",
                    default = "mean",
                    help = "Gene score: 'sum', 'mean', 'wmean'")

parser.add_argument('-o', '--output',
                    dest = "outfile",
                    action = "store",
                    default = None,
                    required=True,
                    help = "Output file prefix without extension")

args = parser.parse_args()

sys.stderr.write("OncodriveCLUSTL testing\n")
sys.stderr.write("\tGene mutations threshold (GM):\t%s\n" %args.GM)
sys.stderr.write("\tCluster mutations threshold (CM): %s\n" %args.CM)
sys.stderr.write("\tSmoothing window (SW):\t%s\n" %args.SW)
sys.stderr.write("\tClustering window (CW):\t%s\n" %args.CW)
sys.stderr.write("\tCluster score mode (CS): %s\n" %args.clusterscore)
sys.stderr.write("\tGene score mode (GS):\t%s\n" %args.genescore)
outfile_prefix = re.sub('\..*', '', str(args.outfile)) # Remove any extension added
outfile_file = outfile_prefix+"_"+str(args.GM)+"_"+str(args.CM)+"_"+str(args.SW)+"_"+str(args.CW)+"_"+str(args.clusterscore)+"_"+str(args.genescore)
sys.stderr.write("\tOutput file:\t\t%s\n" %outfile_file)

# Parse regions
input_regions = '/home/carnedo/projects/oncodriveclustl/inputs/regions/02_cds.regions.gz'
trees, regions_d, strands_d = odf.regions(input_regions)
sys.stderr.write('Regions parsed\n')

# Read mutations, intersect with regions
input_mutations = '/home/carnedo/projects/oncodriveclustl/inputs/mutations/pancanatlas/SKCM.txt'
mutations_d = odf.read_mutations(input_mutations, trees)
sys.stderr.write('Mutations read\n')


# Smooth
total_regions = defaultdict(dict)
for gene in regions_d.keys():
    # Analyze only genes with mutations >= threshold
    if len(mutations_d[gene]) >= args.GM:
        region_info = odf.smoothing(symbol=gene, regions=regions_d, mutations=mutations_d[gene], window=args.SW)
        # dictionary with regions of all genes
        total_regions[region_info['symbol']] = region_info
sys.stderr.write('Smoothing finished. Total genes for clustering: %s\n' %len(total_regions.keys()))

# Clusters
total_clusters = defaultdict(dict)
for gene in total_regions.keys():
    clusters = odf.clustering(regions=total_regions[gene], mutations=mutations_d[gene], window=args.CW, cutoff=args.CM, mode=args.clusterscore)
    if len(clusters.keys()) != 0:
        total_clusters[gene] = clusters
sys.stderr.write('Clustering finished. Total genes analyzed: %s\n' %len(total_clusters.keys()))

# Score genes
genes_results = defaultdict()
for gene, value in total_clusters.items():
    print(gene, value)
    result = odf.score_genes(clusters=value, mode=args.genescore)
    genes_results[gene] = result

# Output
CGC = set([line.split('\t')[0] for line in open('/home/carnedo/projects/oncodriveclustl/inputs/CGC/CGCMay17_cancer_types_TCGA.tsv', 'r')])
out = '/home/carnedo/projects/oncodriveclustl/outputs/20170731'

genes_table = out + outfile_file

with open(genes_table, 'w') as fd:
    fd.write('Gene\tLength\tClusters\tMutations\tScore\n')
    for gene, result in genes_results.items():
        fd.write("%s\t%s\t%s\t%s\t%s\n" % (gene, len(total_regions[gene]['genomic']), result['n_clusters'], len(mutations_d[gene]), result['score']))

df = pd.read_csv(genes_table, sep='\t', header=0)
df['CGC'] = df['Gene'].map(lambda x: x in CGC)
df.sort_values(by=['Score', 'CGC'], ascending=[False, False], inplace=True)
df.to_csv(path_or_buf=genes_table, sep='\t', na_rep='', index=False)


def run(gene):
    smooting()
    clustering()
    score()
    write_output()
    return gene, result


def main(cores=4, genes):
    """Main function
    :param cores: int, number of CPUs to use
    :param genes: lis, list of genes to analize
    :return: None
    """
    results = dict()
    with Pool(max_workers=cores) as executor, tqdm(total=len(genes)) as pbar:
        for geneid, result in executor.map(run, genes):
            pbar.update(1)
            results[geneid] = result
            # save(results, output_file)