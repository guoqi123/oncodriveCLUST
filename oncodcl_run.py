# Oncodriveclustl run
import sys
from collections import defaultdict

from concurrent.futures import ProcessPoolExecutor as Pool
from tqdm import tqdm
import pandas as pd

import oncodcl_funct as odf
import oncodcl_funct2 as odf2
import signature as sign

# Global variables
OUTPUT_FILE = '/home/carnedo/projects/oncodriveclustl/outputs/SKCM.txt'


def write_results(results, output_file):
    """Save results to the output file
    :param results: dict, dictionary of results, keys are gene's names
    :param output_file: path, path of the output file
    :return: None
    """
    header = ['SYMBOL', 'SCORE_OBS', 'SCORE_SIM', 'CGC']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for gene_name, values in results.items():
            score_obs, score_sim, cgc = values
            fd.write('{}\t{}\t{}\t{}\n'.format(
                gene_name, score_obs, score_sim, cgc))

    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['SCORE_OBS', 'CGC'], ascending=[False, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='')


def main(regions_d, chromosomes_d, mutations_d, cores=4):
    """Main function
    :param regions_d: dictionary of dictionaries containing genomic regions
    :param chromosomes_d: dict, keys are genes, values are chromosomes
    :param mutations_d: dictionary of lists, 'gene' : [mutations]
    :param cores: int, number of CPUs to use
    :return: None
    """
    results_d = defaultdict()

    # Analyze only genes with >= 2 mutations
    genes = [(g, regions_d[g], chromosomes_d[g], mutations_d[g]) for g in regions_d.keys() if len(mutations_d[g]) >= 2]
    CGC_genes = set([line.split('\t')[0] for line in
                     open('/home/carnedo/projects/oncodriveclustl/inputs/CGC/CGCMay17_cancer_types_TCGA.tsv', 'r')])

    with Pool(max_workers=4) as executor, tqdm(total=len(genes)) as pbar:
        for gene, scores in executor.map(odf2.run_region, genes):
            pbar.update(1)
            CGC = gene in CGC_genes
            results_d[gene] = (scores['obs'], scores['sim'], CGC)

    sys.stderr.write('Results calculated\n')
    write_results(results_d, OUTPUT_FILE)

# Parse regions
input_regions = '/home/carnedo/projects/oncodriveclustl/inputs/regions/02_cds.regions.gz'
trees, regions_d, chromosomes_d = odf.regions(input_regions)
sys.stderr.write('Regions parsed\n')

# Read mutations, intersect with regions
input_mutations = '/home/carnedo/projects/oncodriveclustl/inputs/mutations/pancanatlas/SKCM.txt'
mutations_d = odf.read_mutations(input_mutations, trees)
sys.stderr.write('Mutations read\n')

# Calculate signatures for cancer type dataset
obj = sign.Signature(start_at_0=True)
obj.calculate(input_mutations)
signatures = obj.signatures['probabilities']
sys.stderr.write('Signatures calculated\n')


if __name__ == '__main__':
    main(regions_d, chromosomes_d, mutations_d)