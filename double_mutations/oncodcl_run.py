# Oncodriveclustl run
import os
import sys
import logging
from collections import defaultdict

import click
import daiquiri
from concurrent.futures import ProcessPoolExecutor as Pool
from tqdm import tqdm
import pandas as pd

import signature as sign
import oncodcl_funct as odf
import oncodcl_funct2 as odf2



# Global variables
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def write_results(results, output_file):
    """Save results to the output file
    :param results: dict, dictionary of results, keys are gene's names
    :param output_file: path, path of the output file
    :return: None
    """
    header = ['SYMBOL', 'SCORE_OBS', 'EP-VALUE', 'AP-VALUE', 'CGC']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for gene_name, values in results.items():
            score_obs, epval, apval, cgc = values
            fd.write('{}\t{}\t{}\t{}\t{}\n'.format(
                gene_name, score_obs, epval, apval, cgc))

    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['EP-VALUE', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='')


def write_sims(results, output_file):
    """Save results to the output file
    :param results: dict, dictionary of results, keys are gene's names
    :param output_file: path, path of the output file
    :return: None
    """

    with open(output_file, 'w') as fd:
        for gene_name, values in results.items():
            score_obs, sims = values
            for score in sims:
                fd.write('{}\n'.format(score))

    #df = pd.read_csv(output_file, sep='\t', header=0)
    #df.sort_values(by=['P-VALUE', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    #df.to_csv(path_or_buf=output_file, sep='\t', na_rep='')



def run(regions_d, chromosomes_d, mutations_d,output_file, n_simulations, smooth_window, cluster_window, cores):
    """Main function
    :param regions_d: dictionary of dictionaries containing genomic regions
    :param chromosomes_d: dict, keys are genes, values are chromosomes
    :param mutations_d: dictionary of lists, 'gene' : [mutations]
    :param output_file: output file
    :param n_simulations: int, number of simulations
    :param smooth_window: int, smoothing window
    :param cluster_window: int, clustering window
    :param cores: int, number of CPUs to use
    :return: None
    """
    results_d = defaultdict()

    # Analyze only genes with >= 2 mutations
    genes = [(g, regions_d[g], chromosomes_d[g], mutations_d[g], n_simulations, smooth_window, cluster_window, cores)
             for g in regions_d.keys() if len(mutations_d[g]) >= 2]

    # Read CGC
    CGC_genes = set([line.split('\t')[0] for line in
                     open('../inputs/CGC/CGCMay17_cancer_types_TCGA.tsv', 'r')])

    logger.info('Calculating results {} genes...'.format(len(genes)))
    with tqdm(total=len(genes)) as pbar:
        for gene, obs, epval, apval in map(odf2.run_region, genes):
            pbar.update(1)
            CGC = gene in CGC_genes
            results_d[gene] = (obs, epval, apval, CGC)

    logger.info('Results calculated')
    write_results(results_d, output_file)
    logger.info('Finished')



@click.command()
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations')
@click.option('-o', '--output-file', default=None, required=True,
              help='File to save with clusters')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='File with the genomic regions to analyze')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbol of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element to analyze')
@click.option('-n', '--n-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-sw', '--smooth-window', type=click.INT, default=50,
              help='smoothing window. Default is 50')
@click.option('-cw', '--cluster-window', type=click.INT, default=50,
              help='cluster window. Default is 50')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file, output_file, regions_file, elements_file, elements,
         n_simulations, smooth_window, cluster_window, cores, log_level):
    """Oncodriveclustl is a program that looks for mutational hotspots"""
    global logger

    daiquiri.setup(level=LOGS[log_level])
    logger = daiquiri.getLogger()

    logger.debug(' '.join([input_file,
                           output_file,
                           regions_file,
                           str(n_simulations),
                           str(smooth_window),
                           str(cluster_window),
                           str(cores)]))

    # Create a list of elements to analyze
    print(elements)
    if elements is not None:
        elements = set(elements)
    if elements_file is not None:
        elements |= set([line.strip().split()[0] for line in open(elements_file, 'r')])
    if elements is None and elements_file is None:
        elements = set([])

    if elements:
        logger.info(
            'The following {} element{} will be analyzed'.format(len(elements), 's' if len(elements) > 1 else '')
        )
        logger.info(', '.join(elements))

    logger.info('Computing region trees...')
    trees, regions_d, chromosomes_d = odf.regions(regions_file, elements)

    logger.info('Parsing mutations...')
    mutations_d = odf.read_mutations(input_file, trees)

    logger.info('Computing signatures...')
    obj = sign.Signature(start_at_0=True)
    obj.calculate(input_file)
    obj.save('signature.pickle')

    run(regions_d, chromosomes_d, mutations_d, output_file, n_simulations, smooth_window, cluster_window, cores)

if __name__ == '__main__':
    main()