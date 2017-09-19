# Oncodriveclustl run
import logging
import os

import click
import daiquiri
import pandas as pd

from oncodriveclustl.utils import signature as sign
from oncodriveclustl.utils import parsing as pars
from oncodriveclustl.utils import run as exp


# Global variables
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def write_element_results(results, output_file):
    """Save results to the output file
    :param results: dict, dictionary of results, keys are element's symbols
    :param output_file: path, path of the output file
    :return: None
    """

    header = ['SYM', 'LEN', 'N_MUT',
              'CLU', 'MEAN_SIM_CLU', 'MEDIAN_SIM_CLU', 'SD_SIM_CLU',
              'SCORE_OBS', 'MEAN_SIM_SCORE', 'MEDIAN_SIM_SCORE', 'SD_SIM_SCORE',
              'E_PVAL', 'A_PVAL', 'CGC']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for gene_name, values in results.items():
            length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, \
            mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc = values
            fd.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\t{}\t{}\n'.format(
               gene_name, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))
    # Sort
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['E_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)


def write_cluster_results(results, output_file, sorter):
    """Save results to the output file
    :param results: dict, dictionary of results, keys are element's symbols
    :param output_file: path, path of the output file
    :param sorter: list, element symbols ranked by elements p-value
    :return: None
    """

    sorterindex = dict(zip(sorter, range(len(sorter))))

    header = ['RANK','SYMBOL', 'CLUSTER', 'N', 'MIN_L', 'MAX', 'MIN_R', 'WIDTH', 'MUT', 'SCORE', 'CGC']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element, values in results.items():
            clustersinfo, cgc = values
            rank = sorterindex[element]
            for c, v in clustersinfo.items():
                if v['min_l'] and v['min_r']:
                    fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rank, element, v['mode'], c, v['min_l'][0], v['max'][0], v['min_r'][0], abs(v['min_r'][0]-v['min_l'][0]), v['n_mutations'], v['score'], cgc))
                elif not v['min_l']:
                    fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rank, element, v['mode'], c, '-', v['max'][0], v['min_r'][0], abs(v['min_r'][0]-v['max'][0]), v['n_mutations'], v['score'], cgc))
                elif not v['min_r']:
                    fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        rank, element, v['mode'], c, v['min_l'][0], v['max'][0], '-', abs(v['max'][0]-v['min_l'][0]), v['n_mutations'], v['score'], cgc))
    # Sort
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['RANK', 'SCORE'], ascending=[True, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)


@click.command()
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations')
@click.option('-o', '--output-file', default=None, required=True,
              help='File to save with clusters')
@click.option('-oc', '--output-file-clusters', default=None, required=True,
              help='File to save with clusters')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='File with the genomic regions to analyze')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbol of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element to analyze')
@click.option('--element-mutations', type=click.INT, default=2,
              help='Cutoff of element mutations. Default is 2')
@click.option('--cluster-mutations', type=click.INT, default=2,
              help='Cutoff of cluster mutations. Default is 2')
@click.option('-sw', '--smooth-window', type=click.INT, default=50,
              help='smoothing window. Default is 50')
@click.option('-cw', '--cluster-window', type=click.INT, default=50,
              help='cluster window. Default is 50')
@click.option('--cluster-score', default='nobias', help='Cluster score formula',
              type=click.Choice(['nobias', 'fmutations']))
@click.option('--element-score', default='mean', help='Gene score formula',
              type=click.Choice(['sum', 'mean']))
@click.option('-n', '--n-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('--seed', type=click.INT, default=None,
              help='seed to use in the simulations')
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file,
         output_file,
         output_file_clusters,
         regions_file,
         elements_file,
         elements,
         element_mutations,
         cluster_mutations,
         smooth_window,
         cluster_window,
         cluster_score,
         element_score,
         n_simulations,
         cores,
         seed,
         log_level):
    """Oncodriveclustl is a program that looks for mutational hotspots
    :param input_file: input file
    :param output_file: output file
    :param output_file_clusters: output file with clusters information
    :param regions_file: path input genomic regions, tab file
    :param elements_file: file containing one element per row
    :param elements: element symbol or file containing elements
    :param element_mutations: int, cutoff of element mutations
    :param cluster_mutations: int, cutoff of cluster mutations
    :param smooth_window: int, smoothing window
    :param cluster_window: int, clustering window
    :param cluster_score: cluster score method
    :param element_score: gene score method
    :param n_simulations: int, number of simulations
    :param cores: int, number of CPUs to use
    :param seed: int, seed
    :param log_level: verbosity of the logger
    :return: None
    """

    daiquiri.setup(level=LOGS[log_level])
    logger = daiquiri.getLogger()
    logger.debug(' '.join([input_file,
                           output_file,
                           output_file_clusters,
                           regions_file,
                           str(element_mutations),
                           str(cluster_mutations),
                           str(smooth_window),
                           str(cluster_window),
                           cluster_score,
                           element_score,
                           str(n_simulations),
                           str(cores)]))

    # Create a list of elements to analyze
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

    # Parse regions and dataset mutations
    logger.info('Parsing input regions and input mutations...')
    regions_d, chromosomes_d, mutations_d = pars.parse(regions_file, elements, input_file)

    # Compute dataset trinucleotide signatures
    logger.info('Computing signatures...')
    obj = sign.Signature(start_at_0=True)
    obj.calculate(input_file)

    # TODO remove this hardcoded file?
    obj.save(os.path.join(os.path.dirname(__file__), 'cache/signature.pickle'))

    # Initialize Experiment class variables and run
    elements_results, clusters_results = exp.Experiment(
                                regions_d, chromosomes_d, mutations_d,
                                element_mutations, cluster_mutations,
                                smooth_window, cluster_window,
                                cluster_score, element_score,
                                n_simulations, cores, seed
                                ).run()

    # Write results
    write_element_results(results=elements_results, output_file=output_file)
    logger.info('Gene results calculated')
    df = pd.read_csv(output_file, sep='\t')
    sorted_list_elements = df['SYM'].tolist()
    write_cluster_results(results=clusters_results, output_file=output_file_clusters, sorter=sorted_list_elements)
    logger.info('Clusters results calculated')
    logger.info('Finished')


if __name__ == '__main__':
    main()
