# Oncodriveclustl run
import logging
import os

import click
import daiquiri
import pandas as pd

from utils import signature as sign
from utils import parsing as pars
from utils import run as exp


# Global variables
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def write_element_results(genome, results, dir, file):
    """Save results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param dir: str, output directory
    :param file: str, output file, if elements in elements file
    :return: None
    """

    header = ['SYM', 'LEN', 'N_MUT',
              'CLU', 'MEAN_SIM_CLU', 'MEDIAN_SIM_CLU', 'SD_SIM_CLU',
              'SCORE_OBS', 'MEAN_SIM_SCORE', 'MEDIAN_SIM_SCORE', 'SD_SIM_SCORE',
              'E_PVAL', 'A_PVAL', 'CGC']

    output_file = dir+'/elements_'+file+'.txt'
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for gene_name, values in results.items():
            length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, \
            mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc = values
            if genome != 'hg19':
                cgc = 'Non Available'
            if obs_clusters and type(obs_clusters) != float:
                fd.write('{}\t{}\t{}\t{}\t{:.4f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\t{}\t{}\n'.format(
                   gene_name, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                    obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))
            else:
                fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gene_name, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                    obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))
    # Sort
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['E_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)

    return output_file


def write_cluster_results(genome, results, dir, file, sorter):
    """Save results to the output file
    :param genome: reference genome
    :param results: dict, dictionary of results, keys are element's symbols
    :param dir: str, output directory
    :param file: str, output file, if elements in elements file
    :param sorter: list, element symbols ranked by elements p-value
    :return: None
    """

    sorterindex = dict(zip(sorter, range(len(sorter))))
    output_file = dir+'/clusters_'+file+'.tsv'
    header = ['RANK','SYMBOL', 'CGC', 'CLUSTER', 'N', 'LEFT_M', 'MAX', 'RIGHT_M', 'WIDTH', 'MUT', 'SCORE']
    with open(output_file, 'w') as fd:
        fd.write('{}\n'.format('\t'.join(header)))
        for element, values in results.items():
            clustersinfo, cgc = values
            if genome != 'hg19':
                cgc = 'Non Available'
            if clustersinfo and type(clustersinfo) != float:
                rank = sorterindex[element]+1
                for c, v in clustersinfo.items():
                    fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            rank, element, cgc, v['mode'], c, v['left_m'][0], v['max'][0], v['right_m'][0], abs(v['right_m'][0]-v['left_m'][0]), v['n_mutations'], v['score']))
    # Sort
    df = pd.read_csv(output_file, sep='\t', header=0)
    df.sort_values(by=['RANK', 'SCORE'], ascending=[True, False], inplace=True)
    df.to_csv(path_or_buf=output_file, sep='\t', na_rep='', index=False)


def write_info(input_file,
    output_directory,
    regions_file,
    genome,
    elements_file,
    elements,
    element_mutations,
    cluster_mutations,
    smooth_window,
    cluster_window,
    cluster_score,
    element_score,
    kmer,
    n_simulations,
    simulation_mode,
    simulation_window,
    cores,
    seed,
    log_level):
    """Write info to output
        :param input_file: input file
        :param output_directory: output directory path
        :param regions_file: path input genomic regions, tab file
        :param genome: genome to use
        :param elements_file: file containing one element per row
        :param elements: element symbol or file containing elements
        :param element_mutations: int, cutoff of element mutations
        :param cluster_mutations: int, cutoff of cluster mutations
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
        :param log_level: verbosity of the logger
        :return: None
    """

    info_file = output_directory + '/' + output_directory.split('/')[-1] + '.info'

    if not os.path.isfile(info_file):
        with open(info_file, 'w') as fd:
            fd.write('input_file: {}\noutput_directory: {}\nregions_file: {}\ngenome: {}\nelements_file: {}\n'
                     'elements: {}\nelement_mutations: {}\ncluster_mutations: {}\nsmooth_window: {}\n'
                     'cluster_window: {}\ncluster_score: {}\nelement_score: {}\n'
                     'kmer: {}\nn_simulations: {}\n'
                     'simulation_mode: {}\nsimulation_window: {}\ncores: {}\nseed: {}\n'
                     'log_level: {}\n'.format(input_file,output_directory,regions_file,genome,
                      elements_file,elements,element_mutations,cluster_mutations,smooth_window,cluster_window,
                      cluster_score,element_score,kmer,n_simulations,simulation_mode, simulation_window,
                      cores,seed,log_level))
    else:
        pass

@click.command()
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations')
@click.option('-o', '--output-directory', default=None, required=True,
              help='Output directory')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='File with the genomic regions to analyze')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbol of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element to analyze')
@click.option('-g', '--genome', default='hg19', type=click.Choice(['hg19', 'mm10', 'c3h']),
              help="genome to use")
@click.option('-emut', '--element-mutations', type=click.INT, default=2,
              help='Cutoff of element mutations. Default is 2')
@click.option('-cmut', '--cluster-mutations', type=click.INT, default=2,
              help='Cutoff of cluster mutations. Default is 2')
@click.option('-sw', '--smooth-window', type=click.INT, default=50,
              help='smoothing window. Default is 50')
@click.option('-cw', '--cluster-window', type=click.INT, default=50,
              help='cluster window. Default is 50')
@click.option('-cs', '--cluster-score', default='nobias', help='Cluster score formula',
              type=click.Choice(['nobias', 'fmutations']))
@click.option('-es', '--element-score', default='mean', help='Element score formula',
              type=click.Choice(['sum', 'mean']))
@click.option('-kmer', '--kmer', default='3', help='Number of nucleotides of the signature',
              type=click.Choice(['3', '5']))
@click.option('-n', '--n-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-sim', '--simulation-mode', default='element', help='Simulation mode',
              type=click.Choice(['segment', 'hotspot', 'element']))
@click.option('-simw', '--simulation-window', type=click.INT, default=20,
              help='Simulation window. Default is 20')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('--seed', type=click.INT, default=None,
              help='seed to use in the simulations')
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file,
         output_directory,
         regions_file,
         elements_file,
         elements,
         genome,
         element_mutations,
         cluster_mutations,
         smooth_window,
         cluster_window,
         cluster_score,
         element_score,
         kmer,
         n_simulations,
         simulation_mode,
         simulation_window,
         cores,
         seed,
         log_level):
    """Oncodriveclustl is a program that looks for mutational hotspots
    :param input_file: input file
    :param output_directory: output directory
    :param regions_file: path input genomic regions, tab file
    :param elements_file: file containing one element per row
    :param elements: element symbol or file containing elements
    :param genome: genome to use
    :param element_mutations: int, cutoff of element mutations
    :param cluster_mutations: int, cutoff of cluster mutations
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
    :param log_level: verbosity of the logger
    :return: None
    """

    # Get output directory
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Get file name
    if elements_file is not None:
        output_file = elements_file.split('/')[-1]
    else:
        output_file = 'results'

    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(filename=output_file+'.log', directory=output_directory)
    ))
    logger = daiquiri.getLogger()
    logger.debug(' '.join([input_file,
                           output_directory,
                           regions_file,
                           str(element_mutations),
                           str(cluster_mutations),
                           str(smooth_window),
                           str(cluster_window),
                           cluster_score,
                           element_score,
                           kmer,
                           str(n_simulations),
                           simulation_mode,
                           str(simulation_window),
                           str(cores)]))

    logger.info('Initializing OncodriveCLUSTL...')

    # Create a list of elements to analyze
    if elements is not None:
        elements = set(elements)
    if elements_file is not None:
        elements |= set([line.strip().split()[0] for line in open(elements_file, 'r')])
    if elements is None and elements_file is None:
        elements = set([])
    if elements:
        logger.info(
            'Input element{}: {}'.format('s' if len(elements) > 1 else '', len(elements))
        )
        logger.info(', '.join(elements))

    # Parse regions and dataset mutations
    logger.info('Parsing input regions and input mutations...')
    regions_d, chromosomes_d, mutations_d, gz = pars.parse(regions_file, elements, input_file)

    mut = 0
    elem = 0
    for k, v in mutations_d.items():
        mut += len(v)
        elem += 1
    logger.info('Validated input elements: {}'.format(len(regions_d.keys())))
    logger.info('Validated mutated elements to analyze: {}'.format(elem))
    logger.info('Substitution mutations to analyze: {}'.format(mut))

    # Compute dataset kmer signatures
    signatures_pickle = input_file.split('/')[-1][:-4] + '_' + kmer + '.pickle'
    path_cache = output_directory + '/cache'
    os.makedirs(path_cache, exist_ok=True)
    path_pickle = path_cache + '/' + signatures_pickle

    if not os.path.isfile(path_pickle):
        logger.info('Computing signatures...')
        obj = sign.Signature(start_at_0=True, genome=genome, kmer=int(kmer), log_level=log_level)
        obj.calculate(input_file, gz)
        obj.save(path_pickle)
        logger.info('Signatures computed')
    else:
        logger.info('Signatures computed')

    # Initialize Experiment class variables and run
    elements_results, clusters_results = exp.Experiment(
                                regions_d, chromosomes_d, mutations_d, genome,
                                path_pickle,
                                element_mutations, cluster_mutations,
                                smooth_window, cluster_window,
                                cluster_score, element_score,
                                int(kmer),
                                n_simulations, simulation_mode, simulation_window,
                                cores, seed
                                ).run()

    # Write results
    elements_path = write_element_results(genome=genome, results=elements_results, dir=output_directory, file=output_file)
    logger.info('Element results calculated')
    df = pd.read_csv(elements_path, sep='\t')
    sorted_list_elements = df['SYM'].tolist()
    write_cluster_results(genome=genome, results=clusters_results, dir=output_directory, file=output_file,
                          sorter=sorted_list_elements)
    logger.info('Clusters results calculated')
    logger.info('Finished')

    # Write info
    write_info(input_file,output_directory,regions_file,genome,elements_file,elements,element_mutations,
         cluster_mutations,smooth_window,cluster_window,cluster_score,element_score,int(kmer),
         n_simulations,simulation_mode,simulation_window,cores,seed,log_level)

if __name__ == '__main__':
    main()
