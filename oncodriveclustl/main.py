# Oncodriveclustl run
import logging
import os
import sys

import click
import daiquiri

from oncodriveclustl.utils import signature as sign
from oncodriveclustl.utils import parsing as pars
from oncodriveclustl.utils import run as exp
from oncodriveclustl.utils import postprocessing as postp

# Global variables
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

@click.command()
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations')
@click.option('-vep', '--vep-file', default=None, required=False, type=click.Path(exists=True),
              help='File containing somatic mutations in vep format')
@click.option('-o', '--output-directory', default=None, required=True,
              help='Output directory')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='File with the genomic regions to analyze')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbol of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element to analyze')
@click.option('-g', '--genome', default='hg19', type=click.Choice(['hg19', 'mm10', 'c3h', 'car', 'cast']),
              help='Genome to use')
@click.option('-emut', '--element-mutations', type=click.INT, default=2,
              help='Cutoff of element mutations. Default is 2')
@click.option('-cmut', '--cluster-mutations', type=click.INT, default=2,
              help='Cutoff of cluster mutations. Default is 2')
@click.option('-sw', '--smooth-window', type=click.INT, default=25,
              help='Smoothing window. Default is 50')
@click.option('-cw', '--cluster-window', type=click.INT, default=50,
              help='Cluster window. Default is 50')
@click.option('-cs', '--cluster-score', default='fmutations', help='Cluster score formula',
              type=click.Choice(['fmutations']))
@click.option('-es', '--element-score', default='sum', help='Element score formula',
              type=click.Choice(['sum']))
@click.option('-kmer', '--kmer', default='3', help='Number of nucleotides of the signature',
              type=click.Choice(['3', '5']))
@click.option('-n', '--n-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-sim', '--simulation-mode', default='exon', help='Simulation mode',  # TODO change exon for region
              type=click.Choice(['exon', 'cds', 'exon_restricted']))
@click.option('-simw', '--simulation-window', type=click.INT, default=60,
              help='Simulation window. Default is 20')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('--seed', type=click.INT, default=None,
              help='seed to use in the simulations')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
@click.option('--gzip', is_flag=True, help='Gzip compress files')
@click.option('--cds', is_flag=True, help='Calculate clustering on coding DNA sequence (cds)',)
@click.option('--conseq', is_flag=True, help='Use mutations consequence type from VEP (CODING)',)
@click.option('--plot', is_flag=True, help='Generate a clustering plot for an element',)
@click.option('--oncohort', is_flag=True, help='Generate output file for OnCohortDrive',)


def main(input_file,
         vep_file,
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
         log_level,
         gzip,
         cds,
         conseq,
         plot,
         oncohort):
    """Oncodriveclustl is a program that looks for mutational hotspots
    :param input_file: input file
    :param vep_file: input vep file, optional
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
    :param gzip: bool, True generates gzip compressed output file
    :param cds: bool, True calculates clustering on cds
    :param conseq: bool, True uses consequence type for cds
    :param plot: bool, True generates a clustering plot for an element
    :param oncohort: bool, True generates output file for OncohortDrive
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
    global logger
    logger = daiquiri.getLogger()

    logger.info('\n'.join([
        '',
        'input_file: {}'.format(input_file),
        'vep: {}'.format(vep_file),
        'output_directory: {}'.format(output_directory),
        'regions_file: {}'.format(regions_file),
        'genome: {}'.format(genome),
        'element_mutations: {}'.format(element_mutations),
        'cluster_mutations: {}'.format(cluster_mutations),
        'cds: {}'.format(cds),
        'smooth_window: {}'.format(smooth_window),
        'cluster_window: {}'.format(cluster_window),
        'cluster_score: {}'.format(cluster_score),
        'element_score: {}'.format(element_score),
        'kmer: {}'.format(kmer),
        'simulation_mode: {}'.format(simulation_mode),
        'simulation_window: {}'.format(simulation_window),
        'n_simulations: {}'.format(n_simulations),
        'cores: {}'.format(cores),
        'gzip: {}'.format(gzip),
        'oncohort: {}'.format(oncohort),
        'VEP conseq: {}'.format(conseq),
        ''
    ]))

    logger.info('Initializing OncodriveCLUSTL...')

    # Check parameters
    if n_simulations < 1000:
        logger.error('Invalid number of simulations: please choose integer greater than 1000')
        sys.exit(1)

    if conseq and cds is False:
        logger.error('Analysis using mutations consequence type requires analysis mode "--cds"'.format(simulation_mode))
        sys.exit(1)

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

    # If --plot, only one element is analyzed
    if plot:
        if len(elements) != 1:
            logger.critical('Plot can only be calculated for one element')
            sys.exit(1)
        if not cds:
            logger.critical('Plots are only available for cds')
            sys.exit(1)

    # Compute dataset kmer signatures
    signatures_pickle = input_file.split('/')[-1][:-4] + '_' + kmer + '.pickle'
    path_cache = output_directory + '/cache'
    os.makedirs(path_cache, exist_ok=True)
    path_pickle = path_cache + '/' + signatures_pickle
    if not os.path.isfile(path_pickle):
        logger.info('Computing signatures...')
        obj = sign.Signature(start_at_0=True, genome=genome, kmer=int(kmer), log_level=log_level)
        obj.calculate(input_file)
        obj.save(path_pickle)
        logger.info('Signatures computed')
    else:
        logger.info('Signatures computed')

    # Parse regions and dataset mutations
    logger.info('Parsing genomic regions and mutations...')
    regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d = pars.parse(regions_file, elements,
                                                                                    input_file, cds, vep_file, conseq)
    mut = 0
    elem = 0
    element_mutations_cutoff = False
    for k, v in mutations_d.items():
        mut += len(v)
        elem += 1
        if not element_mutations_cutoff:
            if len(v) >= element_mutations:
                element_mutations_cutoff = True
    logger.info('Validated elements in genomic regions: {}'.format(len(regions_d.keys())))
    logger.info('Validated elements with mutations: {}'.format(elem))
    logger.info('Total substitution mutations: {}'.format(mut))
    if not element_mutations_cutoff:
        logger.critical('No element with enough mutations to perform analysis')
        sys.exit(1)

    # Initialize Experiment class variables and run
    elements_results, clusters_results = exp.Experiment(
                                regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, genome,
                                path_pickle,
                                element_mutations, cluster_mutations,
                                smooth_window, cluster_window,
                                cluster_score, element_score,
                                int(kmer),
                                n_simulations, simulation_mode, simulation_window,
                                cores, seed, conseq, plot
                                ).run()
    # Write results
    sorted_list_elements = postp.write_element_results(genome=genome, results=elements_results,
                                                       directory=output_directory, file=output_file, gzip=gzip)
    logger.info('Elements results calculated')
    postp.write_cluster_results(genome=genome, results=clusters_results, directory=output_directory, file=output_file,
                                sorter=sorted_list_elements, gzip=gzip, cds_d=cds_d)
    logger.info('Clusters results calculated')
    if oncohort:
        postp.write_oncohortdrive_results(mutations=input_file, directory=output_directory, file=output_file,
                                          regions_d=regions_d)
        logger.info('Oncohortdrive file generated')
    logger.info('Finished')

if __name__ == '__main__':
    main()
