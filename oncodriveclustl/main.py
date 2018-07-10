# Oncodriveclustl run
import logging
import os
import csv

import click
import daiquiri

from oncodriveclustl.utils import signature as sign
from oncodriveclustl.utils import parsing as pars
from oncodriveclustl.utils import run as exp
from oncodriveclustl.utils import preprocessing as prep
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
@click.option('-o', '--output-directory', default=None, required=True,
              help='Output directory to be created')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='GZIP compressed file with the genomic regions to analyze')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbol of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element to analyze')
@click.option('-g', '--genome', default='hg19', type=click.Choice(['hg19', 'mm10', 'c3h', 'car', 'cast', 'f344']),
              help='Genome to use')
@click.option('-emut', '--element-mutations', type=click.INT, default=2,
              help='Cutoff of element mutations. Default is 2')
@click.option('-cmut', '--cluster-mutations', type=click.INT, default=2,
              help='Cutoff of cluster mutations. Default is 2')
@click.option('-sw', '--smooth-window', type=click.INT, default=30,
              help='Smoothing window. Default is 30')
@click.option('-cw', '--cluster-window', type=click.INT, default=30,
              help='Cluster window. Default is 30')
@click.option('-cs', '--cluster-score', default='fmutations', help='Cluster score formula',
              type=click.Choice(['fmutations']))
@click.option('-es', '--element-score', default='sum', help='Element score formula',
              type=click.Choice(['sum']))
@click.option('-kmer', '--kmer', default='3', help='Number of nucleotides of the signature',
              type=click.Choice(['3', '5']))
@click.option('-n', '--n-simulations', type=click.INT, default=10000,
              help='number of simulations. Default is 10000')
@click.option('-sim', '--simulation-mode', default='mutation_centered', help='Simulation mode',
              type=click.Choice(['mutation_centered', 'region_restricted']))
@click.option('-simw', '--simulation-window', type=click.INT, default=45,
              help='Simulation window. Default is 45')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all the available cores.')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
@click.option('--gzip', is_flag=True, help='Gzip compress files')
@click.option('--cds', is_flag=True, help='Calculate clustering on coding DNA sequence (cds)')
@click.option('--conseq', is_flag=True, help='Use mutations consequence type from VEP (CODING)')
@click.option('--plot', is_flag=True, help='Generate a clustering plot for an element')
@click.option('--oncohort', is_flag=True, help='Generate output file for OnCohortDrive')
@click.option('--pancancer', is_flag=True, help='PanCancer cohort analysis')


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
         log_level,
         gzip,
         cds,
         conseq,
         plot,
         oncohort,
         pancancer):
    """OncodriveCLUSTL (MSc version) is a sequence based clustering method to identify cancer drivers across the genome
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
    :param log_level: verbosity of the logger
    :param gzip: bool, True generates gzip compressed output file
    :param cds: bool, True calculates clustering on cds
    :param conseq: bool, True uses consequence type for cds
    :param plot: bool, True generates a clustering plot for an element
    :param oncohort: bool, True generates output file for OncohortDrive
    :param pancancer: bool, True computes PanCancer cohort analysis
    :return: None
    """

    # Get output directory
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Get output files name
    if elements_file is not None:
        output_file = elements_file.split('/')[-1]
    else:
        output_file = 'results'
    if gzip:
        elements_output_file = 'elements_{}.txt.gz'.format(output_file)
        clusters_output_file = 'clusters_{}.tsv.gz'.format(output_file)
    else:
        elements_output_file = 'elements_{}.txt'.format(output_file)
        clusters_output_file = 'clusters_{}.tsv'.format(output_file)


    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(filename=output_file+'.log', directory=output_directory)
    ))
    global logger
    logger = daiquiri.getLogger()

    logger.info('OncodriveCLUSTL')
    logger.info('\n'.join([
        '',
        'input_file: {}'.format(input_file),
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
        'Pancancer: {}'.format(pancancer),
        ''
    ]))
    logger.info('Initializing OncodriveCLUSTL...')

    # Check parameters
    if n_simulations < 1000:
        logger.error('Invalid number of simulations: please choose integer greater than 1000')
        quit(-1)

    if conseq and cds is False:
        logger.error('Analysis using mutations consequence type requires analysis mode "--cds"'.format(simulation_mode))
        quit(-1)

    # If --plot, only one element is analyzed
    if plot:
        if len(elements) != 1:
            logger.critical('Plot can only be calculated for one element')
            quit(-1)
        if not cds:
            logger.critical('Plots are only available for cds')
            quit(-1)

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

    # Check format input file and calculate signature
    read_function, mode, delimiter, cancer_type = prep.check_tabular_csv(input_file)
    path_cache = os.path.join(output_directory, 'cache')
    os.makedirs(path_cache, exist_ok=True)
    obj = sign.Signature(start_at_0=True, genome='hg19', kmer=int(kmer), log_level='info', pancancer=pancancer)
    cohorts_of_analysis = set()

    if pancancer:
        # Check header
        if cancer_type:
            # Read file and get input cohorts
            pancancer_pickles = 0

            with read_function(input_file, mode) as read_file:
                fd = csv.DictReader(read_file, delimiter=delimiter)
                for line in fd:
                    cohorts_of_analysis.add(line['CANCER_TYPE'])
            logger.info('PanCancer analysis for {} cohort{}'.format(len(cohorts_of_analysis), 's' if len(cohorts_of_analysis) >1 else '')) # TODO warning if len == 1?

            # Check if signatures computed
            for cohort in cohorts_of_analysis:
                path_pickle = os.path.join(path_cache, '{}_kmer_{}.pickle'.format(cohort, kmer))
                if os.path.isfile(path_pickle):
                    pancancer_pickles += 1
            if len(cohorts_of_analysis) == pancancer_pickles:
                logger.info('Signatures computed')
            # Calculate signatures
            else:
                logger.info('Computing signatures...')
                obj.calculate(input_file)
                obj.save(path_cache, prefix=None)
                logger.info('Signatures computed')
        else:
            logger.critical('PanCancer analysis requires "CANCER_TYPE" column in input file')
            quit(-1)
    else:
        # Check if signatures computed
        file_prefix = input_file.split('/')[-1].split('.')[0]
        path_pickle = os.path.join(path_cache, '{}_kmer_{}.pickle'.format(file_prefix, kmer))
        if os.path.isfile(path_pickle):
            logger.info('Signatures computed')
        # Calculate signatures
        else:
            logger.info('Computing signatures...')
            obj.calculate(input_file)
            obj.save(path_cache, prefix=file_prefix)
            logger.info('Signatures computed')
        cohorts_of_analysis.add(file_prefix)

    # Parse regions and dataset mutations
    logger.info('Parsing genomic regions and mutations...')
    regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, cohorts_d = pars.parse(
        regions_file,
        elements,
        input_file,
        cds,
        conseq
    )
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
        quit(-1)

    # Initialize Experiment class variables and run
    elements_results, clusters_results = exp.Experiment(
                                regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, genome,
                                path_cache, cohorts_d,
                                element_mutations, cluster_mutations,
                                smooth_window, cluster_window,
                                cluster_score, element_score,
                                int(kmer),
                                n_simulations, simulation_mode, simulation_window,
                                cores, conseq, plot
                                ).run()

    # Write elements results (returns list of ranked elements)
    sorted_list_elements = postp.write_element_results(genome=genome,
                                                       results=elements_results,
                                                       directory=output_directory,
                                                       file=elements_output_file,
                                                       gzip=gzip)
    logger.info('Elements results calculated')
    # Write clusters results
    postp.write_cluster_results(genome=genome,
                                results=clusters_results,
                                directory=output_directory,
                                file=clusters_output_file,
                                sorter=sorted_list_elements,
                                gzip=gzip,
                                cds_d=cds_d)
    logger.info('Clusters results calculated')
    # Write Oncohortdrive results
    if oncohort:
        postp.write_oncohortdrive_results(mutations_file=input_file,
                                          directory=output_directory,
                                          file=clusters_output_file,
                                          regions_d=regions_d)
        logger.info('Oncohortdrive file generated')
    logger.info('Finished')

if __name__ == '__main__':
    main()
