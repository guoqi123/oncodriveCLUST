"""
Contains the command line parsing
"""
import logging
import os
import csv

import click
import daiquiri

from oncodriveclustl.utils import exceptions as excep
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
@click.option('-sign', '--input-signature', default=None, required=False, type=click.Path(exists=True),
              help='File containing input context based mutational probabilities (signature)')
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
              type=click.Choice(['fmutations', 'cmutcorrected']))
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
@click.option('--protein', is_flag=True, help='Analyze clustering in translated protein sequences (CODING)')
@click.option('--plot', is_flag=True, help='Generate a clustering plot for an element')
@click.option('--pancancer', is_flag=True, help='PanCancer cohort analysis')
def main(input_file,
         input_signature,
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
         protein,
         plot,
         pancancer):
    """
    OncodriveCLUSTL is a sequence based clustering method to identify cancer drivers across the genome

    Args:
        input_file (str): path to mutations file
        input_signature (str): path to file containing input context based mutational probabilities (optional).
            By default (when no input signatures), OncodriveCLUSTL will calculate them from the mutations input file.
        output_directory(str): path to output directory. Output files will be generated in it.
        regions_file (str): path to input genomic coordinates file
        elements_file (str): path to file containing one element per row (optional) to analyzed the listed elements.
            By default, OncodriveCLUSTL analyzes all genomic elements contained in `regions_file`.
        elements (str): genomic element symbol (optional). The analysis will be performed only on the specified element.
        genome (str): genome to use: 'hg19', 'mm10', 'c3h', 'car', 'cast' and 'f344'
        element_mutations (int): minimum number of mutations per genomic element to undertake analysis
        cluster_mutations (int): minimum number of mutations to define a cluster
        smooth_window (int): Tukey kernel smoothing window length
        cluster_window (int): clustering window length
        cluster_score (str): cluster score method
        element_score (str): element score method
        kmer (int): context nucleotides to calculate the mutational probabilities (trinucleotides or pentanucleotides)
        n_simulations (int): number of simulations
        simulation_mode (str): simulation mode
        simulation_window (int): window length to simulate mutations
        cores (int): number of CPUs to use
        log_level (str): verbosity of the logger
        gzip (bool): flag to generate GZIP compressed output files
        cds (bool): flag to calculate clustering on collapsed genomic regions (e.g., coding regions in a gene)
        conseq (bool): flag to use only non-synonymous mutations for clustering analysis
        protein (bool): flag to analyze clustering in translated protein sequences
        plot (bool): flag to generate a clustering plot for an element
        pancancer (bool): flag to compute PanCancer cohort analysis

    Returns:
        None

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
        'input_signature: {}'.format(input_signature),
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
        'VEP conseq: {}'.format(conseq),
        'protein clustering: {}'.format(protein),
        'pancancer: {}'.format(pancancer),
        ''
    ]))
    logger.info('Initializing OncodriveCLUSTL...')

    # Check parameters
    if n_simulations < 1000:
        raise excep.UserInputError('Invalid number of simulations: please choose integer greater than 1000')

    if conseq and cds is False:
        raise excep.UserInputError('Analysis using mutations consequence type requires analysis mode "--cds"')

    if protein and cds is False:
        raise excep.UserInputError('Analysis in translated protein sequences requires analysis mode "--cds"')

    # If --plot, only one element is analyzed
    if plot:
        if len(elements) != 1:
            raise excep.UserInputError('Plot can only be calculated for one element')
        if not cds:
            raise excep.UserInputError('Plots are only available for analysis mode "--cds"')

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
    if not input_signature:
        read_function, mode, delimiter, cancer_type = prep.check_tabular_csv(input_file)
        path_cache = os.path.join(output_directory, 'cache')
        os.makedirs(path_cache, exist_ok=True)
        obj = sign.Signature(start_at_0=True, genome=genome, kmer=int(kmer), log_level='info', pancancer=pancancer)
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
                # TODO warning if len == 1?
                # TODO threshold number of mutations per cohort to compute signature?
                logger.info(
                    'PanCancer analysis for {} cohort{}'.format(
                        len(cohorts_of_analysis), 's' if len(cohorts_of_analysis) > 1 else ''
                    )
                )
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
                raise excep.UserInputError('PanCancer analysis requires "CANCER_TYPE" column in input file')
        else:
            # Check if signatures computed
            file_prefix = input_file.split('/')[-1].split('.')[0]
            path_pickle = os.path.join(path_cache, '{}_kmer_{}.pickle'.format(file_prefix, kmer))
            if os.path.isfile(path_pickle):
                logger.info('Signatures computed')  # TODO check format
            # Calculate signatures
            else:
                logger.info('Computing signatures...')
                obj.calculate(input_file)
                obj.save(path_cache, prefix=file_prefix)
                logger.info('Signatures computed')
            cohorts_of_analysis.add(file_prefix)
    else:
        # Check format
        path_cache = input_signature
        error = prep.check_signature(input_signature, kmer)
        if error:
            raise excep.UserInputError('Signatures file could not be read. Please check {}'.format(input_signature))
        if pancancer:
            raise excep.UserInputError('PanCancer analysis computes one signature file for each cancer type present '
                                       'in the mutations input file, according to "CANCER_TYPE" column.')
        logger.info('Input signatures loaded')

    # Parse regions and dataset mutations
    logger.info('Parsing genomic regions and mutations...')
    regions_d, cds_d, chromosomes_d, strands_d, mutations_d, samples_d, cohorts_d = pars.parse(
        regions_file,
        elements,
        input_file,
        cds,
        conseq,
        protein,
        genome
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
    logger.info('Validated {} in genomic regions: {}'.format(
        'elements' if not protein else 'transcripts', len(regions_d.keys())
    ))
    logger.info('Validated {} with mutations: {}'.format('elements' if not protein else 'transcripts', elem))
    logger.info('Total substitution mutations: {}'.format(mut))
    if not element_mutations_cutoff:
        raise excep.UserInputError('No element found with enough mutations to perform analysis')

    # Initialize Experiment class variables and run
    elements_results, clusters_results = exp.Experiment(
        regions_d,
        cds_d,
        chromosomes_d,
        strands_d,
        mutations_d,
        samples_d,
        genome,
        cohorts_d,
        path_cache,
        element_mutations,
        cluster_mutations,
        smooth_window,
        cluster_window,
        cluster_score,
        element_score,
        int(kmer),
        n_simulations,
        simulation_mode,
        simulation_window,
        cores,
        conseq,
        protein,
        plot
    ).run()

    # Write elements results (returns list of ranked elements)
    sorted_list_elements = postp.write_element_results(genome=genome,
                                                       results=elements_results,
                                                       directory=output_directory,
                                                       file=elements_output_file,
                                                       is_gzip=gzip)
    logger.info('Elements results calculated')

    # Write clusters results
    postp.write_cluster_results(genome=genome,
                                results=clusters_results,
                                directory=output_directory,
                                file=clusters_output_file,
                                sorter=sorted_list_elements,
                                regions_d = regions_d,
                                is_gzip=gzip,
                                is_protein=protein)
    logger.info('Clusters results calculated')
    logger.info('Finished')


if __name__ == '__main__':
    main()
