# Import modules
import os
import re

import click
import pandas as pd


def read_write_elements(directory, efile, afile):
    """
    Read elements results and merge them in a single file.
    :param directory: output directory
    :param efile: output file sorted by empirical p-value
    :param afile: output file sorted by analytical p-value
    :return: None
    """

    header = ['SYM', 'LEN', 'N_MUT',
              'CLU', 'MEAN_SIM_CLU', 'MEDIAN_SIM_CLU', 'SD_SIM_CLU',
              'SCORE_OBS', 'MEAN_SIM_SCORE', 'MEDIAN_SIM_SCORE', 'SD_SIM_SCORE',
              'E_PVAL', 'A_PVAL', 'CGC']
    with open(efile, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for file in os.listdir(directory):
            if file.endswith(".txt"):
                file = directory + '/' + file
                with open(file, 'r') as fd:
                    next(fd)  # don't read header
                    for line in fd:
                        sym, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters, obs_score, \
                        mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc = line.strip().split('\t')
                        ofd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            sym, length, muts, obs_clusters, mean_sim_clusters, median_sim_clusters, std_sim_clusters,
                            obs_score, mean_sim_score, median_sim_score, std_sim_score, epval, apval, cgc))

    # Format empirical p-value
    df = pd.read_csv(efile, sep='\t', header=0)
    df.sort_values(by=['E_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    df.to_csv(path_or_buf=efile, sep='\t', na_rep='', index=False)

    # Format analytical p-value
    df = pd.read_csv(efile, sep='\t', header=0)
    df.sort_values(by=['A_PVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)
    df.to_csv(path_or_buf=afile, sep='\t', na_rep='', index=False)


def read_write_clusters(directory, elements_file, clusters_file):
    """
    Read cluster results and merge them in a single file according to elements ranking of empirical or
    analytical p-value.
    :param directory: path, results directory
    :param elements_file: input elements file
    :param clusters_file: output clusters file sorted according to elements_file
    :return: None
    """

    # Get sorter
    df = pd.read_csv(elements_file, sep='\t')
    sorter = df['SYM'].tolist()
    sorterindex = dict(zip(sorter, range(len(sorter))))

    header = ['RANK','SYMBOL', 'CLUSTER', 'N', 'LEFT_M', 'MAX', 'RIGHT_M', 'WIDTH', 'MUT', 'SCORE', 'CGC']
    with open(clusters_file, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for file in os.listdir(directory):
            if file.endswith(".tsv"):
                file = directory + '/' + file
                with open(file, 'r') as fd:
                    next(fd)  # don't read header
                    for line in fd:
                        _, sym, cluster, n, left_m, max, right_m, width, mut, score, cgc = line.strip().split('\t')
                        rank = sorterindex[sym]
                        ofd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            rank, sym, cluster, n, left_m, max, right_m, width, mut, score, cgc))

    # Sort
    df = pd.read_csv(clusters_file, sep='\t', header=0)
    df.sort_values(by=['RANK', 'SCORE'], ascending=[True, False], inplace=True)
    df.to_csv(path_or_buf=clusters_file, sep='\t', na_rep='', index=False)


def merge(directory):
    """Merge results from a given directory
    :param directory: directory to merge results into .tab files"""

    # Files to write
    experiment = directory + '/' + directory.split('/')[-1]
    ef_empirical = experiment + '_elements_empirical.tab'
    ef_analytical = experiment + '_elements_analytical.tab'
    cf_empirical = experiment + '_clusters_empirical.tab'
    cf_analytical = experiment + '_clusters_analytical.tab'

    read_write_elements(directory=directory, efile=ef_empirical, afile=ef_analytical)
    read_write_clusters(directory=directory, elements_file=ef_empirical, clusters_file=cf_empirical)
    read_write_clusters(directory=directory, elements_file=ef_analytical, clusters_file=cf_analytical)

def get_tree(directory):
    """Return path and subdirectories"""

    rtab = re.compile('.*tab')
    rtxt = re.compile('.*txt')

    for entry in os.scandir(directory):
        if entry.is_dir(follow_symlinks=False):
            if not entry.name.startswith('.') and entry.path.split('/')[-1] != 'cache':

                # Get list of content in path
                content = os.listdir(entry.path)

                # Search for .tab and .txt files
                tab_files = list(filter(rtab.match, content))
                txt_files = list(filter(rtxt.match, content))

                if txt_files and not tab_files:
                    print('Results merged in: ', entry.path)
                    merge(entry.path)

                get_tree(entry.path)

@click.command()
@click.option('-o', '--output-directory', default=None, required=True, type=click.Path(exists=True),
              help='Directory containing result files')
def main(output_directory):
    """Merge results from oncodriveclustl
    :param output_directory: output directory path"""

    get_tree(output_directory)
    #merge(output_directory)

if __name__ == '__main__':
    main()
