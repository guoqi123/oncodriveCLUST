# Import modules
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import click
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt


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


def tab_files_iterator():
    """Return path and subdirectories"""

    rtab = re.compile('.*tab')
    rtxt = re.compile('.*txt')

    for entry in os.scandir(self.directory):
        if entry.is_dir(follow_symlinks=False):
            if not entry.name.startswith('.') and entry.path.split('/')[-1] != 'cache':

                # Get list of content in path
                content = os.listdir(entry.path)

                # Search for .tab and .txt files
                tab_files = list(filter(rtab.match, content))
                txt_files = list(filter(rtxt.match, content))

                if txt_files and not tab_files:
                    merge(entry.path)
                    print('Results merged in: ', entry.path)

                tab_files_iterator(entry.path)


def mtc(p_value):
    """
    Calculate q-value
    :param p_value: array
    :return: q_value
    """
    return mlpt(p_value, alpha=0.05, method='fdr_bh')[1]


def qval_format(input, output, pvalue):
    """
    File format adding q-value
    :param file: file to format in place
    :param pvalue: reference p-value, empirical or analytical
    :return:
    """

    df = pd.read_csv(input, sep='\t', header=0)
    q_values = pd.DataFrame(mtc(df.loc[:, pvalue]))
    df['QVAL'] = q_values
    df.sort_values(by=['QVAL', 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=False)
    df.to_csv(path_or_buf=output, sep='\t', na_rep='', index=False)


def get_qvalue(directory):
    """

    :param directory:
    :return:
    """

    rtab = re.compile('.*tab')

    for entry in os.scandir(directory):
        if entry.is_dir(follow_symlinks=False):
            if not entry.name.startswith('.') and entry.path.split('/')[-1] != 'cache':

                # Get list of content in path
                content = os.listdir(entry.path)

                # Search for .tab files
                tab_files = list(filter(rtab.match, content))

                if tab_files:
                    for file in tab_files:
                        #  Get if its empirical or analytical pvalue
                        if file.endswith('elements_empirical.tab'):
                            pval = 'E_PVAL'
                        elif file.endswith('elements_analytical.tab'):
                            pval = 'A_PVAL'
                        else:
                            pval = ''

                        if pval != '':
                            # Get paths
                            input_path = entry.path + '/' + file
                            output_path = entry.path + '/' + file[:-4] + '_qvalues.tab'
                            # Calculate q-values
                            qval_format(input=input_path, output=output_path, pvalue=pval)
                            print('Q-values calculated in: ', output_path)

                get_figures(entry.path)


def qqplot(df, col_values, top=10, output=None):
    """
    """
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot2grid((1, 1), (0, 0))

    df = df[np.isfinite(df[col_values])].copy()
    df.sort_values(by=[col_values, 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)

    min_pvalue = min([x for x in df[col_values].tolist() if x > 0])
    obs = list(df[col_values].map(lambda x: -np.log10(x) if x > 0 else -np.log10(min_pvalue)))

    data = pd.DataFrame({
        'SYM': df['SYM'],
        'SCORE_OBS': df['SCORE_OBS'],
        'CGC': df['CGC'],
        'obs': obs,
    })

    data.sort_values(by=['obs', 'SCORE_OBS', 'CGC'], ascending=[True, True, False], inplace=True)
    exp = -np.log10(np.arange(1, len(data) + 1) / len(data))
    exp.sort()
    data['exp'] = exp

    # Scatter of the genes
    ax.scatter(data['exp'], data['obs'], color='blue', alpha=0.5)

    # Null hypothesis
    ax.plot(data['exp'], data['exp'], color='red', linestyle='--')

    # Significant genes
    sign = data[-top:].copy()

    x_text = ax.get_xlim()[1] * 1.1
    y_text = max(data['obs'])
    delta = y_text / 20

    sign.sort_values(by='exp', ascending=False, inplace=True)
    for count, line in sign.iterrows():
        color = 'red' if line['CGC'] is True else 'black'
        ax.annotate(line['SYM'], xy=(line['exp'], line['obs']),
                    xytext=(x_text, y_text), color=color,
                    arrowprops=dict(color='grey', shrink=0.05, width=0.5, headwidth=5, alpha=0.20), )
        y_text -= delta

    ax.set_xlabel(r"Expected pvalues $(-\log_{10})$", fontsize=16)
    ax.set_ylabel(r"Observed pvalues $(-\log_{10})$", fontsize=16)
    ax.set_title(r"Empirical raw p-values", fontsize=18) if col_values == 'E_PVAL' else ax.set_title(r"Analytical raw p-values", fontsize=18)

    if output:
        plt.savefig(output, bbox_inches='tight')


def histogram(df, col_values, output=None):
    """
    """
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot2grid((1, 1), (0, 0))

    df = df[np.isfinite(df[col_values])].copy()
    df.sort_values(by=[col_values, 'SCORE_OBS', 'CGC'], ascending=[True, False, False], inplace=True)

    # Histogram
    plt.hist(df[col_values], edgecolor='black', linewidth=0.5)
    ax.set_xlabel(r"p-values", fontsize=16)
    ax.set_ylabel(r"Frequency", fontsize=16)
    ax.set_title(r"Empirical raw p-values", fontsize=18) if col_values == 'E_PVAL' else ax.set_title(r"Analytical raw p-values", fontsize=18)

    if output:
        plt.savefig(output, bbox_inches='tight')


def get_figures(directory):
    """Return path and subdirectories"""

    rtab = re.compile('.*tab')

    for entry in os.scandir(directory):
        if entry.is_dir(follow_symlinks=False):
            if not entry.name.startswith('.') and entry.path.split('/')[-1] != 'cache':

                # Get list of content in path
                content = os.listdir(entry.path)

                # Search for .tab files
                tab_files = list(filter(rtab.match, content))
                figures = 'figures' in content

                if not figures and tab_files:
                #if tab_files:
                    figures_path = entry.path + '/figures'
                    os.makedirs(figures_path)

                    for file in tab_files:
                        #  Get if its empirical or analytical pvalue
                        if file.endswith('elements_empirical.tab'):
                            pval = 'E_PVAL'
                        elif file.endswith('elements_analytical.tab'):
                            pval = 'A_PVAL'
                        else:
                            pval = ''

                        if pval != '':
                            # Get paths
                            input_path = entry.path+'/'+file
                            output_path = figures_path+'/'+file[:-4]
                            # Save figures
                            df = pd.read_csv(input_path, sep='\t', header=0)
                            qqplot(df, col_values=pval, output=output_path+'_qqplot.png')
                            histogram(df, col_values=pval, output=output_path+'_histogram.png')

                            print('Figures saved in: ', output_path)

                get_figures(entry.path)

@click.command()
@click.option('-o', '--parent-directory', default=None, required=True, type=click.Path(exists=True),
              help='Parent directory containing result files')
def main(parent_directory):
    """Merge results from oncodriveclustl
    :param parent_directory: parent output directory path"""

    get_tab_files(parent_directory)
    #get_qvalue(parent_directory)
    get_figures(parent_directory)


if __name__ == '__main__':
    main()
