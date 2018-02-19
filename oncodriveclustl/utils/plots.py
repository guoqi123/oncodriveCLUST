# Import modules
import logging
import os
import gzip

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import click
import daiquiri

from oncodriveclustl.utils import preprocessing as prep

# Global variables
logs = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def histogram(df, col_values, output=None):
    """
    """
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot2grid((1, 1), (0, 0))

    # Drop nan and sort by p-value and score
    df = df[np.isfinite(df[col_values])].copy()
    df.sort_values(by=[col_values, 'SCORE_OBS'], ascending=[True, False], inplace=True)

    # Histogram
    plt.hist(df[col_values], edgecolor='black', linewidth=0.5)
    ax.set_xlabel(r"p-values", fontsize=16)
    ax.set_ylabel(r"Frequency", fontsize=16)
    plt.suptitle("Empirical raw p-values", fontsize=18, fontweight='bold') if col_values == 'E_PVAL' else \
        plt.suptitle("Analytical raw p-values", fontsize=18, fontweight='bold')

    if output:
        # Add file name to figure
        ax.set_title(output.split('/')[-1], fontsize=16)

        # Define name
        if col_values == 'E_PVAL':
            output += 'histogram_empirical'
        else:
            output += 'histogram_analytical'
        # Save figure
        plt.savefig(output, bbox_inches='tight')


def qqplot(df, col_values, top=10, output=None):
    """
    QQ plot
    """
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot2grid((1, 1), (0, 0))

    df = df[np.isfinite(df[col_values])].copy()
    df.sort_values(by=[col_values, 'SCORE_OBS'], ascending=[True, False], inplace=True)

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

    # Scatter of the elements
    ax.scatter(data['exp'], data['obs'], color='blue', alpha=0.5)

    # Null hypothesis
    ax.plot(data['exp'], data['exp'], color='red', linestyle='--')

    # Significant elements
    sign = data[-top:].copy()

    x_text = ax.get_xlim()[1] * 1.05
    y_text = max(data['obs'])
    delta = y_text / 18

    sign.sort_values(by='exp', ascending=False, inplace=True)
    for count, line in sign.iterrows():
        ax.annotate(line['SYM'], xy=(line['exp'], line['obs']),
                    xytext=(x_text, y_text), color='black', fontsize=16,
                    arrowprops=dict(color='grey', shrink=0.05, width=0.5, headwidth=5, alpha=0.20), )
        y_text -= delta

    ax.set_xlabel(r"Expected pvalues $(-\log_{10})$", fontsize=16)
    ax.set_ylabel(r"Observed pvalues $(-\log_{10})$", fontsize=16)
    plt.suptitle("Empirical raw p-values", fontsize=18, fontweight='bold') if col_values == 'E_PVAL' else \
        plt.suptitle("Analytical raw p-values", fontsize=18, fontweight='bold')

    if output:
        # Add file name to figure
        ax.set_title(output.split('/')[-1], fontsize=16)

        # Define name
        if col_values == 'E_PVAL':
            output += 'qqplot_empirical'
        else:
            output += 'qqplot_analytical'
        plt.savefig(output, bbox_inches='tight')


@click.command()
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='Input path to file containing result')
@click.option('-o', '--output-directory', default=None, required=True,
              help='Output directory')
@click.option('--log-level', default='info', help='verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file, output_directory, log_level):
    """Generate plots from oncodriveclustl results
    :param input_file: input file
    :param output_directory: output directory
    :param log_level: verbosity of the logger
    """

    # Configure logger
    global logger
    daiquiri.setup(level=logs[log_level])
    logger = daiquiri.getLogger()

    # Get output directory
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Get output file name prefix
    output_prefix = input_file.split('.')[0]
    logger.info('Plots will be generated for:')
    logger.info(output_prefix.split('/')[-1])

    # Read df
    df = pd.read_csv(input_file, sep='\t', header=0)

    # Plots for analytical and empirical p-values
    logger.info('Starting plots...')
    for pval in ['A_PVAL', 'E_PVAL']:
        qqplot(df, col_values=pval, output=output_prefix)
        histogram(df, col_values=pval, output=output_prefix)

    logger.info('Plots finished')

if __name__ == '__main__':
    main()
