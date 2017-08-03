# Oncodriveclustl run
import sys
from concurrent.futures import ProcessPoolExecutor as Pool
from tqdm import tqdm
from collections import defaultdict

import oncodcl_funct as odf
import oncodcl_funct2 as odf3

def main(regions_d, mutations_d, cores=4):
    """Main function
    :param regions_d: dictionary of dictionaries containing genomic regions
    :param mutations_d: dictionary of lists, 'gene' : [mutations]
    :param cores: int, number of CPUs to use
    :return: None
    """
    results_d = defaultdict()

    # Analyze only genes with >= 2 mutations
    genes = [(g, regions_d[g], mutations_d[g]) for g in regions_d.keys() if len(mutations_d[g]) >= 2]

    with Pool(max_workers=cores) as executor, tqdm(total=len(genes)) as pbar:
        for result in executor.map(odf3.run_region, genes):
            pbar.update(1)
            # 'gene' : (length, n clusters, total mutations, score)
            results_d[result[0]] = (result[1], result[2], result[3], result[4])

        print(len(results_d.keys()))
        sys.stderr.write('Results calculated\n')


# Parse regions
input_regions = '/home/carnedo/projects/oncodriveclustl/inputs/regions/02_cds.regions.gz'
trees, regions_d = odf.regions(input_regions)
sys.stderr.write('Regions parsed\n')

# Read mutations, intersect with regions
input_mutations = '/home/carnedo/projects/oncodriveclustl/inputs/mutations/pancanatlas/SKCM.txt'
mutations_d = odf.read_mutations(input_mutations, trees)
sys.stderr.write('Mutations read\n')

if __name__ == '__main__':
    main(regions_d, mutations_d)
