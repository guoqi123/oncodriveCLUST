# Import modules

import tabix
from concurrent.futures import ProcessPoolExecutor as Pool
from functools import partial


class Query:
    tb = tabix.open('/workspace/datasets/phd_snp_g/input_files_cds/vep_canonical.tsv.gz')

    def query_tabix(self, position, chromosome):
        res = []
        for i in self.tb.query(chromosome, position - 1, position):
            res.append(i)
        return res

    def parallel_queries(self, chromosome, positions):
        with Pool(2) as pool:
            fx = partial(self.query_tabix, chromosome=chromosome)
            for i in pool.map(fx, positions):
                print(i)