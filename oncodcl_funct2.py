# Oncodriveclustl functions 2
import oncodcl_funct as odf


def run_region(arguments):
    """
    Given gene, its genomic regions and its mutations, calculate the smooth, clusters
    and gene score.
    :param arguments: tuple, (gen, genomic regions, mutations)
            element: str
            genomic regions: dict
            mutations: list
    :return: tuple, (gene, length, n clusters, total mutations, score)
    """
    gene, regions, mutations = arguments

    region_info = odf.smoothing(symbol=gene, regions=regions, mutations=mutations, window=50)
    clusters = odf.clustering(regions=region_info, mutations=mutations, window=50)
    n_clusters, gene_score = odf.score_gene(clusters=clusters, cutoff=2)
    return (gene, len(region_info['genomic']), n_clusters, len(mutations), gene_score)
