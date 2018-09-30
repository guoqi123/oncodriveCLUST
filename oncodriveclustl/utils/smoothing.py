"""
Contains function to apply kde smooth of mutations
"""
import os
from intervaltree import IntervalTree
from collections import namedtuple

import numpy as np
import bgreference as bg
import json

Mutation = namedtuple('Mutation', 'position, region, alt, muttype, sample, cancertype')


def smooth_nucleotide(regions, concat_regions_d, mutations, tukey_filter, simulation_window):
    """Generate a smoothing curve for a list of element's mutations in the nucleotide sequence

    Args:
        regions (IntervalTree): IntervalTree with genomic positions of an element
        concat_regions_d (dict): keys are start genomic regions, values are positions (index) relative to the start
        mutations (list): list of mutations formatted as namedtuple
        tukey_filter (numpy.ndarray): kde array, length equals smoothing window.
        simulation_window (int): simulation window

    Returns:
        final_smooth_tree (IntervalTree): interval are genomic regions or indexes (concatenate mode),
            data np.array of smoothing score by position
        mutations_in (list): list of mutations in regions
    """
    first_smooth_tree = IntervalTree()
    final_smooth_tree = IntervalTree()
    mutations_in = []

    # Generate smoothing arrays for regions
    for interval in regions:
        # Add extra bases for smoothing of simulated mutations that fall outside regions and tukey_filter
        first_smooth_tree.addi(interval.begin, interval.end,
                         np.zeros((interval.end - interval.begin + len(tukey_filter) + simulation_window - 2)))

    if not concat_regions_d:
        # Smooth
        for mutation in mutations:
            for interval in first_smooth_tree[mutation.region[0]]:
                # Get index of mutation in region
                new_begin = interval.begin - (simulation_window + len(tukey_filter) - 2)//2   # always integer
                index = mutation.position - new_begin
                tukey_begin = index - (len(tukey_filter) - 1)//2
                # Smooth mutations
                interval.data[tukey_begin: tukey_begin + len(tukey_filter)] += tukey_filter
            # Get mutations inside regions
            if regions[mutation.position]:
                mutations_in.append(mutation)

        # Remove extra bp
        for interval in first_smooth_tree:
            begin = interval.begin
            end = interval.end
            slicer = (simulation_window + len(tukey_filter) - 2)//2
            final_smooth_tree.addi(begin, end, interval.data[slicer: - slicer])

    else:
        # Smooth simulated mutations outside regions
        for mutation in mutations:
            if mutation.muttype == 1:
                if not first_smooth_tree[mutation.position]:
                    for interval in first_smooth_tree[mutation.region[0]]:
                        new_begin = interval.begin - (simulation_window + len(tukey_filter) - 2) // 2  # always integer
                        index = mutation.position - new_begin
                        tukey_begin = index - (len(tukey_filter) - 1) // 2
                        # Smooth mutations
                        interval.data[tukey_begin: tukey_begin + len(tukey_filter)] += tukey_filter

        # Remove extra bp
        for interval in first_smooth_tree:
            begin = interval.begin
            end = interval.end
            slicer = (simulation_window + len(tukey_filter) - 2) // 2
            final_smooth_tree.addi(begin, end, interval.data[slicer: - slicer])

        # Merge sorted regions (one interval == concatenated sequence) and add tukey//2 to both ends
        concat_tree = IntervalTree()
        concat_array = np.zeros((len(tukey_filter) - 1) // 2)
        for interval in sorted(final_smooth_tree):
            concat_array = np.append(concat_array, interval.data)
        concat_array = np.append(concat_array, np.zeros((len(tukey_filter) - 1) // 2))
        concat_tree.addi(final_smooth_tree.begin(), final_smooth_tree.end(), concat_array)
        final_smooth_tree = IntervalTree()

        # Smooth mutations inside regions
        for mutation in mutations:
            if mutation.muttype == 1:
                if first_smooth_tree[mutation.position]:
                    for interval in concat_tree[mutation.position]:
                        # Get index of mutation in concatenated sequence
                        index = (mutation.position - mutation.region[0]) + concat_regions_d[mutation.region[0]].start
                        # Smooth mutations
                        interval.data[index: (index + len(tukey_filter))] += tukey_filter
                    mutations_in.append(mutation)

        # Remove extra bp
        for interval in concat_tree:
            begin = interval.begin
            end = interval.end
            slicer = (len(tukey_filter) -1) // 2
            final_smooth_tree.addi(begin, end, interval.data[slicer: -slicer])

    return final_smooth_tree, mutations_in


def smooth_aminoacid(regions, chromosome, strand, genome, tukey_filter, cds_d, mutations):
    """Generate a smoothing curve for a list of element's mutations in the aminoacid sequence (non-synonymous)

    Args:
        regions: IntervalTree with genomic positions of an element
        chromosome: str, chromosome
        strand: str, strand
        genome: str, genome
        cds_d: dict, keys are start genomic regions, values are cds positions
        mutations: list, list of mutations formatted as namedtuple
        tukey_filter: numpy array. Length equals smoothing window. The elements sum to 1

    Returns:
        final_smooth_tree, IntervalTree. Interval are genomic regions or cds, data np.array of smoothing score
        by position.
        mutations_in: list of mutations in regions
    """
    cds_aa_tree = IntervalTree()
    mutations_in = []
    reverse_cds_d = IntervalTree()
    for genomic, cds in cds_d.items():
        reverse_cds_d.addi(cds.start, cds.end + 1, genomic)
    # TODO: remove hardcoded file
    with open(os.path.join(os.path.dirname(__file__), '../data/genetic_code_ncbi_20180727_v1.json'), 'rt') as fd:
        genetic_code = json.load(fd)
    reverse_d = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    # Merge sorted regions (one interval == cds) and add tukey//2 to both ends
    cds_nucleotides = 0
    for interval in sorted(regions):
        cds_nucleotides += interval.end - interval.begin  # end is +1
    # Convert to aminoacid and add extra aa for regions in the first and last positions according to tukey_filter
    cds_aa_array = np.zeros((len(tukey_filter) - 1) + cds_nucleotides // 3)

    # Smooth mutations inside cds (aa sequence)
    for mutation in mutations:
        # print(mutation)
        for interval in regions[mutation.position]:
            # Only mutations inside cds region are considered
            if interval:
                nucleotide_cds_position = (mutation.position - mutation.region[0]) + cds_d[mutation.region[0]].start
                aa_protein_position = nucleotide_cds_position // 3  # triplets

                # Get mutation position in triplet
                triplet_position_proxi = round(nucleotide_cds_position / 3, 1)
                if triplet_position_proxi == (0.7 + aa_protein_position):
                    triplet_position = 2
                else:
                    triplet_position = 1 if triplet_position_proxi == (0.3 + aa_protein_position) else 0

                # Get nucleotide cds position for each nucleotide in the triplet
                elements_to_search = [i + (nucleotide_cds_position - triplet_position) for i in [0, 1, 2]]
                ref_triplet = []
                for cds_position in elements_to_search:  # One by one, splicing taken into account
                    # Get reference nucleotide
                    for i in reverse_cds_d[cds_position]:
                        genomic_position = i.data
                        ref_triplet.append(bg.refseq(genome, chromosome, genomic_position, 1))

                # Reverse if negative strand
                if strand == '-':
                    # Reverse
                    ref_triplet.reverse()
                    # Complementary
                    ref_triplet = [reverse_d.get(i, i) for i in ref_triplet]
                    alternate = reverse_d.get(mutation.alt)
                else:
                    alternate = mutation.alt

                # Get alternate triplet
                alt_triplet = ref_triplet.copy()
                alt_triplet[triplet_position] = alternate

                # Translate
                ref_aa = genetic_code[''.join(ref_triplet)][0]
                alt_aa = genetic_code[''.join(alt_triplet)][0]

                # print(nucleotide_cds_position, aa_protein_position, triplet_position, elements_to_search)
                # print('Ref', ref_triplet, ref_aa)
                # print('Alt', alt_triplet, alt_aa)

                # Only non-synonymous mutations are considered
                if alt_aa != ref_aa:
                    cds_aa_array[aa_protein_position: (aa_protein_position + len(tukey_filter))] += tukey_filter
                    aa_mutation = Mutation(aa_protein_position, (0, cds_nucleotides // 3), mutation.alt, 1, mutation.sample, mutation.cancertype)
                    mutations_in.append(aa_mutation)

    # Remove extra bp
    slicer = (len(tukey_filter) - 1) // 2
    cds_aa_array = cds_aa_array[slicer: -slicer]

    # New aa based tree
    cds_aa_tree.addi(0, cds_nucleotides // 3, cds_aa_array)

    return cds_aa_tree, mutations_in