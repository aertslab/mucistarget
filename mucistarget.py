#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Calculate the impact of mutations on the removal or introduction of TF binding sites.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import argparse
import sys
import tempfile
import time

from collections import OrderedDict


def read_genes_filename(genes_filename):
    """
    Read genes from a file and store them in a set.

    :param genes_filename: File with gene names.
    :return: genes_set
    """

    genes_set = set()

    with open(genes_filename, 'r') as fh:
        for gene_line in fh:
            gene = gene_line.rstrip(' \r\n')

            if gene == '' or gene.startswith('#'):
                continue

            genes_set.add(gene)

    return genes_set


def read_motif_ids_filename(motif_ids_filename):
    """
    Read motif IDs from a file and store them in a set.

    :param motif_ids_filename: File with motif IDs.
    :return: motif_ids_set
    """

    motif_ids_set = set()

    with open(motif_ids_filename, 'r') as fh:
        for motif_id_line in fh:
            motif_id = motif_id_line.rstrip(' \r\n')

            if motif_id == '' or motif_id.startswith('#'):
                continue

            motif_ids_set.add(motif_id)

    return motif_ids_set


def read_tfs_filename(tfs_filename):
    """
    Read TFs from a file and store them in a set.

    :param tfs_filename: File with TFs.
    :return: tfs_set
    """

    tfs_set = set()

    with open(tfs_filename, 'r') as fh:
        for tf_line in fh:
            tf = tf_line.rstrip(' \r\n')

            if tf == '' or tf.startswith('#'):
                continue

            tfs_set.add(tf)

    return tfs_set


def get_all_mutations_that_overlap_with_regdoms_of_genes(vcf_mut_iterator, genes_set, keep_all_mutations=False):
    """
    Get all unique mutations that fall in a regulatory domain of one of the genes in the gene set.

    :param vcf_mut_iterator:
        Iterator that yields a VCFmut object:
            - mutations.VCFmut.from_vcf_file()
            - mutations.VCFmut.from_mut_ids_file()
            - mutations.VCFmut.from_bedlike_mut_ids_file()
    :param genes_set:
        Set of genes which regulatory domains will be used to only keep those mutations that fall in those domains.
        If set to None, no filtering will be applied.
    :param keep_all_mutations:
        Keep also mutations that do not overlap with a regulatory domain, when set to True.
    :return:
        (vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
         input_vcf_mut_ids,
         associated_genes_set)

        Where vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict:
            OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
            mutation to the TSS and TSS as values.
        Where input_vcf_mut_ids:
            unique mutation IDs in the vcf_mut_iterator.
        Where associated_genes_set:
            unique gene names which are associated with at least one mutation of the input and which pass the input gene
            set filter if specified.
    """

    input_vcf_mut_ids = set()
    vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict = OrderedDict()
    associated_genes_set = set()

    for vcf_mut in vcf_mut_iterator:
        if vcf_mut.mut_id not in input_vcf_mut_ids:
            # Store all unique mutation IDs.
            input_vcf_mut_ids.add(vcf_mut.mut_id)

            associated_genes_and_distance_to_tss_and_tss_dict = vcf_mut.get_associated_genes_and_distance_to_tss_and_tss()

            if genes_set:
                # If a set of genes was provided, only keep those associated genes for the mutation that appear in this
                # set of genes.
                associated_genes_and_distance_to_tss_and_tss_dict = (
                    {associated_genes_and_distance_to_tss_and_tss_key: associated_genes_and_distance_to_tss_and_tss_dict[
                        associated_genes_and_distance_to_tss_and_tss_key]
                     for associated_genes_and_distance_to_tss_and_tss_key in associated_genes_and_distance_to_tss_and_tss_dict
                     if associated_genes_and_distance_to_tss_and_tss_key in genes_set
                     }
                )

            if associated_genes_and_distance_to_tss_and_tss_dict:
                # Store all mutations (VCFmut object) and associated genes information in a ordered dict.
                vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict[vcf_mut] \
                    = associated_genes_and_distance_to_tss_and_tss_dict

                # Keep track of all associated genes for all mutations.
                associated_genes_set.update(set(associated_genes_and_distance_to_tss_and_tss_dict))
            elif keep_all_mutations:
                # Store the VCFmut object anyway, but without gene and TSS info, if keep_all_mutations is set to True.
                vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict[vcf_mut] = {'': [None, None]}

    return (vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
            input_vcf_mut_ids,
            associated_genes_set)


def write_mut_to_associated_gene_output(mut_to_associated_genes_output_filename,
                                        vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                                        log_fh=sys.stderr):
    """
    Write all mutations and their associated genes and distance of the mutation to the TSS to a file.

    :param mut_to_associated_genes_output_filename:
        Output filename which will contain mutation info and associated genes info.
    :param vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict:
        OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
        mutation to the TSS and TSS as values.
        This dictionary can be made with get_all_mutations_that_overlap_with_regdoms_of_genes.
    :param log_fh:
        File handle to which the progress information is written.
    :return:
    """

    import mutations

    print('Write mutations to associated gene output file: ', end='', file=log_fh)

    mut_to_associated_genes_start_time = time.time()

    with open(mut_to_associated_genes_output_filename, 'w') as mut_to_associated_genes_fh:
        # Sort TADs identifiers.
        tads_ids = (
            sorted(mutations.VCFmut.tads_start_end_array_per_chrom)
            if mutations.VCFmut.tads_start_end_array_per_chrom
            else None
        )

        # Write header to the output file.
        print('# chrom',
              'start',
              'reference',
              'mutation',
              'mutation type',
              'mutation ID',
              'associated gene',
              'distance to TSS',
              'in #TADs of {0:d}'.format(len(tads_ids)) if tads_ids else '',
              *['in ' + tads_id + ' TADs' for tads_id in tads_ids] if tads_ids else '',
              sep='\t',
              file=mut_to_associated_genes_fh)

        for vcf_mut, associated_genes_and_distance_to_tss_and_tss_dict in vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict.items():
            # Write to the output file.
            for associated_gene, distance_to_tss_and_tss in associated_genes_and_distance_to_tss_and_tss_dict.items():
                distance_to_tss, tss = distance_to_tss_and_tss

                if tads_ids:
                    tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list = [
                        1
                        if vcf_mut.tss_of_associated_gene_in_same_tad_as_mutation(
                            tads_id=tads_id,
                            tss=tss
                        )
                        else 0
                        for tads_id in tads_ids
                    ]

                print(vcf_mut,
                      associated_gene,
                      '{0:+}'.format(distance_to_tss) if distance_to_tss else '',
                      sum(tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list) if tads_ids else '',
                      *tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list if tads_ids else '',
                      sep='\t',
                      file=mut_to_associated_genes_fh)

    mut_to_associated_genes_end_time = time.time()

    print('{0:f} seconds.'.format(mut_to_associated_genes_end_time - mut_to_associated_genes_start_time),
          file=log_fh)


def calculate_and_write_clusterbuster_crm_and_motif_delta_scores(vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                                                                 motif_collection_version,
                                                                 motif_ids_set,
                                                                 clusterbuster_output_filename,
                                                                 min_clusterbuster_crm_score_threshold,
                                                                 log_fh=sys.stderr):
    """
    Calculate the MotifLocator scores for the wildtype and mutant FASTA sequence for each mutation and for each motif
    and write the result to a file.

    :param vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict:
        OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
        mutation to the TSS and TSS as values.
        This dictionary can be made with get_all_mutations_that_overlap_with_regdoms_of_genes.
    :param motif_collection_version:
        motif collection version
    :param motif_ids_set:
        set of motif IDs
    :param clusterbuster_output_filename:
        Output filename to which Cluster-Buster CRM and motif delta scores are written.
    :param min_clusterbuster_crm_score_threshold:
        Minimum Cluster-Buster CRM score threshold needed in wildtype or mutant to keep the result.
    :param log_fh:
        File handle to which the progress information is written.
    :return:
        nbr_of_mutations_which_pass_clusterbuster_crm_score_threshold:
            number of mutation that pass the Cluster-Buster CRM score threshold.
    """

    print('Score mutations with Cluster-Buster ...\n', file=log_fh)

    import motifsinfo
    import mutations
    import clusterbuster

    # Fill MotifsInfo class with content for the chosen motif collection version.
    motifsinfo.MotifsInfo.set_motif_collection_version(motif_collection_version=motif_collection_version)

    vcf_mut_ids_passing_clusterbuster_crm_score_threshold = set()

    nbr_motifs = len(motif_ids_set)

    # Sort TADs identifiers.
    tads_ids = (
        sorted(mutations.VCFmut.tads_start_end_array_per_chrom)
        if mutations.VCFmut.tads_start_end_array_per_chrom
        else None
    )

    with open(clusterbuster_output_filename, 'w') as clusterbuster_output_fh:
        # Write header to the output file.
        print('# chrom',
              'start',
              'reference',
              'mutation',
              'mutation type',
              'mutation ID',
              'associated gene',
              'distance to TSS',
              'motif ID',
              'motif name',
              'directly annotated TFs',
              'wildtype Cluster-Buster CRM score',
              'mutant Cluster-Buster CRM score',
              'delta Cluster-Buster CRM score',
              'wildtype Cluster-Buster motif score',
              'mutant Cluster-Buster motif score',
              'delta Cluster-Buster motif score',
              'wildtype consensus sequence',
              'mutant consensus sequence',
              'in #TADs of {0:d}'.format(len(tads_ids)) if tads_ids else '',
              *['in ' + tads_id + ' TADs' for tads_id in tads_ids] if tads_ids else '',
              sep='\t',
              file=clusterbuster_output_fh)

        for motif_idx, motif_id in enumerate(sorted(motif_ids_set)):
            print('  Scoring all mutations with Cluster-Buster for motif "{0:s}" ({1:d} of {2:d}): '.format(
                motif_id,
                motif_idx + 1,
                nbr_motifs),
                end='',
                file=log_fh)

            clusterbuster_start_time = time.time()

            try:
                clusterbuster_delta_scores = \
                    clusterbuster.calculate_clusterbuster_delta_scores(
                        vcf_muts=vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                        motif_id=motif_id,
                        min_crm_score_threshold=min_clusterbuster_crm_score_threshold
                    )
            except ValueError as e:
                print('\n\n  Invalid mutation: ' + e.message + '\n',
                      file=log_fh)
                sys.exit(1)

            # Write to the output file.
            for vcf_mut, clusterbuster_delta_score in clusterbuster_delta_scores.items():
                associated_genes_and_distance_to_tss_and_tss_dict = vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict[vcf_mut]

                for associated_gene, distance_to_tss_and_tss in associated_genes_and_distance_to_tss_and_tss_dict.items():
                    vcf_mut_ids_passing_clusterbuster_crm_score_threshold.add(vcf_mut.mut_id)

                    tfs = motifsinfo.MotifsInfo.get_tfs_for_motif(motif_id)

                    distance_to_tss, tss = distance_to_tss_and_tss

                    if tads_ids:
                        tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list = [
                            1
                            if vcf_mut.tss_of_associated_gene_in_same_tad_as_mutation(
                                tads_id=tads_id,
                                tss=tss
                            )
                            else 0
                            for tads_id in tads_ids
                        ]

                    print(vcf_mut,
                          associated_gene,
                          '{0:+}'.format(distance_to_tss) if distance_to_tss else '',
                          motif_id,
                          motifsinfo.MotifsInfo.get_motif_name(motif_id),
                          ';'.join(sorted(tfs) if tfs else ['']),
                          str(clusterbuster_delta_score.wt_crm_score),
                          str(clusterbuster_delta_score.mut_crm_score),
                          str(clusterbuster_delta_score.crm_delta_score),
                          str(clusterbuster_delta_score.wt_motif_score),
                          str(clusterbuster_delta_score.mut_motif_score),
                          str(clusterbuster_delta_score.motif_delta_score),
                          clusterbuster_delta_score.wt_consensus,
                          clusterbuster_delta_score.mut_consensus,
                          sum(tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list) if tads_ids else '',
                          *tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list if tads_ids else '',
                          sep='\t',
                          file=clusterbuster_output_fh)

            clusterbuster_end_time = time.time()
            print('{0:f} seconds.'.format(clusterbuster_end_time - clusterbuster_start_time),
                  file=log_fh)

            clusterbuster_output_fh.flush()

    # Count the number of input mutations that are associated with genes and that pass the Cluster-Buster CRM score
    # threshold.
    nbr_of_mutations_which_pass_clusterbuster_crm_score_threshold = len(
        vcf_mut_ids_passing_clusterbuster_crm_score_threshold
    )

    return nbr_of_mutations_which_pass_clusterbuster_crm_score_threshold


def calculate_and_write_motiflocator_delta_scores(vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                                                  motif_collection_version,
                                                  motif_ids_set,
                                                  motiflocator_output_filename,
                                                  min_motiflocator_score_threshold,
                                                  log_fh=sys.stderr):
    """
    Calculate the MotifLocator scores for the wildtype and mutant FASTA sequence for each mutation and for each motif
    and write the result to a file.

    :param vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict:
        OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
        mutation to the TSS and TSS as values.
        This dictionary can be made with get_all_mutations_that_overlap_with_regdoms_of_genes.
    :param motif_collection_version:
        motif collection version
    :param motif_ids_set:
        set of motif IDs
    :param motiflocator_output_filename:
        Output filename to which MotifLocator delta scores are written.
    :param min_motiflocator_score_threshold:
        Minimum MotifLocator threshold needed in wildtype or mutant to keep the result.
    :param log_fh:
        File handle to which the progress information is written.
    :return:
        nbr_of_mutations_which_pass_motiflocator_threshold: number of mutation that pass the MotifLocator threshold.
    """

    print('Score mutations with MotifLocator ...\n', file=log_fh)

    import motifsinfo
    import mutations
    import motiflocator

    # Fill MotifsInfo class with content for the chosen motif collection version.
    motifsinfo.MotifsInfo.set_motif_collection_version(motif_collection_version=motif_collection_version)

    nbr_of_mutations_which_pass_motiflocator_threshold = 0

    # Sort TADs identifiers.
    tads_ids = (
        sorted(mutations.VCFmut.tads_start_end_array_per_chrom)
        if mutations.VCFmut.tads_start_end_array_per_chrom
        else None
    )

    with tempfile.NamedTemporaryFile() as inclusive_matrix_max_motif_size_15_fh, \
            tempfile.NamedTemporaryFile() as inclusive_matrix_min_motif_size_16_max_motif_size_25_fh, \
            tempfile.NamedTemporaryFile() as inclusive_matrix_min_motif_size_26_fh, \
            open(motiflocator_output_filename, 'w') as motiflocator_output_fh:

        # Store all FilterINCLUSiveMotifsOnLength objects for the requested motifs.
        # Set different bp_upstream and bp_downstream settings (how much reference sequences to add before and after
        # the mutation position) based on the motif length so not too much excess nucleotides need to be scored which
        # won't affect the delta motif score around the mutation position.
        filtered_inclusive_motifs_on_length_dict = dict()

        # Get all motifs with a motif length smaller than or equal to 15 in a FilterINCLUSiveMotifsOnLength() object.
        filtered_inclusive_motifs_on_length_dict['max_motif_size_15'] = motifsinfo.FilterINCLUSiveMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=20,
            bp_downstream=20,
            min_motif_length=None,
            max_motif_length=15,
            header=True,
            inclusive_matrix_fh=inclusive_matrix_max_motif_size_15_fh
        )

        # Get all motifs with a motif length greater than 15 and smaller than or equal to 25  in a
        # FilterINCLUSiveMotifsOnLength() object.
        filtered_inclusive_motifs_on_length_dict['min_motif_size_16_max_motif_size_25'] = motifsinfo.FilterINCLUSiveMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=30,
            bp_downstream=30,
            min_motif_length=16,
            max_motif_length=25,
            header=True,
            inclusive_matrix_fh=inclusive_matrix_min_motif_size_16_max_motif_size_25_fh
        )

        # Get all motifs with a motif length greater than 25 in a FilterINCLUSiveMotifsOnLength() object.
        filtered_inclusive_motifs_on_length_dict['min_motif_size_26'] = motifsinfo.FilterINCLUSiveMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=60,
            bp_downstream=60,
            min_motif_length=26,
            max_motif_length=None,
            header=True,
            inclusive_matrix_fh=inclusive_matrix_min_motif_size_26_fh
        )

        # Store the keys in a separate list.
        filtered_inclusive_motifs_on_length_keys = [
            'max_motif_size_15',
            'min_motif_size_16_max_motif_size_25',
            'min_motif_size_26',
        ]

        # Loop over a copy of the list, else after deleting an element of the list,
        # the next iteration of the for loop skips an element.
        for filtered_inclusive_motifs_on_length_key in filtered_inclusive_motifs_on_length_keys[:]:
            if filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key].has_motif_ids:
                # Write temporary file with PWMs in INCLUSive format.
                filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key].write_matrix_file()
            else:
                # Remove key from dictionary if the FilterINCLUSiveMotifsOnLength() object did not retain any motif IDs.
                del filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key]

                # Also remove the entry from the filtered_inclusive_motifs_on_length_keys list.
                filtered_inclusive_motifs_on_length_keys.remove(filtered_inclusive_motifs_on_length_key)

        # Write header to the output file.
        print('# chrom',
              'start',
              'reference',
              'mutation',
              'mutation type',
              'mutation ID',
              'associated gene',
              'distance to TSS',
              'motif ID',
              'motif name',
              'directly annotated TFs',
              'wildtype MotifLocator score',
              'mutant MotifLocator score',
              'delta MotifLocator score',
              'wildtype consensus sequence',
              'mutant consensus sequence',
              'in #TADs of {0:d}'.format(len(tads_ids)) if tads_ids else '',
              *['in ' + tads_id + ' TADs' for tads_id in tads_ids] if tads_ids else '',
              sep='\t',
              file=motiflocator_output_fh)

        nbr_mutations = len(vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict)

        for mutation_idx, (vcf_mut, associated_genes_and_distance_to_tss_and_tss_dict) in enumerate(
                vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict.items()):
            print('  Scoring mutation "{0:s}" ({1:d} of {2:d}) with MotifLocator: '.format(vcf_mut.mut_id,
                                                                                           mutation_idx + 1,
                                                                                           nbr_mutations),
                  end='',
                  file=log_fh)

            motiflocators_start_time = time.time()

            motiflocator_delta_scores = dict()

            for filtered_inclusive_motifs_on_length_key in filtered_inclusive_motifs_on_length_keys:
                try:
                    fasta_string = vcf_mut.make_fasta_for_wt_and_mut(
                        bp_upstream=filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key].bp_upstream,
                        bp_downstream=filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key].bp_downstream
                    )
                except ValueError as e:
                    print('\n\n  Invalid mutation: ' + e.message + '\n',
                          file=log_fh)
                    sys.exit(1)

                # Score mutation with MotifLocator with settings provided in FilterMotifsOnLength() object.
                motiflocator_delta_scores.update(
                    motiflocator.calculate_motiflocator_delta_scores(
                        fasta_string=fasta_string,
                        inclusive_matrix_filename=filtered_inclusive_motifs_on_length_dict[filtered_inclusive_motifs_on_length_key].inclusive_matrix_fh.name,
                        min_score_threshold=min_motiflocator_score_threshold
                    )
                )

            # Count the number of input mutations that are associated with genes and that pass the MotifLocator
            # threshold.
            nbr_of_mutations_which_pass_motiflocator_threshold += (
                1 if len(motiflocator_delta_scores) != 0 else 0
            )

            # Write to the output file.
            for motif_id, motiflocator_delta in motiflocator_delta_scores.items():
                for associated_gene, distance_to_tss_and_tss in associated_genes_and_distance_to_tss_and_tss_dict.items():
                    tfs = motifsinfo.MotifsInfo.get_tfs_for_motif(motif_id)

                    distance_to_tss, tss = distance_to_tss_and_tss

                    if tads_ids:
                        tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list = [
                            1
                            if vcf_mut.tss_of_associated_gene_in_same_tad_as_mutation(
                                tads_id=tads_id,
                                tss=tss
                            )
                            else 0
                            for tads_id in tads_ids
                        ]

                    print(vcf_mut,
                          associated_gene,
                          '{0:+}'.format(distance_to_tss) if distance_to_tss else '',
                          motif_id,
                          motifsinfo.MotifsInfo.get_motif_name(motif_id),
                          ';'.join(sorted(tfs) if tfs else ['']),
                          str(motiflocator_delta.wt_score),
                          str(motiflocator_delta.mut_score),
                          str(motiflocator_delta.delta_score),
                          motiflocator_delta.wt_consensus,
                          motiflocator_delta.mut_consensus,
                          sum(tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list) if tads_ids else '',
                          *tss_of_associated_gene_in_same_tad_as_mutation_for_tads_list if tads_ids else '',
                          sep='\t',
                          file=motiflocator_output_fh)

            motiflocators_end_time = time.time()
            print('{0:f} seconds.'.format(motiflocators_end_time - motiflocators_start_time),
                  file=log_fh)

            motiflocator_output_fh.flush()

    return nbr_of_mutations_which_pass_motiflocator_threshold


def main():
    default_motif_collection_version = 'v9'
    default_min_clusterbuster_crm_score_threshold = 0.00
    default_min_motiflocator_score_threshold = 0.80

    parser = argparse.ArgumentParser(
        description='Calculate the impact of mutations on the removal or introduction of TF binding sites.'
    )

    mutations_input_group = parser.add_argument_group(
        title='Mutations',
        description='Specify mutation format of the mutations which will be scored.'
    )
    mutations_input_exclusive_group = mutations_input_group.add_mutually_exclusive_group(required=True)
    mutations_input_exclusive_group.add_argument(
        '--vcf',
        dest='vcf_filename',
        action='store',
        type=str,
        required=False,
        help='VCF file with mutations (column 1, 2, 4 and 5 are used).'
    )
    mutations_input_exclusive_group.add_argument(
        '--mut-ids',
        dest='mut_ids_filename',
        action='store',
        type=str,
        required=False,
        help='File with mutation IDs: '
             'chr10__100038800__TTTTTG__T__DEL '
             'chr10__10011659__A__AAT__INS '
             'chr10__100061062__C__T__SNV'
    )
    mutations_input_exclusive_group.add_argument(
        '--bedlike-mut-ids',
        dest='bedlike_mut_ids_filename',
        action='store',
        type=str,
        required=False,
        help='File with BED-like mutation IDs: '
             'chr10_100038800_100038801_TTTTG_-----_DEL '
             'chr10_10011659_10011660_--_AT_INS '
             'chr10_100142677_100142678_T_-_INDEL '
             'chr10_100061061_100061062_C_T_SNP'
    )

    filter_output_group = parser.add_argument_group(
        title='Filter mutations',
        description='Specify if mutations should be filtered out if they are not in a regulatory domain of a gene.'
    )
    filter_output_exclusive_group = filter_output_group.add_mutually_exclusive_group(required=False)
    filter_output_exclusive_group.add_argument(
        '--genes',
        dest='genes_filename',
        action='store',
        type=str,
        required=False,
        help='Filename with gene names to use as input. Only mutations in regulatory domains of those genes are scored.'
    )
    filter_output_exclusive_group.add_argument(
        '--keep-all-mutations',
        dest='keep_all_mutations',
        action='store_true',
        required=False,
        help='Score all mutations, even if they do not overlap with regulatory domains of genes.'
    )

    motifs_group = parser.add_argument_group(
        title='Motifs',
        description='Specify which motifs will be used for scoring the mutations.'
    )
    motifs_group.add_argument(
        '--motifs',
        dest='motif_ids_filename',
        action='store',
        type=str,
        required=False,
        help='Filename with motif IDs to score. If not specified and --tfs is also not specified, '
             'scores with all directly annotated motifs.'
    )
    motifs_group.add_argument(
        '--tfs',
        dest='tfs_filename',
        action='store',
        type=str,
        required=False,
        help='Filename with TFs for which directly annotated motif IDs will be scored. If not '
             'specified and --motifs is also not specified, scores with all directly annotated motifs.'
    )
    motifs_group.add_argument(
        '--all-motifs',
        dest='all_motifs',
        action='store_true',
        required=False,
        help='Use all motifs instead of only directly annotated motifs, if --motifs and --tfs are not specified.'
    )
    motifs_group.add_argument(
        '--motifcollection',
        dest='motif_collection_version',
        action='store',
        type=str,
        required=False,
        choices=['v7', 'v8', 'v9'],
        default=default_motif_collection_version,
        help='Motif collection to use: "v7", "v8" or "v9" (default).'
    )

    assembly_group = parser.add_argument_group(
        title='Assembly',
        description='Specify which assembly to use.'
    )
    assembly_group.add_argument(
        '--assembly',
        dest='assembly',
        action='store',
        type=str,
        required=False,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        default='hg19',
        help='Assembly version: "hg19" (TADs and associated genes: default) or "dm3", "dm6", "hg38", "mm9", "mm10" (no TADs and associated genes).'
    )

    mutation_scoring_group = parser.add_argument_group(
        title='Mutation scoring methods',
        description='Specify which mutation scoring method should be used.'
    )
    mutation_scoring_group.add_argument(
        '--clusterbuster',
        dest='clusterbuster_output_filename',
        action='store',
        type=str,
        required=False,
        help='Filename to which the Cluster-Buster CRM and motif delta score output will be written.'
    )
    mutation_scoring_group.add_argument(
        '--motiflocator',
        dest='motiflocator_output_filename',
        action='store',
        type=str,
        required=False,
        help='Filename to which the MotifLocator delta score output will be written.'
    )
    mutation_scoring_group.add_argument(
        '--min-clusterbuster-crm-score-threshold',
        dest='min_clusterbuster_crm_score_threshold',
        action='store',
        type=float,
        required=False,
        default=default_min_clusterbuster_crm_score_threshold,
        help='Minimum Cluster-Buster CRM score threshold (default: {0:f}).'.format(
            default_min_clusterbuster_crm_score_threshold)
    )
    mutation_scoring_group.add_argument(
        '--min-motiflocator-score-threshold',
        dest='min_motiflocator_score_threshold',
        action='store',
        type=float,
        required=False,
        default=default_min_motiflocator_score_threshold,
        help='Minimum MotifLocator score threshold (default: {0:f}).'.format(
            default_min_motiflocator_score_threshold)
    )

    mut2genes_group = parser.add_argument_group(
        title='Mutations to gene mapping',
        description='Save mutations to gene mapping.'
    )
    mut2genes_group.add_argument(
        '--mut2genes',
        dest='mut_to_associated_genes_output_filename',
        action='store',
        type=str,
        required=False,
        help='TSV output file with mutation info and associated gene name and distance of mutation to TSS.'
    )

    tads_group = parser.add_argument_group(
        title='TADs info',
        description='Add TADs info to output files.'
    )
    tads_group.add_argument(
        '--tads',
        dest='enable_tads_info',
        action='store_true',
        required=False,
        help='Indicate for each mutation if mutation and TSS of associated gene are located in same TAD region or not.'
    )

    logging_group = parser.add_argument_group(
        title='Logging',
        description='Specify progress logging and some statistics.'
    )
    logging_group.add_argument(
        '--log',
        dest='log_output_filename',
        action='store',
        type=str,
        required=False,
        help='Write the progress and statistics to a log file instead of standard error.'
    )

    args = parser.parse_args()

    total_start_time = time.time()

    if args.log_output_filename:
        log_fh = open(args.log_output_filename, 'w', buffering=1, encoding='utf-8')
    else:
        log_fh = sys.stderr

    print('\nCommandline:\n------------\n\n{0:s}\n\n'.format(' '.join(sys.argv)), file=log_fh)

    print('μ-cisTarget:\n------------\n', file=log_fh)
    print('Import motifsinfo ...', file=log_fh)
    import motifsinfo
    print('Import mutations ...\n', file=log_fh)
    import mutations

    # Fill MotifsInfo class with content for the chosen motif collection version.
    motifsinfo.MotifsInfo.set_motif_collection_version(motif_collection_version=args.motif_collection_version)

    # Set genomic fasta and assembly for VCFmut class.
    mutations.VCFmut.set_genome_fasta(
        fasta_filename=mutations.fasta_filename_dict[args.assembly],
        assembly=args.assembly
    )

    if args.assembly == 'hg19' or args.assembly == 'hg38' or args.assembly == 'mm9' or args.assembly == 'mm10':
        # Set regulatory domains for VCFmut class.
        mutations.VCFmut.set_reg_doms()

        if args.enable_tads_info:
            # Set TADs for VCFmut class.
            mutations.VCFmut.set_tads()

    if args.vcf_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.vcf_filename), file=log_fh)

        vcf_mut_iterator = mutations.VCFmut.from_vcf_file(args.vcf_filename)
    elif args.mut_ids_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.mut_ids_filename), file=log_fh)

        vcf_mut_iterator = mutations.VCFmut.from_mut_ids_file(args.mut_ids_filename)
    elif args.bedlike_mut_ids_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.bedlike_mut_ids_filename), file=log_fh)

        vcf_mut_iterator = mutations.VCFmut.from_bedlike_mut_ids_file(args.bedlike_mut_ids_filename)

    stats_dict = dict()

    genes_set = None

    if args.genes_filename:
        print('Read gene list from "{0:s}" ...\n'.format(args.genes_filename), file=log_fh)
        genes_set = read_genes_filename(args.genes_filename)

        if len(genes_set) == 0:
            print('An empty gene list file was provided.\n', file=log_fh)
            sys.exit(1)

    motif_ids_set = set()

    if args.motif_ids_filename:
        print('Read motif IDs from "{0:s}" ...\n'.format(args.motif_ids_filename), file=log_fh)

        motif_ids_set = read_motif_ids_filename(args.motif_ids_filename)

        if len(motif_ids_set) == 0:
            print('An empty motif IDs file was provided.\n', file=log_fh)
            sys.exit(1)

        all_motif_ids_set = motifsinfo.MotifsInfo.get_all_motif_ids(directly_annotated_motifs_only=False)

        if not motif_ids_set.issubset(all_motif_ids_set):
            print('The following motif IDs are not recognised:\n  - ',
                  '\n  - '.join(motif_ids_set.difference(all_motif_ids_set)),
                  '\n',
                  sep='',
                  file=log_fh)
            sys.exit(1)

    if args.tfs_filename:
        print('Read TFs from "{0:s}" ...\n'.format(args.tfs_filename), file=log_fh)
        tfs_set = read_tfs_filename(args.tfs_filename)

        tfs_with_directly_annotated_motifs_set = {tf for tf in tfs_set if tf in motifsinfo.MotifsInfo.tf_to_motifs_dict}

        if len(tfs_with_directly_annotated_motifs_set) == 0:
            print('None of the provided TFs have directly annotated motifs.\n', file=log_fh)
            sys.exit(1)

        stats_dict['nbr_of_input_tfs'] = len(tfs_set)
        stats_dict['nbr_of_input_tfs_with_directly_annotated_motifs'] = len(tfs_with_directly_annotated_motifs_set)

        # Add all motif IDs which are directly annotated for a TF.
        [motif_ids_set.update(motifsinfo.MotifsInfo.get_motifs_for_tf(tf))
         for tf in tfs_with_directly_annotated_motifs_set]

    if not args.motif_ids_filename and not args.tfs_filename:
        # Use all directly annotated motif IDs or all motif IDs (if all_motifs is True) if --motifs and --tfs are not specified.
        motif_ids_set = motifsinfo.MotifsInfo.get_all_motif_ids(directly_annotated_motifs_only=not args.all_motifs)

    print('Get all mutations that overlap with the regulatory domains of {0:s}: '.format('the provided gene set'
                                                                                         if args.genes_filename
                                                                                         else 'genes'),
          end='',
          file=log_fh)

    mut_to_associated_genes_start_time = time.time()

    try:
        (vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
         input_vcf_mut_ids,
         associated_genes_set
         ) = get_all_mutations_that_overlap_with_regdoms_of_genes(vcf_mut_iterator, genes_set, args.keep_all_mutations)
    except ValueError as e:
        print('\nInvalid mutation: ' + e.message + '\n',
              file=log_fh)
        sys.exit(1)

    stats_dict['nbr_of_input_mutations'] = len(input_vcf_mut_ids)
    stats_dict['nbr_of_mutations_associated_with_genes'] = len(
        list(filter(None, vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict.values()))
    )
    stats_dict['nbr_of_genes_associated_with_mutations'] = len(associated_genes_set)

    mut_to_associated_genes_end_time = time.time()

    print('{0:f} seconds.\n'.format(mut_to_associated_genes_end_time - mut_to_associated_genes_start_time),
          file=log_fh)

    if args.mut_to_associated_genes_output_filename:
        # Write all mutations and their associated genes and distance of the mutation to the TSS to a file.
        write_mut_to_associated_gene_output(
            mut_to_associated_genes_output_filename=args.mut_to_associated_genes_output_filename,
            vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict=vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
            log_fh=log_fh
        )

    if args.clusterbuster_output_filename:
        # Calculate the Cluster-Buster CRM and motif scores for the wildtype and mutant FASTA sequence for each
        # mutation and for each motif and write the result to a file.
        stats_dict['nbr_of_mutations_which_pass_clusterbuster_crm_score_threshold'] = \
            calculate_and_write_clusterbuster_crm_and_motif_delta_scores(
                vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict=vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                motif_collection_version=args.motif_collection_version,
                motif_ids_set=motif_ids_set,
                clusterbuster_output_filename=args.clusterbuster_output_filename,
                min_clusterbuster_crm_score_threshold=args.min_clusterbuster_crm_score_threshold,
                log_fh=log_fh
        )

    if args.motiflocator_output_filename:
        # Calculate the MotifLocator scores for the wildtype and mutant FASTA sequence for each mutation and for each
        # motif and write the result to a file.
        stats_dict['nbr_of_mutations_which_pass_motiflocator_threshold'] = \
            calculate_and_write_motiflocator_delta_scores(
                vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict=vcf_mut_to_associated_genes_and_distance_to_tss_and_tss_dict,
                motif_collection_version=args.motif_collection_version,
                motif_ids_set=motif_ids_set,
                motiflocator_output_filename=args.motiflocator_output_filename,
                min_motiflocator_score_threshold=args.min_motiflocator_score_threshold,
                log_fh=log_fh
        )

    total_end_time = time.time()

    print('\n\nTotal time: {0:f} seconds.\n'.format(total_end_time - total_start_time),
          file=log_fh)

    # Print some statistics.
    print(
        '\nStatistics:',
        '-----------\n',
        'Number of mutations in input file:\t{0:d}'.format(
            stats_dict['nbr_of_input_mutations']),
        'Number of mutations associated with genes:\t{0:d}'.format(
            stats_dict['nbr_of_mutations_associated_with_genes']),
        'Number of genes associated with mutations:\t{0:d}'.format(
            stats_dict['nbr_of_genes_associated_with_mutations']),
        sep='\n',
        file=log_fh
    )

    if args.clusterbuster_output_filename:
        print(
            'Number of mutations which pass Cluster-Buster CRM score threshold:\t{0:d}'.format(
                stats_dict['nbr_of_mutations_which_pass_clusterbuster_crm_score_threshold']),
            sep='\n',
            file=log_fh
        )

    if args.motiflocator_output_filename:
        print(
            'Number of mutations which pass MotifLocator threshold:\t{0:d}'.format(
                stats_dict['nbr_of_mutations_which_pass_motiflocator_threshold']),
            sep='\n',
            file=log_fh
        )

    print(
        'Number of motifs used for scoring: {0:d}'.format(len(motif_ids_set)),
        sep='\n',
        file=log_fh
    )

    if args.tfs_filename:
        print(
            'Number of input TFs:\t{0:d}'.format(stats_dict['nbr_of_input_tfs']),
            'Number of input TFs with annotated motifs:\t{0:d}'.format(
                stats_dict['nbr_of_input_tfs_with_directly_annotated_motifs']),
            sep='\n',
            file=log_fh
        )

    print('', file=log_fh)

    log_fh.close()


if __name__ == "__main__":
    main()
