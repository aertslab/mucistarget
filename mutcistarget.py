#!/usr/bin/env python

"""
Purpose :      Calculate the impact of mutations on the removal or introduction of TF binding sites.

Copyright (C): 2016 - Gert Hulselmans
"""

from __future__ import print_function

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


def get_all_mutations_that_overlap_with_regdoms_of_genes(vcf_mut_iterator, genes_set):
    """
    Get all unique mutations that fall in a regulatory domain of one of the genes in the gene set.

    :param vcf_mut_iterator:
        Iterator that yields a VCFmut object:
            - mutations.VCFmut.from_vcf_file()
            - mutations.VCFmut.from_mut_ids_file()
            - mutations.VCFmut.from_bedlike_mut_ids_file()
    :param genes_set:
        Gene set which regulatory domains will be used to only keep those mutations that fall in those domains.
    :return:
        (vcf_mut_to_associated_genes_and_distance_to_tss_dict,
         input_vcf_mut_ids)

        Where vcf_mut_to_associated_genes_and_distance_to_tss_dict:
            OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
            mutation to the TSS as values.
        Where input_vcf_mut_ids:
            unique mutation IDs in the vcf_mut_iterator.
    """

    input_vcf_mut_ids = set()
    vcf_mut_to_associated_genes_and_distance_to_tss_dict = OrderedDict()

    for vcf_mut in vcf_mut_iterator:
        if vcf_mut.mut_id not in input_vcf_mut_ids:
            # Store all unique mutation IDs.
            input_vcf_mut_ids.add(vcf_mut.mut_id)

            associated_genes_and_distance_to_tss_dict = vcf_mut.get_associated_genes_and_distance_to_tss()

            if not set(associated_genes_and_distance_to_tss_dict).isdisjoint(genes_set):
                # Store all mutations (VCFmut object) and associated genes information in a ordered dict.
                vcf_mut_to_associated_genes_and_distance_to_tss_dict[vcf_mut] \
                    = associated_genes_and_distance_to_tss_dict

    return (vcf_mut_to_associated_genes_and_distance_to_tss_dict,
            input_vcf_mut_ids)


def write_mut_to_associated_gene_output(mut_to_associated_genes_output_filename,
                                        vcf_mut_to_associated_genes_and_distance_to_tss_dict):
    """
    Write all mutations and their associated genes and distance of the mutation to the TSS to a file.

    :param mut_to_associated_genes_output_filename:
        Output filename which will contain mutation info and associated genes info.
    :param vcf_mut_to_associated_genes_and_distance_to_tss_dict:
        OrderedDict with VCFmut objects as keys and as values a dictionary of gene names as keys and distance of the
        mutation to the TSS as values.
        This dictionary can be made with get_all_mutations_that_overlap_with_regdoms_of_genes.
    :return:
    """

    with open(mut_to_associated_genes_output_filename, 'w') as mut_to_associated_genes_fh:
        # Write header to the output file.
        print('# chrom',
              'start',
              'reference',
              'mutation',
              'mutation type',
              'mutation ID',
              'associated gene',
              'distance to TSS',
              sep='\t',
              file=mut_to_associated_genes_fh)

        for vcf_mut, associated_genes_and_distance_to_tss_dict in vcf_mut_to_associated_genes_and_distance_to_tss_dict.iteritems():
            # Write to the output file.
            for associated_gene, distance_to_tss in associated_genes_and_distance_to_tss_dict.iteritems():
                print(vcf_mut,
                      associated_gene,
                      '{0:+}'.format(distance_to_tss),
                      sep='\t',
                      file=mut_to_associated_genes_fh)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate the impact of mutations on the removal or introduction of TF binding sites.'
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--vcf',
                       dest='vcf_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='VCF file with mutations (column 1, 2, 3 and 5 are used).')
    group.add_argument('--mut-ids',
                       dest='mut_ids_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='File with mutation IDs: '
                            'chr10__100038800__TTTTTG__T__DEL '
                            'chr10__10011659__A__AAT__INS '
                            'chr10__100061062__C__T__SNV')
    group.add_argument('--bedlike-mut-ids',
                       dest='bedlike_mut_ids_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='File with BED-like mutation IDs: '
                            'chr10_100038800_100038801_TTTTG_-----_DEL '
                            'chr10_10011659_10011660_--_AT_INS '
                            'chr10_100142677_100142678_T_-_INDEL '
                            'chr10_100061061_100061062_C_T_SNP')
    parser.add_argument('--genes',
                        dest='genes_filename',
                        action='store',
                        type=str,
                        required=True,
                        help='Filename with gene names to use as input.'
                        )
    parser.add_argument('--motifs',
                        dest='motif_ids_filename',
                        action='store',
                        type=str,
                        required=False,
                        help='Filename with motif IDs to score. If not specified and --tfs is also not specified, '
                             'scores with all motifs.'
                        )
    parser.add_argument('--tfs',
                        dest='tfs_filename',
                        action='store',
                        type=str,
                        required=False,
                        help='Filename with TFs for which directly annotated motif IDs will be scored. If not '
                             'specified and --motifs is also not specified, scores with all motifs.'
                        )
    parser.add_argument('--motiflocator',
                        dest='motiflocator_output_filename',
                        action='store',
                        type=str,
                        required=True,
                        help='Filename to which the MotifLocator delta score output will be written.'
                        )
    parser.add_argument('--mut2genes',
                        dest='mut_to_associated_genes_output_filename',
                        action='store',
                        type=str,
                        required=False,
                        help='TSV output file with mutation info and associated gene name and distance of mutation to '
                             'TSS.'
                        )

    args = parser.parse_args()

    print('\nImport motiflocator ...', file=sys.stderr)
    import motiflocator
    print('Import motifsinfo ...', file=sys.stderr)
    import motifsinfo
    print('Import mutations ...\n', file=sys.stderr)
    import mutations

    if args.vcf_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.vcf_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_vcf_file(args.vcf_filename)
    elif args.mut_ids_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.mut_ids_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_mut_ids_file(args.mut_ids_filename)
    elif args.bedlike_mut_ids_filename:
        print('Using mutation filename: "{0:s}"\n'.format(args.bedlike_mut_ids_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_bedlike_mut_ids_file(args.bedlike_mut_ids_filename)

    print('Read gene list from "{0:s}" ...\n'.format(args.genes_filename), file=sys.stderr)
    genes_set = read_genes_filename(args.genes_filename)

    motif_ids_set = set()

    if args.motif_ids_filename:
        print('Read motif IDs from "{0:s}" ...\n'.format(args.motif_ids_filename), file=sys.stderr)
        motif_ids_set.update(read_motif_ids_filename(args.motif_ids_filename))
    if args.tfs_filename:
        print('Read TFs from "{0:s}" ...\n'.format(args.tfs_filename), file=sys.stderr)
        tfs_set = read_tfs_filename(args.tfs_filename)

        # Add all motif IDs which are directly annotated for a TF.
        [motif_ids_set.update(motifsinfo.MotifsInfo.get_motifs_for_tf(tf)) for tf in tfs_set]

    if not args.motif_ids_filename and not args.tfs_filename:
        # Use all motif IDs if --motifs and --tfs are not specified.
        motif_ids_set = set(motifsinfo.MotifsInfo.motif_id_to_filename_dict)

    mutations_stats = dict()

    print('Get all mutations that overlap with the regulatory domains of the provided gene set: ',
          end='',
          file=sys.stderr)

    mut_to_associated_genes_start_time = time.time()

    (vcf_mut_to_associated_genes_and_distance_to_tss_dict,
     input_vcf_mut_ids
     ) = get_all_mutations_that_overlap_with_regdoms_of_genes(vcf_mut_iterator, genes_set)

    mutations_stats['nbr_of_input_mutations'] = len(input_vcf_mut_ids)
    mutations_stats['nbr_of_mutations_associated_with_genes'] = len(vcf_mut_to_associated_genes_and_distance_to_tss_dict)

    mut_to_associated_genes_end_time = time.time()

    print('{0:f} seconds.\n'.format(mut_to_associated_genes_end_time - mut_to_associated_genes_start_time),
          file=sys.stderr)

    if args.mut_to_associated_genes_output_filename:
        # Write all mutations and their associated genes and distance of the mutation to the TSS to a file.
        write_mut_to_associated_gene_output(
            mut_to_associated_genes_output_filename=args.mut_to_associated_genes_output_filename,
            vcf_mut_to_associated_genes_and_distance_to_tss_dict=vcf_mut_to_associated_genes_and_distance_to_tss_dict
        )

    print('Score mutations with MotifLocator ...\n', file=sys.stderr)

    mutations_stats['nbr_of_mutations_which_pass_motiflocator_threshold'] = 0

    with tempfile.NamedTemporaryFile() as matrix_max_motif_size_15_fh, \
            tempfile.NamedTemporaryFile() as matrix_min_motif_size_16_max_motif_size_25_fh, \
            tempfile.NamedTemporaryFile() as matrix_min_motif_size_26_fh, \
            open(args.motiflocator_output_filename, 'w') as motiflocator_output_fh:

        filtered_motifs_on_length_dict = dict()

        # Get all motifs with a motif length smaller than or equal to 15 in a FilterMotifsOnLength() object.
        filtered_motifs_on_length_dict['max_motif_size_15'] = motifsinfo.FilterMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=20,
            bp_downstream=20,
            min_motif_length=None,
            max_motif_length=15,
            header=True,
            matrix_fh=matrix_max_motif_size_15_fh
        )

        # Get all motifs with a motif length greater than 15 and smaller than or equal to 25  in a
        # FilterMotifsOnLength() object.
        filtered_motifs_on_length_dict['min_motif_size_16_max_motif_size_25'] = motifsinfo.FilterMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=30,
            bp_downstream=30,
            min_motif_length=16,
            max_motif_length=25,
            header=True,
            matrix_fh=matrix_min_motif_size_16_max_motif_size_25_fh
        )

        # Get all motifs with a motif length greater than 25 in a FilterMotifsOnLength() object.
        filtered_motifs_on_length_dict['min_motif_size_26'] = motifsinfo.FilterMotifsOnLength(
            motif_ids=motif_ids_set,
            bp_upstream=60,
            bp_downstream=60,
            min_motif_length=26,
            max_motif_length=None,
            header=True,
            matrix_fh=matrix_min_motif_size_26_fh
        )

        # Store the keys in a separate list.
        filtered_motifs_on_length_keys = [
            'max_motif_size_15',
            'min_motif_size_16_max_motif_size_25',
            'min_motif_size_26',
        ]

        # Loop over the list.
        for filtered_motifs_on_length_key in filtered_motifs_on_length_keys:
            if filtered_motifs_on_length_dict[filtered_motifs_on_length_key].has_motif_ids:
                # Write temporary file with PWMs in INCLUSive format.
                filtered_motifs_on_length_dict[filtered_motifs_on_length_key].write_matrix_file()
            else:
                # Remove key from dictionary if the FilterMotifsOnLength() object did not retain any motif IDs.
                del filtered_motifs_on_length_dict[filtered_motifs_on_length_key]

                # Also remove the entry from the filtered_motifs_on_length_keys list.
                filtered_motifs_on_length_keys.remove(filtered_motifs_on_length_key)

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
              sep='\t',
              file=motiflocator_output_fh)

        for vcf_mut in vcf_mut_iterator:
            associated_genes_and_distance_to_tss_dict = vcf_mut.get_associated_genes_and_distance_to_tss()

            if not set(associated_genes_and_distance_to_tss_dict).isdisjoint(genes_set):
                # Only consider mutations which have associated genes which appear in our input set.
                print('  Scoring mutation "{0:s}" with MotifLocator: '.format(vcf_mut.mut_id),
                      end='',
                      file=sys.stderr)

                motiflocators_start_time = time.time()

                motiflocator_delta_scores = dict()

                for filtered_motifs_on_length_key in filtered_motifs_on_length_keys:
                    # Score mutation with MotifLocator with settings provided in FilterMotifsOnLength() object.
                    motiflocator_delta_scores.update(
                        motiflocator.calculate_motiflocator_delta_scores(
                            fasta_string=vcf_mut.make_fasta_for_wt_and_mut(
                                bp_upstream=filtered_motifs_on_length_dict[filtered_motifs_on_length_key].bp_upstream,
                                bp_downstream=filtered_motifs_on_length_dict[filtered_motifs_on_length_key].bp_downstream
                            ),
                            matrix_filename=filtered_motifs_on_length_dict[filtered_motifs_on_length_key].matrix_fh.name
                        )
                    )

                # Count the number of input mutations that are associated with genes and that pass the MotifLocator
                # threshold.
                mutations_stats['nbr_of_mutations_which_pass_motiflocator_threshold'] += (
                    1 if len(motiflocator_delta_scores) != 0 else 0
                )

                # Write to the output file.
                for motif_id, motiflocator_delta in motiflocator_delta_scores.iteritems():
                    for associated_gene, distance_to_tss in associated_genes_and_distance_to_tss_dict.iteritems():
                        print(vcf_mut,
                              associated_gene,
                              '{0:+}'.format(distance_to_tss),
                              '\t'.join([motif_id,
                                         motifsinfo.MotifsInfo.get_motif_name(motif_id),
                                         ';'.join(motifsinfo.MotifsInfo.get_tfs_for_motif(motif_id)),
                                         str(motiflocator_delta.wt_score),
                                         str(motiflocator_delta.mut_score),
                                         str(motiflocator_delta.delta_score),
                                         motiflocator_delta.wt_consensus,
                                         motiflocator_delta.mut_consensus,
                                         ]
                                        ),
                              sep='\t',
                              file=motiflocator_output_fh)

                motiflocators_end_time = time.time()
                print('{0:f} seconds.'.format(motiflocators_end_time - motiflocators_start_time),
                      file=sys.stderr)

                motiflocator_output_fh.flush()

    # Print some statistics about the number of mutations.
    print(
        '\nNumber of mutations in input file: {0:d}'.format(
            mutations_stats['nbr_of_input_mutations']),
        'Number of mutations associated with genes: {0:d}'.format(
            mutations_stats['nbr_of_mutations_associated_with_genes']),
        'Number of mutations which pass MotifLocator threshold: {0:d}'.format(
            mutations_stats['nbr_of_mutations_which_pass_motiflocator_threshold']),
        'Number of motifs used for scoring: {0:d}\n'.format(len(motif_ids_set)),
        sep='\n',
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
