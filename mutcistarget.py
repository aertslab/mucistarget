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
                        help='Filename with motif IDs to score. If not specified, scores with all motifs.'
                        )
    parser.add_argument('--output',
                        dest='output_filename',
                        action='store',
                        type=str,
                        required=True,
                        help='Filename to which the output will be written.'
                        )

    args = parser.parse_args()

    print('\nImport motiflocator ...', file=sys.stderr)
    import motiflocator
    print('Import motifsinfo ...', file=sys.stderr)
    import motifsinfo
    print('Import mutations ...\n', file=sys.stderr)
    import mutations

    if args.vcf_filename:
        print('Using mutation filename: "{0:s}"'.format(args.vcf_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_vcf_file(args.vcf_filename)
    elif args.mut_ids_filename:
        print('Using mutation filename: "{0:s}"'.format(args.mut_ids_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_mut_ids_file(args.mut_ids_filename)
    elif args.bedlike_mut_ids_filename:
        print('Using mutation filename: "{0:s}"'.format(args.bedlike_mut_ids_filename), file=sys.stderr)

        vcf_mut_iterator = mutations.VCFmut.from_bedlike_mut_ids_file(args.bedlike_mut_ids_filename)

    print('Read gene list from "{0:s}" ...\n'.format(args.genes_filename), file=sys.stderr)
    genes_set = read_genes_filename(args.genes_filename)

    if args.motif_ids_filename:
        print('Read motif IDs from "{0:s}" ...\n'.format(args.motif_ids_filename), file=sys.stderr)
        motif_ids_set = read_motif_ids_filename(args.motif_ids_filename)
    else:
        motif_ids_set = None

    print('Score mutations with MotifLocator ...\n', file=sys.stderr)

    with tempfile.NamedTemporaryFile() as matrix_max_motif_size_15_fh, \
            tempfile.NamedTemporaryFile() as matrix_min_motif_size_16_max_motif_size_25_fh, \
            tempfile.NamedTemporaryFile() as matrix_min_motif_size_26_fh, \
            open(args.output_filename, 'w') as output_fh:

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
              file=output_fh)

        mutations_stats = {
            'nbr_of_input_mutations': 0,
            'nbr_of_mutations_associated_with_genes': 0,
            'nbr_of_mutations_which_pass_motiflocator_threshold': 0,
        }

        for vcf_mut in vcf_mut_iterator:
            # Count the number of input mutations.
            mutations_stats['nbr_of_input_mutations'] += 1

            associated_genes_and_distance_to_tss_dict = vcf_mut.get_associated_genes_and_distance_to_tss()

            if not set(associated_genes_and_distance_to_tss_dict).isdisjoint(genes_set):
                # Count the number of input mutations that are associated with genes.
                mutations_stats['nbr_of_mutations_associated_with_genes'] += 1

                # Only consider mutations which have associated genes which appear in our input set.
                print('Scoring mutation "{0:s}" with MotifLocator: '.format(vcf_mut.mut_id),
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
                              file=output_fh)

                motiflocators_end_time = time.time()
                print('{0:f} seconds.'.format(motiflocators_end_time - motiflocators_start_time),
                      file=sys.stderr)

                output_fh.flush()

        # Print some statistics about the number of mutations.
        print(
            '\nNumber of mutations in input file: {0:d}'.format(
                mutations_stats['nbr_of_input_mutations']),
            'Number of mutations associated with genes: {0:d}'.format(
                mutations_stats['nbr_of_mutations_associated_with_genes']),
            'Number of mutations which pass MotifLocator threshold: {0:d}'.format(
                mutations_stats['nbr_of_mutations_which_pass_motiflocator_threshold']),
            sep='\n',
            file=sys.stderr
        )


if __name__ == "__main__":
    main()
