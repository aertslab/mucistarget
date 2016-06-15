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

        # Write temporary file with PWMs in INCLUSive format which length is smaller than or equal to 15.
        matrix_max_motif_size_15_fh.write(
            motifsinfo.MotifsInfo.get_pwms(
                motif_ids=motif_ids_set,
                min_motif_length=None,
                max_motif_length=15,
                header=True
            )
        )
        matrix_max_motif_size_15_fh.flush()

        # Write temporary file with PWMs in INCLUSive format which length is greater than 15 and smaller than or equal to 25.
        matrix_min_motif_size_16_max_motif_size_25_fh.write(
            motifsinfo.MotifsInfo.get_pwms(
                motif_ids=motif_ids_set,
                min_motif_length=16,
                max_motif_length=25,
                header=True
            )
        )
        matrix_min_motif_size_16_max_motif_size_25_fh.flush()

        # Write temporary file with PWMs in INCLUSive format which length is greater than 25.
        matrix_min_motif_size_26_fh.write(
            motifsinfo.MotifsInfo.get_pwms(
                motif_ids=motif_ids_set,
                min_motif_length=26,
                max_motif_length=None,
                header=True
            )
        )
        matrix_min_motif_size_26_fh.flush()

        print('# chrom',
              'start',
              'reference',
              'mutation',
              'mutation ID',
              'associated genes',
              'motif ID',
              'motif name',
              'directly annotated TFs',
              'wildtype MotifLocator score',
              'mutant MotifLocator score',
              'delta MotifLocator score',
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

            associated_genes_set = vcf_mut.get_associated_genes()

            if not associated_genes_set.isdisjoint(genes_set):
                # Count the number of input mutations that are associated with genes.
                mutations_stats['nbr_of_mutations_associated_with_genes'] += 1

                # Only consider mutations which have associated genes which appear in our input set.
                print('Scoring mutation "{0:s}" with MotifLocator: '.format(vcf_mut.mut_id),
                      end='',
                      file=sys.stderr)

                motiflocators_start_time = time.time()

                # Score mutation with MotifLocator with motifs which length is smaller than or equal to 15.
                motiflocator_delta_scores_max_motif_size_15 = motiflocator.calculate_motiflocator_delta_scores(
                    fasta_string=vcf_mut.make_fasta_for_wt_and_mut(bp_upstream=20, bp_downstream=20),
                    matrix_filename=matrix_max_motif_size_15_fh.name
                )

                # Score mutation with MotifLocator with motifs which length is greater than 15 and smaller than or equal to 25.
                motiflocator_delta_scores_min_motif_size_16_max_motif_size_25 = motiflocator.calculate_motiflocator_delta_scores(
                    fasta_string=vcf_mut.make_fasta_for_wt_and_mut(bp_upstream=30, bp_downstream=30),
                    matrix_filename=matrix_min_motif_size_16_max_motif_size_25_fh.name
                )

                # Score mutation with MotifLocator with motifs which length is greater than 25.
                motiflocator_delta_scores_min_motif_size_26 = motiflocator.calculate_motiflocator_delta_scores(
                    fasta_string=vcf_mut.make_fasta_for_wt_and_mut(bp_upstream=60, bp_downstream=60),
                    matrix_filename=matrix_min_motif_size_26_fh.name
                )

                # Combine all individual dictionaries with MotifLocator delta scores.
                motiflocator_delta_scores = dict()
                motiflocator_delta_scores.update(motiflocator_delta_scores_max_motif_size_15)
                motiflocator_delta_scores.update(motiflocator_delta_scores_min_motif_size_16_max_motif_size_25)
                motiflocator_delta_scores.update(motiflocator_delta_scores_min_motif_size_26)

                # Count the number of input mutations that are associated with genes and that pass the MotifLocator
                # threshold.
                mutations_stats['nbr_of_mutations_which_pass_motiflocator_threshold'] += (
                    1 if len(motiflocator_delta_scores) != 0 else 0
                )

                # Write to the output file.
                for motif_id, motiflocator_delta in motiflocator_delta_scores.iteritems():
                    print(vcf_mut,
                          ';'.join(associated_genes_set),
                          '\t'.join([motif_id,
                                     motifsinfo.MotifsInfo.get_motif_name(motif_id),
                                     ';'.join(motifsinfo.MotifsInfo.get_tfs_for_motif(motif_id)),
                                     str(motiflocator_delta.wt_score),
                                     str(motiflocator_delta.mut_score),
                                     str(motiflocator_delta.delta_score)
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
