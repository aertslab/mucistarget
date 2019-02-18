#!/usr/bin/env python3

"""
Purpose :      Create VCF files (first 5 columns in VCF format) from various input formats.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Create VCF files (first 5 columns in VCF format) from various input formats.'
    )

    input_group = parser.add_argument_group(
        title='Input file types',
        description='Supported input file types.'
    )
    input_mutually_exclusive_group = input_group.add_mutually_exclusive_group(required=True)
    input_mutually_exclusive_group.add_argument('--mut-ids',
                       dest='mut_ids_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='File with mutation IDs: '
                            'chr10__100038800__TTTTTG__T__DEL '
                            'chr10__10011659__A__AAT__INS '
                            'chr10__100061062__C__T__SNV')
    input_mutually_exclusive_group.add_argument('--bedlike-mut-ids',
                       dest='bedlike_mut_ids_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='File with BED-like mutation IDs: '
                            'chr10_100038800_100038801_TTTTG_-----_DEL '
                            'chr10_10011659_10011660_--_AT_INS '
                            'chr10_100142677_100142678_T_-_INDEL '
                            'chr10_100061061_100061062_C_T_SNP '
                            'chr10_100038800_100038801_TTTTG_----- '
                            'chr10_100038800_100038801_TTTTG_- '
                            'chr10_10011659_10011660_--_AT '
                            'chr10_10011659_10011660_-_AT '
                            'chr10_100142677_100142678_T_- '
                            'chr10_100061061_100061062_C_T')
    input_mutually_exclusive_group.add_argument('--complete-genomics',
                       dest='complete_genomics_filename',
                       action='store',
                       type=str,
                       required=False,
                       help='File with Complete Genomics mutation calls.')

    parser.add_argument('--column-number',
                       dest='column_number',
                       action='store',
                       type=int,
                       required=False,
                       default=1,
                       help='Column number which contains the (BED-like) mutation ID (default: 1). ')

    output_group = parser.add_argument_group(
        title='Output file types',
        description='Supported output file types.'
    )
    output_group.add_argument('--vcf',
                       dest='vcf_filename',
                       action='store',
                       type=str,
                       required=True,
                       help='VCF output file with mutations (first 5 columns are written in VCF format, '
                            'rest of the line will be the content of the original line).')

    args = parser.parse_args()

    import mutations

    if args.complete_genomics_filename:
        with open(args.vcf_filename, 'w') as vcf_fh, open(args.complete_genomics_filename, 'r') as gc_fh:
            for line in gc_fh:
                line = line.rstrip('\r\n')

                if line == '' or line.startswith('#') or line.startswith('>'):
                    continue

                columns = line.split('\t')

                if len(columns) >= 9:
                    chrom = columns[3]
                    start = columns[4]
                    mut_type = columns[6]
                    ref = columns[7]
                    mut = columns[8]

                    if mut_type in ('snp', 'del', 'ins'):
                        # Only process SNPs, deletions and insertions.
                        vcf_mut = mutations.VCFmut.from_zero_based_no_ref_specified(chrom, start, ref, mut)

                        print(vcf_mut.chrom, vcf_mut.start, vcf_mut.mut_id, vcf_mut.ref, vcf_mut.mut, line,
                              sep='\t',
                              file=vcf_fh)
    else:
        if args.mut_ids_filename:
            vcf_mut_iterator = mutations.VCFmut.from_mut_ids_file(
                mut_ids_filename=args.mut_ids_filename,
                mut_ids_column=args.column_number,
                return_input_line=True)
        elif args.bedlike_mut_ids_filename:
            vcf_mut_iterator = mutations.VCFmut.from_bedlike_mut_ids_file(
                bedlike_mut_ids_filename=args.bedlike_mut_ids_filename,
                bedlike_mut_ids_column=args.column_number,
                return_input_line=True
            )

        with open(args.vcf_filename, 'w') as vcf_fh:
            for vcf_mut, line in vcf_mut_iterator:
                print(vcf_mut.chrom, vcf_mut.start, vcf_mut.mut_id, vcf_mut.ref, vcf_mut.mut, line,
                      sep='\t',
                      file=vcf_fh)


if __name__ == "__main__":
    main()
