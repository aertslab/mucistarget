#!/usr/bin/env bash
#
# Purpose :      Split big VCF files in multiple smaller VCF files with
#                (almost) the same number of mutations in each file.
#
# Copyright (C): 2019 - Gert Hulselmans



split_vcf_file_in_equal_parts () {
    if [ "${#@}" -ne 4 ] ; then
        printf 'Usage: %s vcf_input_filename vcf_output_filename_prefix vcf_output_filename_postfix nbr_splits\n' "${0}";
        exit 1;
    fi

    local vcf_input_filename="${1}";
    local vcf_output_filename_prefix="${2}";
    local vcf_output_filename_postfix="${3}";
    local -i nbr_splits="${4}";

    if [ ${nbr_splits} -eq 0 ] ; then
        printf 'Error: Number of splits was not a valid number or zero.\n';
        exit 1;
    fi

    awk -F '\t' \
        -v vcf_output_filename_prefix="${vcf_output_filename_prefix}" \
        -v vcf_output_filename_postfix="${vcf_output_filename_postfix}" \
        -v nbr_splits="${nbr_splits}" \
        '
        {
            if ($0 ~ /^#/) {
                # Write headers and comments to all VCF output files.
                for (i = 0; i < nbr_splits; i++) {
                    vcf_output_filename = sprintf("%s%03d%s", vcf_output_filename_prefix, i, vcf_output_filename_postfix);
                    print $0 > vcf_output_filename;
                }
            } else {
                # Write a line to a specific VCF output file depending on the line number.
                vcf_output_filename = sprintf("%s%03d%s", vcf_output_filename_prefix, NR % nbr_splits, vcf_output_filename_postfix);
                print $0 > vcf_output_filename;
            }
        } END {
            # Close all VCF output file handles.
            for (i = 0; i < nbr_splits; i++) {
                vcf_output_filename = sprintf("%s%03d%s", vcf_output_filename_prefix, i, vcf_output_filename_postfix);
                close(vcf_output_filename);
            }
        }
        ' "${vcf_input_filename}";
}



split_vcf_file_in_equal_parts "${@}";