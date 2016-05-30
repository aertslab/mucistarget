#!/bin/bash
#
# Purpose:       Convert Cluster Buster motif files to INCLUSive (MotifLocator) motif format.
#
# Copyright (C): 2015-2016 - Gert Hulselmans


usage () {
    printf '\nUsage:   %s [-h|--help] <motif1.cb> <motif2.cb> ... <motifn.cb> > motifs.INCLUSive.txt\n\n' "${0}";
    printf 'Purpose: Convert Cluster Buster motif files to INCLUSive (MotifLocator) motif format.\n\n';
    
    exit 1;
}


if [ ${#@} -eq 0 ] ; then
    usage;
elif [ "${1}" = "-h" ] ; then
    usage;
elif [ "${1}" = "--help" ] ; then
    usage;
fi

for cluster_buster_file in "${@}" ; do
    if [ ! -f "${cluster_buster_file}" ] ; then
        printf '\nERROR: "%s" could not be found or is not a file.\n\n' "${cluster_buster_file}";
        exit 1;
    fi
done


# Convert all Cluster Buster matrices to INCLUSive (MotifLocator) motif format.
awk -F "\t" '
    BEGIN {
        # Print INCLUSive header.
        print "#INCLUSive Motif Model\n";
    }
    {
        if ($0 ~ /^>/) {
            # Do not print the motif when this is the first motif name.
            if (motif_name != "") {
                # Print the motif in INCLUSive (MotifLocator) format when it is
                # not the first time a line starts with ">".
                print "#ID = " motif_name "\n#W = " motif_length motif_matrix "\n";
            };
            
            # Extract motif name from the Cluster Buster ">motif_name" line.
            motif_name = substr($0, 2);
                    
            # Set motif length back to zero for the next motif.
            motif_length = 0;
            
            # Set motif matrix back to an empty string.
            motif_matrix = "";
        } else {
            # Skip empty lines.
            if ($0 != "") {
                # Read Cluster Buster matrix lines.
                
                # Add 1 to the current motif length.
                motif_length += 1;
                
                # Add motif matrix line to already saved motif matrix lines.
                motif_matrix = motif_matrix "\n" $0;
            }
        }
    }
    END {
        # Do not print the motif when this is the first motif name.
            if (motif_name != "") {
                # Print the last motif in INCLUSive (MotifLocator) format.
                print "#ID = " motif_name "\n#W = " motif_length motif_matrix "\n";
            };
    }' \
    "${@}";

exit $?;
