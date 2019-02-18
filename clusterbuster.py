"""
Purpose :      Run Cluster-Buster and get Cluster-Buster CRM and motif delta scores.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

from __future__ import print_function

from io import BytesIO

import command
import mutations
import motifsinfo


default_clusterbuster_path = 'cbust'

default_min_crm_score_threshold = 0.00
default_min_delta_crm_score_threshold = 0.00

default_min_motif_score_threshold = 0.00
default_min_delta_motif_score_threshold = 0.00


class ClusterBusterDeltaScore:
    """
    Create ClusterBusterDeltaScore object.
    """

    def __init__(self, wt_crm_score, mut_crm_score, wt_motif_score, mut_motif_score, wt_consensus='', mut_consensus=''):
        """
        Create ClusterBusterDeltaScore object.

        :param wt_crm_score: wildtype Cluster-Buster CRM score.
        :param mut_crm_score: mutant Cluster-Buster CRM score.
        :param wt_motif_score: wildtype Cluster-Buster motif score.
        :param mut_motif_score: mutant Cluster-Buster motif score.
        :param wt_consensus: wildtype consensus sequence.
        :param mut_consensus: mutant consensus sequence.
        :return:
        """

        if wt_crm_score < 0.0:
            raise ValueError(
                'Invalid wildtype Cluster-Buster CRM score {0:f}. The score should be greater or equal to 0.0.'.format(
                    wt_crm_score)
            )
        if mut_crm_score < 0.0:
            raise ValueError(
                'Invalid mutant Cluster-Buster CRM score {0:f}. The score should be greater or equal to 0.0.'.format(
                    mut_crm_score)
            )
        if wt_motif_score < 0.0:
            raise ValueError(
                'Invalid wildtype Cluster-Buster motif score {0:f}. The score should be greater or equal to 0.0.'.format(
                    wt_motif_score)
            )
        if mut_motif_score < 0.0:
            raise ValueError(
                'Invalid mutant Cluster-Buster motif score {0:f}. The score should be greater or equal to 0.0.'.format(
                    mut_motif_score)
            )

        self.wt_crm_score = wt_crm_score
        self.mut_crm_score = mut_crm_score
        self.wt_motif_score = wt_motif_score
        self.mut_motif_score = mut_motif_score
        self.crm_delta_score = mut_crm_score - wt_crm_score
        self.motif_delta_score = mut_motif_score - wt_motif_score
        self.wt_consensus = wt_consensus
        self.mut_consensus = mut_consensus

    def __repr__(self):
        return '<ClusterBusterDeltaScore>\n' \
               '  wt_crm_score: {0:f}\n' \
               '  mut_crm_score: {1:f}\n' \
               '  wt_motif_score: {2:f}\n' \
               '  mut_motif_score: {3:f}\n' \
               '  crm_delta_score: {4:f}\n' \
               '  motif_delta_score: {5:f}\n' \
               '  wt_consensus: {6:s}\n' \
               '  mut_consensus: {7:s}\n'.format(self.wt_crm_score,
                                                 self.mut_crm_score,
                                                 self.wt_motif_score,
                                                 self.mut_motif_score,
                                                 self.crm_delta_score,
                                                 self.motif_delta_score,
                                                 self.wt_consensus,
                                                 self.mut_consensus)

    def is_crm_gain(self,
                    min_mut_crm_score_threshold=default_min_crm_score_threshold,
                    min_delta_crm_score_threshold=default_min_delta_crm_score_threshold):
        """
        Is the Cluster-Buster delta CRM score a CRM gain?

        :param min_mut_crm_score_threshold:
            Minimum mutant Cluster-Buster CRM score.
        :param min_delta_crm_score_threshold:
            Minimum delta CRM score needed before the mutation is considered to be a CRM gain mutation.

        :return: True or False
        """
        return (self.crm_delta_score >= min_delta_crm_score_threshold
                if self.mut_crm_score >= min_mut_crm_score_threshold
                else False)

    def is_crm_loss(self,
                    min_mut_crm_score_threshold=default_min_crm_score_threshold,
                    min_delta_crm_score_threshold=default_min_delta_crm_score_threshold):
        """
        Is the Cluster-Buster delta CRM score a CRM loss?

        :param min_mut_crm_score_threshold:
             Minimum wildtype Cluster-Buster CRM score.
        :param min_delta_crm_score_threshold:
             Minimum delta CRM score needed before the mutation is considered to be a CRM loss mutation.

        :return: True or False
        """
        return ((- self.crm_delta_score) >= min_delta_crm_score_threshold
                if self.wt_crm_score >= min_mut_crm_score_threshold
                else False)

    def is_motif_gain(self,
                      min_mut_motif_score_threshold=default_min_motif_score_threshold,
                      min_delta_motif_score_threshold=default_min_delta_motif_score_threshold):
        """
        Is the Cluster-Buster delta motif score a motif gain?

        :param min_mut_motif_score_threshold:
            Minimum mutant Cluster-Buster motif score.
        :param min_delta_motif_score_threshold:
            Minimum delta motif score needed before the mutation is considered to be a motif gain mutation.

        :return: True or False
        """
        return (self.motif_delta_score >= min_delta_motif_score_threshold
                if self.mut_motif_score >= min_mut_motif_score_threshold
                else False)

    def is_motif_loss(self,
                      min_mut_motif_score_threshold=default_min_motif_score_threshold,
                      min_delta_motif_score_threshold=default_min_delta_motif_score_threshold):
        """
        Is the Cluster-Buster delta motif score a motif loss?

        :param min_mut_motif_score_threshold:
            Minimum wildtype Cluster-Buster motif score.
        :param min_delta_motif_score_threshold:
            Minimum delta motif score needed before the mutation is considered to be a motif loss mutation.

        :return: True or False
        """
        return ((- self.motif_delta_score) >= min_delta_motif_score_threshold
                if self.wt_motif_score >= min_mut_motif_score_threshold
                else False)


def calculate_clusterbuster_delta_scores(vcf_muts,
                                         motif_id,
                                         clusterbuster_path=default_clusterbuster_path,
                                         min_crm_score_threshold=default_min_crm_score_threshold,
                                         min_motif_score_threshold=default_min_motif_score_threshold):
    """
    Calculate Cluster-Buster delta CRM and delta motif scores between wildtype and mutant FASTA sequences, which are
    automatically generated based on the VCFmut objects, for a certain motif.

    :param vcf_muts:
        List of VCFmut objects.
    :param motif_id:
        Motif ID.
    :param clusterbuster_path:
        Path to Cluster-Buster.
    :param min_crm_score_threshold:
        Minimum Cluster-Buster CRM score threshold for wildtype or mutant to keep the result.
    :param min_motif_score_threshold:
        Minimum Cluster-Buster motif score threshold for wildtype or mutant to keep the result.

    :return: Dictionary with VCFmut objects as keys and ClusterBusterDeltaScore objects as values.
    """

    clusterbuster_matrix_filename = motifsinfo.MotifsInfo.get_clusterbuster_motif_filename(motif_id)
    motif_length = motifsinfo.MotifsInfo.get_motif_length(motif_id)

    bp_extension_size = 500 + motif_length

    # Make wildtype and mutant FASTA sequence for each mutation.
    try:
        fasta_string = '\n'.join([vcf_mut.make_fasta_for_wt_and_mut(bp_upstream=bp_extension_size,
                                                                    bp_downstream=bp_extension_size)
                                  for vcf_mut in vcf_muts]).encode('utf-8')
    except ValueError as e:
        raise e

    # Score wildtype and mutant FASTA sequence for each mutation with Cluster-Buster for a specific motif.
    clusterbuster_command = [clusterbuster_path,
                             '-f', '0',
                             '-c', '0.0',
                             '-m', '0.0',
                             '-r', '10000',
                             clusterbuster_matrix_filename,
                             '/dev/stdin']

    (clusterbuster_command_stdout_data, clusterbuster_command_stderr_data) = command.run_cmd(
        cmd=clusterbuster_command,
        stdin=fasta_string
    )

    clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut = dict()

    consider_current_crm = False

    # Parse Cluster-Buster output.
    with BytesIO(clusterbuster_command_stdout_data) as clusterbuster_fh:
        for clusterbuster_line in clusterbuster_fh:
            clusterbuster_line = clusterbuster_line.rstrip()

            if clusterbuster_line.startswith('>'):
                fasta_seq_id = clusterbuster_line.split(' ')[0][1:]
                (vcf_mut,
                 bp_upstream,
                 bp_downstream,
                 is_wt) = mutations.VCFmut.from_fasta_seq_id(fasta_seq_id)

                # The lowest possible motif start position we care about is the same for both wildtype and mutant
                # sequences and starts at that position where the last position of the motif overlaps with the first
                # position of the reference/mutation where the mutation is introduced.
                lowest_motif_start_pos_to_consider = (
                    bp_upstream + 1 -
                    (motif_length - 1)
                )

                if is_wt:
                    # For wildtype sequence.

                    # The highest possible motif end position we care about ends at that position where the first
                    # position of the motif overlaps with the last position of the reference sequence with will be
                    # altered by the mutation.
                    highest_motif_end_pos_to_consider = (
                        (bp_upstream + 1) +
                        (motif_length - 1) +
                        (len(vcf_mut.ref) - 1)
                    )
                else:
                    # For mutated sequence.

                    # The highest possible motif end position we care about ends at that position where the first
                    # position of the motif overlaps with the last position of the mutant sequence with will be altered
                    # by the mutation.
                    highest_motif_end_pos_to_consider = (
                        (bp_upstream + 1) +
                        (motif_length - 1) +
                        (len(vcf_mut.mut) - 1)
                    )

                continue
            elif clusterbuster_line.startswith('Location:'):
                (crm_start_pos_str, crm_end_pos_str) = clusterbuster_line.split(" ")[1:4:2]

                # One-based coordinates.
                crm_start_pos = int(crm_start_pos_str)
                crm_end_pos = int(crm_end_pos_str)

                if (crm_start_pos <= highest_motif_end_pos_to_consider and
                            crm_end_pos >= lowest_motif_start_pos_to_consider):
                    consider_current_crm = True
                else:
                    consider_current_crm = False

                continue
            if consider_current_crm:
                if clusterbuster_line.startswith('Score:'):
                    crm_score = float(clusterbuster_line.split(" ")[1])

                    # Set initial values for wildtype and mutant Cluster-Buster CRM socre, motif score and consensus
                    # sequence for each mutation.
                    clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut.setdefault(
                        vcf_mut,
                        [[0.0, 0.0, ''], [0.0, 0.0, '']]
                    )

                    if is_wt:
                        if clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][0][0] < crm_score:
                            clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][0][0] \
                                = crm_score
                    else:
                        if clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][1][0] < crm_score:
                            clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][1][0] \
                                = crm_score

                    continue
                else:
                    clusterbuster_columns = clusterbuster_line.split('\t')

                    if len(clusterbuster_columns) == 6:
                        # This is a motif line.
                        (motif_id,
                         motif_start_pos_str,
                         motif_end_pos_str,
                         strand,
                         motif_score_str,
                         consensus) = clusterbuster_columns

                        motif_start_pos = int(motif_start_pos_str)
                        motif_end_pos = int(motif_end_pos_str)
                        motif_score = float(motif_score_str)

                        if (motif_start_pos >= lowest_motif_start_pos_to_consider and
                                    motif_end_pos <= highest_motif_end_pos_to_consider):

                            if is_wt:
                                if clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][0][1] < motif_score:
                                    clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][0][1:2] \
                                        = [motif_score, consensus]
                            else:
                                if clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][1][1] < motif_score:
                                    clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[vcf_mut][1][1:2] \
                                        = [motif_score, consensus]

    clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut_delta_above_threshold = dict()

    for vcf_mut, clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut in clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut.iteritems():
        if (clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[0][0] >= min_crm_score_threshold
                or clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[1][0] >= min_crm_score_threshold):
            clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut_delta_above_threshold[vcf_mut] \
                = ClusterBusterDeltaScore(
                    wt_crm_score=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[0][0],
                    mut_crm_score=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[1][0],
                    wt_motif_score=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[0][1],
                    mut_motif_score=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[1][1],
                    wt_consensus=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[0][2],
                    mut_consensus=clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut[1][2],
            )

    return clusterbuster_max_crm_scores_and_max_motif_scores_and_consensus_for_wt_mut_delta_above_threshold
