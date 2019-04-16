"""
Purpose :      Run MotifLocator and get MotifLocator delta scores.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import os.path

from io import BytesIO, TextIOWrapper

import command
import mutations
import motifsinfo


default_motiflocator_path = 'MotifLocator'

default_background_motiflocator_filename = os.path.join(os.path.dirname(__file__),
                                                        'data',
                                                        'motifs',
                                                        'hg19.noncoding_genome.fasta.motifLoc322.bg')

default_min_score_threshold = 0.80
default_min_delta_score_threshold = 0.05


class MotifLocatorDeltaScore:
    """
    Create MotifLocatorDeltaScore object.
    """

    def __init__(self, wt_score, mut_score, wt_consensus='', mut_consensus=''):
        """
        Create MotifLocatorDeltaScore object.

        :param wt_score: wildtype MotifLocator score.
        :param mut_score: mutant MotifLocator score.
        :param wt_consensus: wildtype consensus sequence.
        :param mut_consensus: mutant consensus sequence.
        :return:
        """

        if wt_score < 0.0 or wt_score > 1.0:
            raise ValueError(
                'Invalid wildtype MotifLocator score {0:f}. The score should be between 0.0 and 1.0.'.format(wt_score)
            )
        if mut_score < 0.0 or mut_score > 1.0:
            raise ValueError(
                'Invalid mutant MotifLocator score {0:f}. The score should be between 0.0 and 1.0.'.format(mut_score)
            )

        self.wt_score = wt_score
        self.mut_score = mut_score
        self.delta_score = mut_score - wt_score
        self.wt_consensus = wt_consensus
        self.mut_consensus = mut_consensus

    def __repr__(self):
        return '<MotifLocatorDeltaScore>\n' \
               '  wt_score: {0:f}\n' \
               '  mut_score: {1:f}\n' \
               '  delta_score: {2:f}\n' \
               '  wt_consensus: {3:s}\n' \
               '  mut_consensus: {4:s}\n'.format(self.wt_score,
                                                 self.mut_score,
                                                 self.delta_score,
                                                 self.wt_consensus,
                                                 self.mut_consensus)

    def is_motif_gain(self,
                      min_mut_score_threshold=default_min_score_threshold,
                      min_delta_score_threshold=default_min_delta_score_threshold):
        """
        Is the MotifLocator delta score a motif gain?

        :param min_mut_score_threshold:
            Minimum mutant MotifLocator score.
        :param min_delta_score_threshold:
            Minimum delta score needed before the mutation is considered to be a motif gain mutation.

        :return: True or False
        """
        return (self.delta_score >= min_delta_score_threshold
                if self.mut_score >= min_mut_score_threshold
                else False)

    def is_motif_loss(self,
                      min_wt_score_threshold=default_min_score_threshold,
                      min_delta_score_threshold=default_min_delta_score_threshold):
        """
        Is the MotifLocator delta score a motif loss?

        :param min_wt_score_threshold:
            Minimum wildtype MotifLocator score.
        :param min_delta_score_threshold:
            Minimum delta score needed before the mutation is considered to be a motif loss mutation.

        :return: True or False
        """
        return ((- self.delta_score) >= min_delta_score_threshold
                if self.wt_score >= min_wt_score_threshold
                else False)


def calculate_motiflocator_delta_scores(fasta_string,
                                        inclusive_matrix_filename,
                                        motiflocator_path=default_motiflocator_path,
                                        background_motiflocator_filename=default_background_motiflocator_filename,
                                        min_score_threshold=default_min_score_threshold):
    """
    Calculate MotifLocator delta scores between wildtype and mutant FASTA sequences.

    :param fasta_string: String with wildtype and mutant FASTA sequence constructed by VCFmut.make_fasta_for_wt_and_mut().
    :param inclusive_matrix_filename: Filename with motifs in INCLUsive format to score with MotifLocator.
    :param motiflocator_path: Path to MotifLocator.
    :param background_motiflocator_filename: Background file used by MotifLocator.
    :param min_score_threshold: Minimum MotifLocator score threshold for wildtype or mutant to keep the result.

    :return: Dictionary with motifs as keys and  MotifLocatorDeltaScore objects as values.
    """

    motiflocator_command = [motiflocator_path,
                            '-f', '/dev/stdin',
                            '-b', background_motiflocator_filename,
                            '-m', inclusive_matrix_filename,
                            '-t', '0.0',
                            '-s', '1',
                            '-a', '0']

    (motiflocator_command_stdout_data, motiflocator_command_stderr_data) = command.run_cmd(
        cmd=motiflocator_command,
        stdin=fasta_string.encode('utf-8')
    )

    motiflocator_max_scores_and_consensus_for_wt_mut = dict()
    fasta_seq_id_to_from_fasta_seq_id_output_dict = dict()

    with TextIOWrapper(buffer=BytesIO(motiflocator_command_stdout_data), encoding='utf-8') as gff_fh:
        for gff_line in gff_fh:
            columns = gff_line.rstrip().split('\t')

            if len(columns) == 9:
                # FASTA sequence ID constructed by VCFmut.make_fasta_for_wt_and_mut().
                fasta_seq_id = columns[0]

                # Motif start and end positions (one-based).
                motif_start_pos = int(columns[3])
                motif_end_pos = int(columns[4])

                # MotifLocator score.
                score = float(columns[5])

                # Get motif name and consensus sequence.
                motif_name, consensus = columns[8].split('"')[1:4:2]

                # Motif ID.
                motif_id = motifsinfo.MotifsInfo.get_motif_id(motif_name)

                if fasta_seq_id not in fasta_seq_id_to_from_fasta_seq_id_output_dict:
                    # Extract info from the FASTA sequence ID constructed by VCFmut.make_fasta_for_wt_and_mut() only
                    # when we see an new FASTA sequence ID.
                    fasta_seq_id_to_from_fasta_seq_id_output_dict[fasta_seq_id] = \
                        mutations.VCFmut.from_fasta_seq_id(fasta_seq_id)

                    (vcf_mut,
                     bp_upstream,
                     bp_downstream,
                     is_wt) = fasta_seq_id_to_from_fasta_seq_id_output_dict[fasta_seq_id]

                # Set initial values for wildtype and mutant MotifLocator score and consensus sequence for each motif.
                motiflocator_max_scores_and_consensus_for_wt_mut.setdefault(motif_id, [[0.0, ''], [0.0, '']])

                # The lowest possible motif start position we care about is the same for both wildtype and mutant
                # sequences and starts at that position where the last position of the motif overlaps with the first
                # position of the reference/mutation where the mutation is introduced.
                lowest_motif_start_pos_to_consider = (
                    bp_upstream + 1 -
                    (motifsinfo.MotifsInfo.get_motif_length(motif_id) - 1)
                )

                if is_wt:
                    # For wildtype sequence.

                    # The highest possible motif end position we care about, ends at that position where the first
                    # position of the motif overlaps with the last position of the reference sequence which will be
                    # altered by the mutation.
                    highest_motif_end_pos_to_consider = (
                        (bp_upstream + 1) +
                        (motifsinfo.MotifsInfo.get_motif_length(motif_id) - 1) +
                        (len(vcf_mut.ref) - 1)
                    )

                    if (motif_start_pos >= lowest_motif_start_pos_to_consider and
                                motif_end_pos <= highest_motif_end_pos_to_consider):
                        if motiflocator_max_scores_and_consensus_for_wt_mut[motif_id][0][0] < score:
                            # Update maximum MotifLocator score and consensus sequence for wildtype sequence for the
                            # current motif.
                            motiflocator_max_scores_and_consensus_for_wt_mut[motif_id][0] = [score, consensus]
                else:
                    # For mutated sequence.

                    # The highest possible motif end position we care about, ends at that position where the first
                    # position of the motif overlaps with the last position of the mutant sequence which will be altered
                    # by the mutation.
                    highest_motif_end_pos_to_consider = (
                        (bp_upstream + 1) +
                        (motifsinfo.MotifsInfo.get_motif_length(motif_id) - 1) +
                        (len(vcf_mut.mut) - 1)
                    )

                    if (motif_start_pos >= lowest_motif_start_pos_to_consider and
                                motif_end_pos <= highest_motif_end_pos_to_consider):
                        if motiflocator_max_scores_and_consensus_for_wt_mut[motif_id][1][0] < score:
                            # Update maximum MotifLocator score and consensus sequence for mutant sequence for the
                            # current motif.
                            motiflocator_max_scores_and_consensus_for_wt_mut[motif_id][1] = [score, consensus]

    motiflocator_max_scores_wt_mut_delta_above_threshold = dict()

    for motif_id, max_scores_and_consencus_for_wt_mut in motiflocator_max_scores_and_consensus_for_wt_mut.items():
        if (max_scores_and_consencus_for_wt_mut[0][0] >= min_score_threshold
                or max_scores_and_consencus_for_wt_mut[1][0] >= min_score_threshold):
            motiflocator_max_scores_wt_mut_delta_above_threshold[motif_id] = MotifLocatorDeltaScore(
                wt_score=max_scores_and_consencus_for_wt_mut[0][0],
                mut_score=max_scores_and_consencus_for_wt_mut[1][0],
                wt_consensus=max_scores_and_consencus_for_wt_mut[0][1],
                mut_consensus=max_scores_and_consencus_for_wt_mut[1][1],
            )

    return motiflocator_max_scores_wt_mut_delta_above_threshold
