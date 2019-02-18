#!/usr/bin/env python3

"""
Purpose :      When importing this file:
                 MotifsInfo class for retrieving info related to motifs:
                   - Get length of a motif.
                   - Get motif filename for a motif.
                   - Get directly annotated TFs for a motif.
                   - Get motifs directly annotated for a certain TF.

               When running as standalone script write 4 TSV files:
                 Create motif ID to TFs association TSV files:
                   - with each motif ID and all associated TFs on one line.
                   - with each motif ID and TF pair on a separate line.
                 Create TF to motif IDs association TSV files:
                   - with each TF and all associated motif IDs on one line .
                   - with each TF and TF pair on a separate line.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import argparse
import os


class MotifPaths:
    # Class variables need to be initialised by calling:
    #     MotifPaths.set_motif_collection_version(motif_collection_version)
    motif_collection_version = None
    default_motif_to_tf_filename = None
    clusterbuster_directly_annotated_motifs_dir = None
    clusterbuster_all_motifs_dir = None
    inclusive_directly_annotated_motifs_dir = None
    inclusive_all_motifs_dir = None

    # Extension used for a motif filename in Cluster-Buster format.
    default_clusterbuster_motifs_extension = '.cb'

    # Extension used for a motif filename in INCLUSive format.
    default_inclusive_motifs_extension = '.INCLUsive.txt'

    @staticmethod
    def set_motif_collection_version(motif_collection_version):
        MotifPaths.motif_collection_version = motif_collection_version

        # Filename with motif to TF annotation (snapshot from motif to TF database).
        MotifPaths.default_motif_to_tf_filename = os.path.join(
            os.path.dirname(__file__),
            'data',
            'motifs',
            motif_collection_version,
            'directly_annotated_motifs',
            'motifs-{0:s}-nr.hgnc-m0.001-o0.0.tbl'.format(motif_collection_version))

        # Directory with directly annotated motifs in Cluster-Buster format.
        MotifPaths.clusterbuster_directly_annotated_motifs_dir = os.path.join(
            os.path.dirname(__file__),
            'data',
            'motifs',
            motif_collection_version,
            'directly_annotated_motifs',
            'clusterbuster')

        # Directory with all motifs in Cluster-Buster format.
        MotifPaths.clusterbuster_all_motifs_dir = os.path.join(
            os.path.dirname(__file__),
            'data',
            'motifs',
            motif_collection_version,
            'all_motifs',
            'clusterbuster')

        # Directory with directly annotated motifs in INCLUSive format.
        MotifPaths.inclusive_directly_annotated_motifs_dir = os.path.join(
            os.path.dirname(__file__),
            'data',
            'motifs',
            motif_collection_version,
            'directly_annotated_motifs',
            'inclusive')

        # Directory with all motifs in INCLUSive format.
        MotifPaths.inclusive_all_motifs_dir = os.path.join(
            os.path.dirname(__file__),
            'data',
            'motifs',
            motif_collection_version,
            'all_motifs',
            'inclusive')


def get_motif_name_and_motif_filenames_and_motif_lengths_and_pwms(
        clusterbuster_motifs_dir=MotifPaths.clusterbuster_directly_annotated_motifs_dir,
        clusterbuster_motifs_extension=MotifPaths.default_clusterbuster_motifs_extension,
        inclusive_motifs_dir=MotifPaths.inclusive_directly_annotated_motifs_dir,
        inclusive_motifs_extension=MotifPaths.default_inclusive_motifs_extension):
    """
    Get motif names and motif filenames and motif lengths and PWMs.

    :param inclusive_motifs_dir: Directory with motifs in INCLUsive format.
    :param inclusive_motifs_extension: Extension used for a motif filename in INCLUSive format.

    :return: (motif_id_to_motif_name_dict,
              motif_name_to_motif_id_dict,
              motif_id_to_motif_length_dict,
              motif_id_to_clusterbuster_filename_dict,
              motif_id_to_inclusive_filename_dict,
              motif_id_to_inclusive_pwm_dict)
    """

    motif_id_to_motif_name_dict = dict()
    motif_name_to_motif_id_dict = dict()
    motif_id_to_motif_length_dict = dict()
    motif_id_to_clusterbuster_filename_dict = dict()
    motif_id_to_inclusive_filename_dict = dict()
    motif_id_to_inclusive_pwm_dict = dict()

    # Get all motif files in Cluster-Buster format from the Cluster-Buster directory.
    for folder, subdirs, filenames in os.walk(clusterbuster_motifs_dir):
        for filename in filenames:
            if filename.endswith(clusterbuster_motifs_extension):
                clusterbuster_motif_filename = os.path.join(clusterbuster_motifs_dir, filename)

                motif_id = filename[0:- len(clusterbuster_motifs_extension)]

                motif_id_to_clusterbuster_filename_dict[motif_id] = clusterbuster_motif_filename

    # Get all motif files in INCLUSive format from the INCLUSive directory.
    for folder, subdirs, filenames in os.walk(inclusive_motifs_dir):
        for filename in filenames:
            if filename.endswith(inclusive_motifs_extension):
                inclusive_motif_filename = os.path.join(inclusive_motifs_dir, filename)

                motif_id = filename[0:- len(inclusive_motifs_extension)]

                inclusive_pwm = ''

                with open(inclusive_motif_filename, 'r') as fh:
                    for full_line in fh:
                        line = full_line.rstrip('\n')

                        if line != '#INCLUSive Motif Model' and line != '':
                            inclusive_pwm += full_line

                            columns = line.split(' ')

                            if len(columns) == 3 and columns[1] == '=':
                                if columns[0] == '#ID':
                                    motif_name = columns[2]

                                    motif_id_to_motif_name_dict[motif_id] = motif_name
                                    motif_name_to_motif_id_dict[motif_name] = motif_id
                                elif columns[0] == '#W':
                                    motif_length = columns[2]

                                    motif_id_to_inclusive_filename_dict[motif_id] = inclusive_motif_filename
                                    motif_id_to_motif_length_dict[motif_id] = int(motif_length)

                    motif_id_to_inclusive_pwm_dict[motif_id] = inclusive_pwm

    if set(motif_id_to_clusterbuster_filename_dict) != set(motif_id_to_inclusive_filename_dict):
        raise ValueError(
            'The Cluster-Buster and INCLUSive motif directory does not contain the same motifs. '
            'Motif IDs which are unique for Cluster-Buster: {0:s}. '
            'Motif IDs which are unique for INCLUSive: {1:s}'.format(
                ' '.join(
                    set(motif_id_to_clusterbuster_filename_dict).difference(
                        set(motif_id_to_inclusive_filename_dict))
                ),
                ' '.join(
                    set(motif_id_to_inclusive_filename_dict).difference(
                        set(motif_id_to_clusterbuster_filename_dict)
                    )
                )
            )
        )

    return (motif_id_to_motif_name_dict,
            motif_name_to_motif_id_dict,
            motif_id_to_motif_length_dict,
            motif_id_to_clusterbuster_filename_dict,
            motif_id_to_inclusive_filename_dict,
            motif_id_to_inclusive_pwm_dict)


def get_direct_motif_to_tf_annotation(motif_to_tf_filename=MotifPaths.default_motif_to_tf_filename,
                                      motif_ids_to_consider=None):
    """
    Get motif to TF annotation for motifs which are directly annotated.

    :param motif_to_tf_filename:
                Motif to TF filename.
    :param motif_ids_to_consider:
                if not None, only add direct motif to TF annotation
                if motif ID occurs is in this list.

    :return: motif_to_tfs_dict, tf_to_motifs_dict
    """

    motif_to_tfs_dict = dict()
    tf_to_motifs_dict = dict()

    with open(motif_to_tf_filename, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            columns = line.rstrip('\n').split('\t')

            if len(columns) == 4:
                if columns[3] == 'gene is directly annotated':
                    motif_id = columns[0]
                    tf = columns[2]

                    if motif_ids_to_consider is None:
                        # Add all motifs which have direct motif to TF annotation.
                        motif_to_tfs_dict.setdefault(motif_id, set()).add(tf)
                        tf_to_motifs_dict.setdefault(tf, set()).add(motif_id)
                    elif motif_id in motif_ids_to_consider:
                        # Only add motifs which have direct motif to TF annotation if it appears in the list of motifs.
                        motif_to_tfs_dict.setdefault(motif_id, set()).add(tf)
                        tf_to_motifs_dict.setdefault(tf, set()).add(motif_id)
            elif len(columns) == 13:
                if columns[12] == 'gene is directly annotated':
                    motif_id = columns[0]
                    tf = columns[5]

                    if motif_ids_to_consider is None:
                        # Add all motifs which have direct motif to TF annotation.
                        motif_to_tfs_dict.setdefault(motif_id, set()).add(tf)
                        tf_to_motifs_dict.setdefault(tf, set()).add(motif_id)
                    elif motif_id in motif_ids_to_consider:
                        # Only add motifs which have direct motif to TF annotation if it appears in the list of motifs.
                        motif_to_tfs_dict.setdefault(motif_id, set()).add(tf)
                        tf_to_motifs_dict.setdefault(tf, set()).add(motif_id)

    return motif_to_tfs_dict, tf_to_motifs_dict


class MotifsInfo:
    """
    Class for retrieving info related to motifs:
      - Get motif name.
      - Get length of a motif.
      - Get motif filename for a motif in Cluster-Buster and INCLUSive format.
      - Get directly annotated TFs for a motif.
      - Get motifs directly annotated for a certain TF.
    """

    # Class variables need to be initialised by calling:
    #     MotifsInfo.set_motif_collection_version(motif_collection_version)
    motif_collection_version = None
    motif_id_to_motif_name_dict = dict()
    motif_name_to_motif_id_dict = dict()
    motif_id_to_motif_length_dict = dict()
    motif_id_to_clusterbuster_filename_dict = dict()
    motif_id_to_inclusive_filename_dict = dict()
    motif_id_to_inclusive_pwm_dict = dict()

    directly_annotated_motifs_set = set()

    motif_to_tfs_dict = dict()
    tf_to_motifs_dict = dict()

    @staticmethod
    def set_motif_collection_version(motif_collection_version):
        if MotifsInfo.motif_collection_version != motif_collection_version:
            MotifPaths.set_motif_collection_version(motif_collection_version=motif_collection_version)

            # Create dictionaries with information about all motifs.
            (MotifsInfo.motif_id_to_motif_name_dict,
             MotifsInfo.motif_name_to_motif_id_dict,
             MotifsInfo.motif_id_to_motif_length_dict,
             MotifsInfo.motif_id_to_clusterbuster_filename_dict,
             MotifsInfo.motif_id_to_inclusive_filename_dict,
             MotifsInfo.motif_id_to_inclusive_pwm_dict) = get_motif_name_and_motif_filenames_and_motif_lengths_and_pwms(
                clusterbuster_motifs_dir=MotifPaths.clusterbuster_all_motifs_dir,
                clusterbuster_motifs_extension=MotifPaths.default_clusterbuster_motifs_extension,
                inclusive_motifs_dir=MotifPaths.inclusive_all_motifs_dir,
                inclusive_motifs_extension=MotifPaths.default_inclusive_motifs_extension
            )

            # Get directly annotated motif names, from motif_id_to_inclusive_filename_dict by getting all INCLUSive
            # motif files in inclusive_directly_annotated_motifs_dir.
            MotifsInfo.directly_annotated_motifs_set = set(
                get_motif_name_and_motif_filenames_and_motif_lengths_and_pwms(
                    clusterbuster_motifs_dir=MotifPaths.clusterbuster_directly_annotated_motifs_dir,
                    clusterbuster_motifs_extension=MotifPaths.default_clusterbuster_motifs_extension,
                    inclusive_motifs_dir=MotifPaths.inclusive_directly_annotated_motifs_dir,
                    inclusive_motifs_extension=MotifPaths.default_inclusive_motifs_extension
                )[4]
            )

            MotifsInfo.motif_to_tfs_dict, MotifsInfo.tf_to_motifs_dict = get_direct_motif_to_tf_annotation(
                MotifPaths.default_motif_to_tf_filename,
                motif_ids_to_consider=MotifsInfo.motif_id_to_inclusive_filename_dict
            )

    @staticmethod
    def get_all_motif_ids(directly_annotated_motifs_only):
        """
        Get all motifs IDs (directly_annotated_motifs_only=False) or
        get all directly annotated motifs IDS (directly_annotated_motifs_only=True).

        :param directly_annotated_motifs_only: True or False
        :return: set of motif IDs
        """
        if directly_annotated_motifs_only:
            return MotifsInfo.directly_annotated_motifs_set
        else:
            return set(MotifsInfo.motif_id_to_inclusive_filename_dict)

    @staticmethod
    def get_motif_name(motif_id):
        """
        Get motif name from a motif ID.

        :param motif_id: motif ID.
        :return: motif name.
        """
        return MotifsInfo.motif_id_to_motif_name_dict.get(motif_id, None)

    @staticmethod
    def get_motif_id(motif_name):
        """
        Get motif ID from motif name.

        :param motif_name: motif name.
        :return: motif ID.
        """
        return MotifsInfo.motif_name_to_motif_id_dict.get(motif_name, None)

    @staticmethod
    def get_motif_length(motif_id):
        """
        Get length of a motif.

        :param motif_id: motif ID.
        :return: length of the motif.
        """
        return MotifsInfo.motif_id_to_motif_length_dict.get(motif_id, None)

    @staticmethod
    def get_clusterbuster_motif_filename(motif_id):
        """
        Get motif filename for a motif in Cluster-Buster format.

        :param motif_id: motif ID.
        :return: Cluster-Buster motif filename.
        """
        return MotifsInfo.motif_id_to_clusterbuster_filename_dict.get(motif_id, None)

    @staticmethod
    def get_inclusive_motif_filename(motif_id):
        """
        Get motif filename for a motif in INCLUSive format.

        :param motif_id: motif ID.
        :return: INCLUSive motif filename.
        """
        return MotifsInfo.motif_id_to_inclusive_filename_dict.get(motif_id, None)

    @staticmethod
    def get_tfs_for_motif(motif_id):
        """
        Get directly annotated TFs for a motif.

        :param motif_id: motif ID.
        :return: set of transcription factor names.
        """
        return MotifsInfo.motif_to_tfs_dict.get(motif_id, None)

    @staticmethod
    def get_motifs_for_tf(tf):
        """
        Get motifs directly annotated for a certain TF.

        :param tf: transcription factor name.
        :return: set of motif IDs.
        """
        return MotifsInfo.tf_to_motifs_dict.get(tf, None)

    @staticmethod
    def get_inclusive_pwm(motif_id, header=False):
        """
        Get PWM in INCLUSive format for a motif.

        :param motif_id: motif ID.
        :param header: Add "#INCLUSive Motif Model" header line (True) or not (False).
        :return: string with PWM in INCLUSive format for motif ID.
        """
        if header:
            return '#INCLUSive Motif Model\n\n' + MotifsInfo.motif_id_to_inclusive_pwm_dict[motif_id] + '\n'
        else:
            return MotifsInfo.motif_id_to_inclusive_pwm_dict[motif_id] + '\n'

    @staticmethod
    def get_inclusive_pwms(motif_ids, min_motif_length=None, max_motif_length=None, header=True):
        """
        Get PWMs in INCLUSive format for a list of motifs.

        :param motif_ids: list of motif IDs.
                          If set to None, all motif IDs in MotifsInfo.motif_id_to_filename_dict are used.
        :param min_motif_length: Only include PWMs which have a minimum motif length of min_motif_length
                                 or do not use this restriction if this parameter is set to None.
        :param max_motif_length: Only include PWMs which have a maximum motif length of max_motif_length
                                 or do not use this restriction if this parameter is set to None.
        :param header: Add "#INCLUSive Motif Model" header line (True) or not (False).
        :return: string with PWMs in INCLUSive format for motif IDs,
                 list of motif IDs that passed the filtering,
                 are there motif IDs that passed the filtering (True or False)
        """

        if not motif_ids:
            # Set motif_ids to sorted list of all motifs if motifs_ids is None.
            motif_ids = sorted(MotifsInfo.motif_id_to_inclusive_filename_dict)

        selected_motif_ids = [
            motif_id
            for motif_id in motif_ids
            if (min_motif_length is None or MotifsInfo.get_motif_length(motif_id) >= min_motif_length) and
               (max_motif_length is None or MotifsInfo.get_motif_length(motif_id) <= max_motif_length)
        ]

        if header:
            inclusive_pwms = '#INCLUSive Motif Model\n\n'
        else:
            inclusive_pwms = ''

        inclusive_pwms += ''.join([MotifsInfo.get_inclusive_pwm(motif_id)
                                   for motif_id in selected_motif_ids])

        # Return:
        #   - string with PWMs in INCLUSive format.
        #   - list of motif IDs that passed the filtering step.
        #   - are there motif IDs that passed the filtering (True or False).
        return inclusive_pwms, selected_motif_ids, True if len(selected_motif_ids) else False

    @staticmethod
    def write_motif_id_and_name_to_tfs_annotation_filename(motif_id_and_name_to_tfs_annotation_filename,
                                                           separate_motifs_by='tab'):
        with open(motif_id_and_name_to_tfs_annotation_filename, 'w') as fh:
            if separate_motifs_by == 'tab':
                print('# motif ID\tmotif name\tTFs', file=fh)

                for motif_id in sorted(MotifsInfo.motif_to_tfs_dict):
                    print(motif_id,
                          MotifsInfo.get_motif_name(motif_id),
                          '\t'.join(sorted(MotifsInfo.motif_to_tfs_dict[motif_id])),
                          sep='\t',
                          file=fh)
            elif separate_motifs_by == 'line':
                print('# motif ID\tmotif name\tTF', file=fh)

                for motif_id in sorted(MotifsInfo.motif_to_tfs_dict):
                    for tf in sorted(MotifsInfo.motif_to_tfs_dict[motif_id]):
                        print(motif_id,
                              MotifsInfo.get_motif_name(motif_id),
                              tf,
                              sep='\t',
                              file=fh)

    @staticmethod
    def write_tf_to_motif_ids_annotation_filename(tf_to_motif_ids_annotation_filename, separate_tfs_by='tab'):
        with open(tf_to_motif_ids_annotation_filename, 'w') as fh:
            if separate_tfs_by == 'tab':
                print('# TF\tmotif IDs', file=fh)

                for tf in sorted(MotifsInfo.tf_to_motifs_dict):
                    print(tf,
                          '\t'.join(sorted(MotifsInfo.tf_to_motifs_dict[tf])),
                          sep='\t',
                          file=fh)
            elif separate_tfs_by == 'line':
                print('# TF\tmotif ID', file=fh)

                for tf in sorted(MotifsInfo.tf_to_motifs_dict):
                    for motif_id in sorted(MotifsInfo.tf_to_motifs_dict[tf]):
                        print(tf,
                              motif_id,
                              sep='\t',
                              file=fh)

    @staticmethod
    def write_tf_to_motif_names_annotation_filename(tf_to_motif_names_annotation_filename, separate_tfs_by='tab'):
        with open(tf_to_motif_names_annotation_filename, 'w') as fh:
            if separate_tfs_by == 'tab':
                print('# TF\tmotif names', file=fh)

                for tf in sorted(MotifsInfo.tf_to_motifs_dict):
                    print(tf,
                          '\t'.join(
                              sorted([MotifsInfo.motif_id_to_motif_name_dict[motif_id]
                                      for motif_id in MotifsInfo.tf_to_motifs_dict[tf]])
                          ),
                          sep='\t',
                          file=fh)
            elif separate_tfs_by == 'line':
                print('# TF\tmotif name', file=fh)

                for tf in sorted(MotifsInfo.tf_to_motifs_dict):
                    for motif_id in sorted(MotifsInfo.tf_to_motifs_dict[tf]):
                        print(tf,
                              MotifsInfo.motif_id_to_motif_name_dict[motif_id],
                              sep='\t',
                              file=fh)


class FilterINCLUSiveMotifsOnLength:
    def __init__(self, motif_ids, bp_upstream=0, bp_downstream=0, min_motif_length=None, max_motif_length=None,
                 header=True, inclusive_matrix_fh=None):
        """
        Filter list of motif IDs based on their motif length.
        Get PWMs in INCLUSive format for a list of motifs.

        :param motif_ids: list of motif IDs.
                          If set to None, all motif IDs in MotifsInfo.motif_id_to_filename_dict are used.
        :param bp_upstream: Set number of base pairs upstream of a mutation start position.
        :param bp_upstream: Set number of base pairs downstream of a mutation end position.
        :param min_motif_length: Only include PWMs which have a minimum motif length of min_motif_length
                                 or do not use this restriction if this parameter is set to None.
        :param max_motif_length: Only include PWMs which have a maximum motif length of max_motif_length
                                 or do not use this restriction if this parameter is set to None.
        :param header: Add "#INCLUSive Motif Model" header line (True) or not (False).
        :param inclusive_matrix_fh: file handle to which the string with PMWs in INCLUSive format will be written.
        :return:
        """

        self.inclusive_pwms, self.motif_ids, self.has_motif_ids = MotifsInfo.get_inclusive_pwms(
            motif_ids=motif_ids,
            min_motif_length=min_motif_length,
            max_motif_length=max_motif_length,
            header=header
        )

        self.bp_upstream = bp_upstream
        self.bp_downstream = bp_downstream
        self.min_motif_length = min_motif_length
        self.max_motif_length = max_motif_length
        self.header = header
        self.inclusive_matrix_fh = inclusive_matrix_fh

    def write_matrix_file(self):
        """
        Write PWMs in INCLUSive format to a file handle.

        :return:
        """
        self.inclusive_matrix_fh.write(self.inclusive_pwms)
        self.inclusive_matrix_fh.flush()


def main():
    parser = argparse.ArgumentParser(
        description='Generate motif (ID and name) to TF and TF to motif (ID and name) TSV files.'
    )

    parser.add_argument(
        '--motifcollection',
        dest='motif_collection_version',
        action='store',
        type=str,
        required=True,
        choices=['v7', 'v8'],
        help='Motif collection to use: "v7" or "v8".'
    )

    args = parser.parse_args()

    # Fill MotifsInfo class with content for the chosen motif collection version.
    MotifsInfo.set_motif_collection_version(motif_collection_version=args.motif_collection_version)

    motif_id_and_name_to_tfs_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'motif_id_and_name_to_tfs_annotation.tab.tsv'
    )
    motif_id_and_name_to_tfs_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'motif_id_and_name_to_tfs_annotation.line.tsv'
    )
    tf_to_motif_ids_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'tf_to_motif_ids_annotation.tab.tsv'
    )
    tf_to_motif_ids_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'tf_to_motif_ids_annotation.line.tsv'
    )
    tf_to_motif_names_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'tf_to_motif_names_annotation.tab.tsv'
    )
    tf_to_motif_names_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'motifs',
        args.motif_collection_version,
        'directly_annotated_motifs',
        'tf_to_motif_names_annotation.line.tsv'
    )

    print('Write motif ID and name to TFs association file "{0:s}" with each motif ID and all associated TFs on one line.'.format(
        motif_id_and_name_to_tfs_annotation_tab_filename)
    )
    MotifsInfo.write_motif_id_and_name_to_tfs_annotation_filename(
        motif_id_and_name_to_tfs_annotation_filename=motif_id_and_name_to_tfs_annotation_tab_filename,
        separate_motifs_by='tab'
    )

    print('Write motif ID and name to TFs association file "{0:s}" with each motif ID and TF pair on a separate line.'.format(
        motif_id_and_name_to_tfs_annotation_line_filename)
    )
    MotifsInfo.write_motif_id_and_name_to_tfs_annotation_filename(
        motif_id_and_name_to_tfs_annotation_filename=motif_id_and_name_to_tfs_annotation_line_filename,
        separate_motifs_by='line'
    )

    print('Write TF to motif IDs association file "{0:s}" with each TF and all associated motif IDs on one line.'.format(
        tf_to_motif_ids_annotation_tab_filename)
    )
    MotifsInfo.write_tf_to_motif_ids_annotation_filename(
        tf_to_motif_ids_annotation_filename=tf_to_motif_ids_annotation_tab_filename,
        separate_tfs_by='tab'
    )

    print('Write TF to motif IDs association file "{0:s}" with each TF and motif ID pair on a separate line.'.format(
        tf_to_motif_ids_annotation_line_filename)
    )
    MotifsInfo.write_tf_to_motif_ids_annotation_filename(
        tf_to_motif_ids_annotation_filename=tf_to_motif_ids_annotation_line_filename,
        separate_tfs_by='line'
    )

    print('Write TF to motif names association file "{0:s}" with each TF and all associated motif names on one line.'.format(
        tf_to_motif_names_annotation_tab_filename)
    )
    MotifsInfo.write_tf_to_motif_names_annotation_filename(
        tf_to_motif_names_annotation_filename=tf_to_motif_names_annotation_tab_filename,
        separate_tfs_by='tab'
    )

    print('Write TF to motif names association file "{0:s}" with each TF and motif name pair on a separate line.'.format(
        tf_to_motif_names_annotation_line_filename)
    )
    MotifsInfo.write_tf_to_motif_names_annotation_filename(
        tf_to_motif_names_annotation_filename=tf_to_motif_names_annotation_line_filename,
        separate_tfs_by='line'
    )


if __name__ == "__main__":
    main()
