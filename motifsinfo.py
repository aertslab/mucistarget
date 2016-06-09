#!/usr/bin/env python

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

Copyright (C): 2016 - Gert Hulselmans
"""

from __future__ import print_function

import os


# Filename with motif to TF annotation (snapshot from motif to TF database).
default_motif_to_tf_filename = os.path.join(os.path.dirname(__file__),
                                            'data',
                                            'directly_annotated_motifs',
                                            'motifs-v7-nr.hgnc-m0.001-o0.0.tbl')

# Directory with all motifs in MotifLocator (INCLUSive) format.
default_motif_locator_motifs_dir = os.path.join(os.path.dirname(__file__),
                                                'data',
                                                'directly_annotated_motifs',
                                                'motif_locator')

# Extension used for a motif filename in MotifLocator (INCLUSive) format.
default_motif_locator_motifs_extension = '.INCLUsive.txt'


def get_motif_name_and_motif_filenames_and_motif_lengths(
        motif_locator_motifs_dir=default_motif_locator_motifs_dir,
        motif_locator_motifs_extension=default_motif_locator_motifs_extension):
    """
    Get motif names and motif filenames and motif lengths.

    :param motif_locator_motifs_dir: Directory with motifs in MotifLocator (INCLUsive) format.
    :param motif_locator_motifs_extension: Extension used for a motif filename in MotifLocator (INCLUSive) format.

    :return: motif_id_to_motif_name_dict, motif_name_to_motif_id_dict, motif_id_to_filename_dict, motif_id_to_motif_length_dict
    """

    motif_id_to_motif_name_dict = dict()
    motif_name_to_motif_id_dict = dict()
    motif_id_to_filename_dict = dict()
    motif_id_to_motif_length_dict = dict()

    # Get all motif files in MotifLocator (INCLUSive) format from the MotifLocator directory.
    for folder, subdirs, filenames in os.walk(motif_locator_motifs_dir):
        for filename in filenames:
            if filename.endswith(motif_locator_motifs_extension):
                motif_locator_motif_filename = os.path.join(motif_locator_motifs_dir, filename)

                motif_id = filename[0:- len(motif_locator_motifs_extension)]

                with open(motif_locator_motif_filename, 'r') as fh:
                    for line in fh:
                        columns = line.rstrip('\n').split(' ')

                        if len(columns) == 3 and columns[1] == '=':
                            if columns[0] == '#ID':
                                motif_name = columns[2]

                                motif_id_to_motif_name_dict[motif_id] = motif_name
                                motif_name_to_motif_id_dict[motif_name] = motif_id
                            elif columns[0] == '#W':
                                motif_length = columns[2]

                                motif_id_to_filename_dict[motif_id] = motif_locator_motif_filename
                                motif_id_to_motif_length_dict[motif_id] = int(motif_length)

    return (motif_id_to_motif_name_dict,
            motif_name_to_motif_id_dict,
            motif_id_to_filename_dict,
            motif_id_to_motif_length_dict)


def get_direct_motif_to_tf_annotation(motif_to_tf_filename=default_motif_to_tf_filename,
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

    return motif_to_tfs_dict, tf_to_motifs_dict


class MotifsInfo:
    """
    Class for retrieving info related to motifs:
      - Get motif name.
      - Get length of a motif.
      - Get motif filename for a motif.
      - Get directly annotated TFs for a motif.
      - Get motifs directly annotated for a certain TF.
    """

    (motif_id_to_motif_name_dict,
     motif_name_to_motif_id_dict,
     motif_id_to_filename_dict,
     motif_id_to_motif_length_dict) = get_motif_name_and_motif_filenames_and_motif_lengths(
        default_motif_locator_motifs_dir,
        default_motif_locator_motifs_extension
    )

    motif_to_tfs_dict, tf_to_motifs_dict = get_direct_motif_to_tf_annotation(
        default_motif_to_tf_filename,
        motif_ids_to_consider=motif_id_to_filename_dict
    )

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
    def get_motif_filename(motif_id):
        """
        Get motif filename for a motif.

        :param motif_id: motif ID.
        :return: motif filename.
        """
        return MotifsInfo.motif_id_to_filename_dict.get(motif_id, None)

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
                    for tf in MotifsInfo.motif_to_tfs_dict[motif_id]:
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
                    for motif_id in MotifsInfo.tf_to_motifs_dict[tf]:
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
                    for motif_id in MotifsInfo.tf_to_motifs_dict[tf]:
                        print(tf,
                              MotifsInfo.motif_id_to_motif_name_dict[motif_id],
                              sep='\t',
                              file=fh)


def main():
    motif_id_and_name_to_tfs_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'directly_annotated_motifs',
        'motif_id_and_name_to_tfs_annotation.tab.tsv'
    )
    motif_id_and_name_to_tfs_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'directly_annotated_motifs',
        'motif_id_and_name_to_tfs_annotation.line.tsv'
    )
    tf_to_motif_ids_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'directly_annotated_motifs',
        'tf_to_motif_ids_annotation.tab.tsv'
    )
    tf_to_motif_ids_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'directly_annotated_motifs',
        'tf_to_motif_ids_annotation.line.tsv'
    )
    tf_to_motif_names_annotation_tab_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
        'directly_annotated_motifs',
        'tf_to_motif_names_annotation.tab.tsv'
    )
    tf_to_motif_names_annotation_line_filename = os.path.join(
        os.path.dirname(__file__),
        'data',
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
