"""
Purpose :      Provide objects and methods to work with mutations in an easy way.

Copyright (C): 2016 - Gert Hulselmans
"""

from __future__ import print_function

import numpy
import os.path

import pyfasta

import create_regulatory_domains

fasta_filename = os.path.join(os.path.dirname(__file__),
                              'data',
                              'genomic_fasta',
                              'Homo_sapiens_assembly19_sorted.fasta')


class GenomicFasta:
    fasta_sequences = pyfasta.Fasta(fasta_filename)

    # Calculate chromosome sizes from file index positions for end and start in flattened FASTA file.
    chrom_sizes_dict = {
        chrom: index_pos[1] - index_pos[0] for chrom, index_pos in fasta_sequences.index.iteritems()
    }

    @staticmethod
    def is_chromosome_name(chrom):
        """
        Check if the passed chromosome name is in this assembly.

        :param chrom: Chromosome name.
        :return: True or False
        """
        return chrom in GenomicFasta.chrom_sizes_dict

    @staticmethod
    def chromosome_size(chrom):
        """
        Get chromosome size for the requested chromosome.

        :param chrom: Chromosome name.
        :return: Chromosome size (or 0 if chromsome name was not found).
        """
        return GenomicFasta.chrom_sizes_dict.get(chrom, 0)


class VCFmut:
    # Create a list of GeneTSS objects sorted by chromosome name, TSS, strand and gene name from a TAB-separated file.
    genes_tss_list = create_regulatory_domains.GenesTSSList.load_genes_tss_file(
        genes_tss_filename=create_regulatory_domains.default_genes_tss_filename
    )

    # Calculate the regulatory domains for each gene.
    # See "create_basal_plus_extension_regdoms" for more information.
    reg_doms_list_per_chrom = create_regulatory_domains.create_basal_plus_extension_regdoms(
        genes_tss_list=genes_tss_list,
        maximum_extension=create_regulatory_domains.default_maximum_extension,
        basal_up=create_regulatory_domains.default_basal_up,
        basal_down=create_regulatory_domains.default_basal_down,
        chrom_sizes=create_regulatory_domains.ChromSizes()
    )

    # Store start and end position, tss and strand of each regulatory domain in a per chromosome numpy array.
    reg_doms_start_end_tss_strand_array_per_chrom = {
        chrom: numpy.array(
            [
                [reg_dom.chrom_start, reg_dom.chrom_end, reg_dom.tss, 1 if reg_dom.strand == '+' else -1]
                for reg_dom in reg_doms
            ]
        )
        for chrom, reg_doms in reg_doms_list_per_chrom.iteritems()
    }

    # Store gene name for each regulatory domain in a per chromosome numpy array.
    reg_doms_genes_array_per_chrom = {
        chrom: numpy.array(
            [
                reg_dom.name
                for reg_dom in reg_doms
            ]
        )
        for chrom, reg_doms in reg_doms_list_per_chrom.iteritems()
    }

    @staticmethod
    def from_mut_id(mut_id):
        """
        Create a VCFmut object from a mutation ID

        Examples:
            chr10__100038800__TTTTTG__T__DEL
            chr10__10011659__A__AAT__INS
            chr10__100061062__C__T__SNV

        :param mut_id: mutation ID.
        :return: VCFmut object.
        """

        if mut_id.startswith('#'):
            return None

        try:
            chrom, start, ref, mut = mut_id.split('__')[0:4]
        except ValueError:
            raise ValueError(
                'Mutation ID "{0:s}" is not valid.'.format(
                    mut_id
                )
            )

        return VCFmut(chrom, start, ref, mut)

    @staticmethod
    def from_mut_ids_file(mut_ids_filename):
        """
        Create VCFmut objects from a file which contains mutation IDs.

        :param mut_ids_filename: Filename which contains a list of mutation IDs.
        :return: yield VCFmut objects for each mutation ID in the mut_ids_filename.
        """

        with open(mut_ids_filename, 'r') as mut_ids_fh:
            for line in mut_ids_fh:
                mut_id = line.rstrip('\r\n')

                if mut_id.startswith('#'):
                    continue

                yield VCFmut.from_mut_id(mut_id)

    @staticmethod
    def from_bedlike_mut_id(bedlike_mut_id):
        """
        Create a VCFmut object from a BED-like mutation ID.

        Examples:
            chr10_100038800_100038801_TTTTG_-----_DEL
            chr10_10011659_10011660_--_AT_INS
            chr10_100142677_100142678_T_-_INDEL
            chr10_100061061_100061062_C_T_SNP

        :param bedlike_mut_id: BED-like mutation ID.
        :return: VCFmut object.
        """

        if bedlike_mut_id.startswith('#'):
            return None

        try:
            chrom, _, start, ref, mut, mut_type = bedlike_mut_id.split('_')[0:6]
        except ValueError:
            raise ValueError(
                'BED-like mutation ID "{0:s}" is not valid.'.format(
                    bedlike_mut_id
                )
            )

        if mut_type == 'DEL' or (mut_type == 'INDEL' and mut[0] == '-'):
            # Fix start position and add additional start base for deletion and remove dashes.
            start = int(start) - 1
            ref = 'N' + ref
            mut = 'N'
            ref_at_first_pos = VCFmut(chrom, start, ref, mut).get_reference_sequence_at_vcfmut()[0]
            ref = ref_at_first_pos + ref[1:]
            mut = ref_at_first_pos
        elif mut_type == 'INS' or (mut_type == 'INDEL' and ref[0] == '-'):
            # Fix start position and add additional start base for insertion and remove dashes.
            start = int(start) - 1
            ref = 'N'
            mut = 'N' + mut
            ref_at_first_pos = VCFmut(chrom, start, ref, mut).get_reference_sequence_at_vcfmut()[0]
            ref = ref_at_first_pos
            mut = ref_at_first_pos + mut[1:]

        return VCFmut(chrom, start, ref, mut)

    @staticmethod
    def from_bedlike_mut_ids_file(bedlike_mut_ids_filename):
        """
        Create VCFmut objects from a file which contains BED-like mutation IDs.

        :param bedlike_mut_ids_filename: Filename which contains a list of BED-like mutation IDs.
        :return: yield VCFmut objects for each mutation ID in the bedlike_mut_ids_filename.
        """

        with open(bedlike_mut_ids_filename, 'r') as bedlike_mut_ids_fh:
            for line in bedlike_mut_ids_fh:
                bedlike_mut_id = line.rstrip('\r\n')

                if bedlike_mut_id.startswith('#'):
                    continue

                yield VCFmut.from_bedlike_mut_id(bedlike_mut_id)

    @staticmethod
    def from_vcf_line(vcf_line):
        """
        Create a VCFmut object from a entry from a VCF file (column 1, 2, 4 and 5 are used).

        :param vcf_line: line from VCF file.
        :return: VCFmut object.
        """

        if vcf_line.startswith('#'):
            return None

        try:
            chrom, start, name, ref, mut = vcf_line.rstrip('\r\n').split('\t')[0:5]
        except ValueError:
            raise ValueError(
                'VCF line "{0:s}" needs to contain at least 5 columns.'.format(
                    vcf_line
                )
            )

        return VCFmut(chrom, start, ref, mut)

    @staticmethod
    def from_vcf_file(vcf_filename):
        """
        Create VCFmut objects from a VCF file.

        :param vcf_filename: VCF filename.
        :return: yield a VCFmut object for each mutation in the VCF file.
        """

        with open(vcf_filename, 'r') as vcf_fh:
            for vcf_line in vcf_fh:
                if vcf_line.startswith('#'):
                    continue

                yield VCFmut.from_vcf_line(vcf_line)

    @staticmethod
    def from_fasta_seq_id(fasta_seq_id):
        """
        Create a VCFmut object from a FASTA sequence ID and get the number of base pairs upstream and downstream added
        to the mutation position and get if the FASTA sequence is for the wildtype or mutant sequence.

        The FASTA sequence IDs accepted by this method are made by VCFmut.make_fasta_for_wt_and_mut().

        Examples:
            chr10__10011659__A__AAT__INS__bp_up_10__bp_down_10__wt
            chr10__10011659__A__AAT__INS__bp_up_10__bp_down_10__mut

        :param fasta_seq_id: FASTA sequence ID in the following format: chrom__start__ref__mut__mut_type__bp_up_X__bp_down_Y__mut_or_wt
        :return: VCFmut object, bp_upstream, bp_downstream, is_wt
        """

        chrom, start, ref, mut, mut_type, bp_upstream_str, bp_downstream_str, wt_or_mut_str = fasta_seq_id.split('__')[0:8]

        if (not bp_upstream_str.startswith('bp_up_') or
                not bp_downstream_str.startswith('bp_down_') or
                (wt_or_mut_str != 'wt' and wt_or_mut_str != 'mut')):
            raise ValueError(
                'FASTA ID "{0:s}" is not in the following format: '
                'chrom__start__ref__mut__mut_type__bp_up_X__bp_down_Y__mut_or_wt.'.format(
                    fasta_seq_id
                )
            )

        # Number of base pairs upstream and downstream of the mutation position.
        bp_upstream = int(bp_upstream_str[6:])
        bp_downstream = int(bp_downstream_str[8:])

        # Is this FASTA sequence ID for the wildtype or mutant sequence?
        is_wt = wt_or_mut_str == 'wt'

        return VCFmut(chrom, start, ref, mut), bp_upstream, bp_downstream, is_wt

    def __init__(self, chrom, start, ref, mut, mut_id=None):
        """
        Create a VCFmut object and do some checking to see if a valid mutation was provided.

        :param chrom: Chromosome name on which the mutation is located.
        :param start: Start position of the mutation (1-based coordinate).
        :param ref: Reference sequence for the mutation.
        :param mut: Mutation sequence for the mutation.
        :param mut_id: Mutation ID (if not set, will be created by using the previous arguments.
        :return:
        """

        if mut_id:
            self.mut_line = '{0:s} {1:s} {2:s} {3:s} {4:s}'.format(chrom, str(start), mut_id, ref, mut)
        else:
            self.mut_line = '{0:s} {1:s} {2:s} {3:s}'.format(chrom, str(start), ref, mut)

        if not GenomicFasta.is_chromosome_name(chrom):
            raise ValueError(
                'Chromosome name "{0:s}" is not valid ({1:s}).'.format(chrom, self.mut_line)
            )

        allowed_nucleotides = {'A', 'C', 'G', 'T', 'N'}

        if not set(ref.upper()).issubset(allowed_nucleotides):
            raise ValueError(
                'Reference nucleotide(s) "{0:s}" contain(s) other characters than: A, C, G, T, or N ({1:s}).'.format(
                    ref,
                    self.mut_line
                )
            )

        if not set(mut.upper()).issubset(allowed_nucleotides):
            raise ValueError(
                'Mutation nucleotide(s) "{0:s}" contain(s) other characters than: A, C, G, T, or N ({1:s}).'.format(
                    mut,
                    self.mut_line
                )
            )

        ref_length = len(ref)
        mut_length = len(mut)

        if ref_length == 0:
            raise ValueError('Reference nucleotide can not be empty ({0:s}).'.format(self.mut_line))

        if mut_length == 0:
            raise ValueError('Mutation nucleotide can not be empty ({0:s}).'.format(self.mut_line))

        self.chrom = chrom

        try:
            self.start = int(start)
        except:
            raise ValueError(
                'Mutation position {0:s} is not an integer ({1:s}).'.format(str(start), self.mut_line)
            )

        if self.start == 0:
            raise ValueError(
                'Mutation positions are one-based and thus cannot be zero ({0:s}).'.format(self.mut_line)
            )
        elif self.start > GenomicFasta.chromosome_size(chrom):
            raise ValueError(
                'Mutation position {0:d} is higher than the chromosome length ({1:d}) '
                'for chromosome "{2:s}" ({3:s}).'.format(self.start,
                                                         GenomicFasta.chromosome_size(chrom),
                                                         self.chrom,
                                                         self.mut_line)
            )

        self.ref = ref
        self.mut = mut

        self.snv = False
        self.mnv = False
        self.deletion = False
        self.insertion = False

        if ref_length == mut_length:
            if ref_length == 1:
                # SNV: single nucleotide variant
                self.snv = True
                self.mut_type = 'SNV'
            else:
                # MNV: multi nucleotide variant (two or more SNVs in succession)
                self.mnv = True
                self.mut_type = 'MNV'
        elif ref_length > mut_length:
            # Deletion.
            self.deletion = True
            self.mut_type = 'DEL'
        elif ref_length < mut_length:
            # Insertion.
            self.insertion = True
            self.mut_type = 'INS'

        self.mut_id = mut_id if mut_id else '{0:s}__{1:d}__{2:s}__{3:s}__{4:s}'.format(
            self.chrom,
            self.start,
            self.ref,
            self.mut,
            self.mut_type
        )

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        return '{0:s}\t{1:d}\t{2:s}\t{3:s}\t{4:s}\t{5:s}'.format(self.chrom,
                                                                 self.start,
                                                                 self.ref,
                                                                 self.mut,
                                                                 self.mut_type,
                                                                 self.mut_id)

    def __repr__(self):
        return '<VCFmut>\n' \
               '  chrom: {0:s}\n' \
               '  start: {1:d}\n' \
               '  ref: {2:s}\n' \
               '  mut: {3:s}\n' \
               '  mut_type: {4:s}\n' \
               '  mut_id: {5:s}\n'.format(self.chrom, self.start, self.ref, self.mut, self.mut_type, self.mut_id)

    @property
    def is_snv(self):
        """
        Check if this mutation is a SNV (single nucleotide variant).

        :return: True or False
        """
        return self.snv

    @property
    def is_mnv(self):
        """
        Check if this mutation is a MNV (multi nucleotide variant: two or more SNVs in succession)

        :return: True or False
        """
        return self.mnv

    @property
    def is_deletion(self):
        """
        Check if this mutation is a deletion.

        :return: True or False
        """
        return self.deletion

    @property
    def is_insertion(self):
        """
        Check if this mutation is an insertion.

        :return: True or False
        """
        return self.insertion

    def generate_bed_region(self):
        """
        Generate a BED region for the mutation (length of regions is always 1 base pair).

        :return: list with chromosome name, start position, end position and  mutation ID
        """
        return self.chrom, self.start - 1, self.start, self.mut_id

    def get_associated_genes_and_distance_to_tss(self):
        """
        Get gene names which regulatory domains contain the mutation and the distance of the mutation to the TSS.

        :return: dictionary of gene names as keys and distance to the TSS as values.
        """

        # Create a boolean array for the regulatory regions located on the same chromosome as on which the mutation is
        # located ( = VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom] ) and for which which their
        # interval contains the mutation start position:
        #   - The start position of the regulatory domain of a gene
        #     ( = VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][:, 0] )
        #     is lower than or equal to the mutation position.
        #   - The end position of the regulatory domain of a gene
        #     ( = VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][:, 1] )
        #     is greater than or equal to the mutation position.

        mutation_overlapping_with_reg_doms_array = (
            (VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][:, 0] <= self.start)
            &
            (VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][:, 1] >= self.start)
        )

        # Return a dictionary with:
        #   - Keys: All associated genes.
        #   - Values:
        #       Calculate distance of mutation to TSS of the associated genes:
        #         = (mutation_start_position - tss) * strand
        #       Where:
        #         - mutation_start_position: self.start
        #         - tss: VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2]
        #         - strand: VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 3]
        return dict(
            zip(
                VCFmut.reg_doms_genes_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array].tolist(),
                ((self.start
                  - VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2])
                  * VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 3]
                 ).tolist()
            )
        )

    def get_reference_sequence_at_vcfmut(self):
        """
        Get reference sequence at the mutation position.

        :return: Reference sequence at the mutation position.
        """

        # Calculate one-based closed interval coordinates and make sure chromosome boundaries are not crossed.
        mut_start = self.start
        mut_end = self.start + (len(self.ref) - 1)

        # Reference sequence at the mutation position.
        ref_seq_at_mut_pos = str(
            GenomicFasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': mut_start,
                 'stop': mut_end,
                 'one_based': True}
            )
        )

        return ref_seq_at_mut_pos

    def get_reference_sequences_around_vcfmut(self, bp_upstream, bp_downstream):
        """
        Get reference sequences before, at and after the mutation position.

        :param bp_upstream: Number of base pairs before mutation start.
        :param bp_downstream: Number of base pairs after mutation end.

        :return: ref_seq_before_mut_pos, ref_seq_at_mut_pos, ref_seq_after_mut_pos
        """

        # Calculate one-based closed interval coordinates and make sure chromosome boundaries are not crossed.
        before_mut_start = max(1, self.start - bp_upstream)
        mut_start = self.start
        mut_end = self.start + (len(self.ref) - 1)
        after_mut_end = min(mut_end + bp_downstream,
                            GenomicFasta.chromosome_size(self.chrom))

        # Reference sequence before mutation position.
        ref_seq_before_mut_pos = str(
            GenomicFasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': before_mut_start,
                 'stop': mut_start - 1,
                 'one_based': True}
            )
        )

        # Reference sequence at the mutation position.
        ref_seq_at_mut_pos = str(
            GenomicFasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': mut_start,
                 'stop': mut_end,
                 'one_based': True}
            )
        )

        # Reference sequence before after mutation position.
        ref_seq_after_mut_pos = str(
            GenomicFasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': mut_end + 1,
                 'stop': after_mut_end,
                 'one_based': True}
            )
        )

        return ref_seq_before_mut_pos, ref_seq_at_mut_pos, ref_seq_after_mut_pos

    def make_fasta_for_wt_and_mut(self, bp_upstream, bp_downstream, allow_first_reference_base_to_be_N=False):
        """
        Make FASTA sequence for wildtype and mutation from bp_upstream to bp_downstream of the mutation.

        :param bp_upstream:
                    Number of base pairs before mutation start.
        :param bp_downstream:
                    Number of base pairs after mutation end.
        :param allow_first_reference_base_to_be_N:
                    If set to True, and the mutation is a deletion or insertion,
                    the first nucleotide of the reference and mutation sequence
                    can be "N" and will be replaced by the correct reference
                    nucleotide.

        :return: FASTA string
        """

        (ref_seq_before_mut_pos,
         ref_seq_at_mut_pos,
         ref_seq_after_mut_pos) = self.get_reference_sequences_around_vcfmut(bp_upstream, bp_downstream)

        # Get reference sequence for mutation position from genome.
        genome_ref_for_mut = self.get_reference_sequence_at_vcfmut()

        # Set mutation sequence, so start "N" can be changed to reference sequence if necessary.
        mut = self.mut

        if (allow_first_reference_base_to_be_N and
                    self.ref.upper()[0] == 'N' and
                    self.mut.upper()[0] == 'N' and
                (self.is_deletion or self.is_insertion)):
            # When allow_first_reference_base_to_be_N is set to False, compare only from the
            # reference sequence from the second nucleotide position onwards with the genomic
            # reference nucleotides for deletions or insertions.
            if genome_ref_for_mut.upper()[1:] != self.ref.upper()[1:]:
                raise ValueError(
                    'Provided reference sequence "{0:s}" is not the same as in the genomic reference sequence "{1:s}" '
                    '({2:s}).'.format(self.ref, genome_ref_for_mut, self.mut_line)
                )

            # Change starting "N" with reference nucleotide.
            mut = genome_ref_for_mut[0] + mut[1:]
        else:
            if genome_ref_for_mut.upper() != self.ref.upper():
                raise ValueError(
                    'Provided reference sequence "{0:s}" is not the same as in the genomic reference sequence "{1:s}" '
                    '({2:s}).'.format(self.ref, genome_ref_for_mut, self.mut_line)
                )

        fasta_string = '>{0:s}__bp_up_{1:d}__bp_down_{2:d}__wt\n'.format(self.mut_id, bp_upstream, bp_downstream) \
                       + ref_seq_before_mut_pos \
                       + ref_seq_at_mut_pos \
                       + ref_seq_after_mut_pos \
                       + '\n' \
                       + '>{0:s}__bp_up_{1:d}__bp_down_{2:d}__mut\n'.format(self.mut_id, bp_upstream, bp_downstream) \
                       + ref_seq_before_mut_pos \
                       + mut \
                       + ref_seq_after_mut_pos \
                       + '\n'

        return fasta_string
