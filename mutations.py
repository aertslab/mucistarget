"""
Purpose :      Provide objects and methods to work with mutations in an easy way.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import glob
import numpy
import os.path

import pyfasta

import create_regulatory_domains


fasta_filename_dict = {
    'dm3': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'dm3.fa'
    ),
    'dm6': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'dm6.fa'
    ),
    'hg19': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'hg19.fa'
    ),
    'hg38': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'hg38.fa'
    ),
    'mm9': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'mm9.fa'
    ),
    'mm10': os.path.join(
        os.path.dirname(__file__),
        'data',
        'genomic_fasta',
        'mm10.fa'
    ),
}


genes_tss_filename_dict = {
    'dm3': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'dm3.tss.tsv'
    ),
    'dm6': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'dm6.tss.tsv'
    ),
    'hg19': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'hg19.tss.tsv'
    ),
    'hg38': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'hg38.tss.tsv'
    ),
    'mm9': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'mm9.tss.tsv'
    ),
    'mm10': os.path.join(
        os.path.dirname(__file__),
        'data',
        'regulatory_domains',
        'mm10.tss.tsv'
    ),
}


class GenomicFasta:
    def __init__(self, fasta_filename, assembly):
        self.fasta_sequences = pyfasta.Fasta(fasta_filename)

        # Calculate chromosome sizes from file index positions from starts and ends in flattened FASTA file.
        self.chrom_sizes_dict = {
            chrom: index_pos[1] - index_pos[0] for chrom, index_pos in self.fasta_sequences.index.items()
        }

        self.assembly = assembly

    def is_chromosome_name(self, chrom):
        """
        Check if the passed chromosome name is in this assembly.

        :param chrom: Chromosome name.
        :return: True or False
        """
        return chrom in self.chrom_sizes_dict

    def chromosome_size(self, chrom):
        """
        Get chromosome size for the requested chromosome.

        :param chrom: Chromosome name.
        :return: Chromosome size (or 0 if chromosome name was not found).
        """
        return self.chrom_sizes_dict.get(chrom, 0)


class TADs:
    @staticmethod
    def load_tads_from_file(tads_filename):
        """
        Load TADs from file and store start and end position of each TAD in a per chromosome numpy array.

        :param tads_filename: Filename with TADs.
        :return: Dictionary with chromosome names as keys and a numpy array with start and end coordinates as values.
        """

        tads_chrom_pos_dict = dict()

        with open(tads_filename, 'r') as fh_tads:
            for line in fh_tads:
                columns = line.rstrip('\n').split('\t')

                if len(columns) >= 3:
                    # Set chromosome name as key and start (add 1 as it is in BED format) and end as value.
                    tads_chrom_pos_dict.setdefault(columns[0], []).append([int(columns[1]) + 1, int(columns[2])])

        # Store start and end position of each TAD in a per chromosome numpy array.
        tads_start_end_array_per_chrom = {
            chrom: numpy.array(
                [
                    [region[0], region[1]]
                    for region in regions
                    ]
            )
            for chrom, regions in tads_chrom_pos_dict.items()
            }

        return tads_start_end_array_per_chrom


class VCFmut:
    # Store GenomicFasta.
    genomic_fasta = None

    # Store start and end position, tss and strand of each regulatory domain in a per chromosome numpy array.
    reg_doms_start_end_tss_strand_array_per_chrom = None

    # Store gene name for each regulatory domain in a per chromosome numpy array.
    reg_doms_genes_array_per_chrom = None

    # Load TADs for primary tissues and cell types and store start and end position of each TAD in a per chromosome
    # numpy array.
    tads_start_end_array_per_chrom = None

    @staticmethod
    def set_genome_fasta(fasta_filename, assembly):
        """
        Set genomic fasta file and assembly version for VCFmut class before making VCFmut instances.

        :param fasta_filename: Genomic FASTA filename.
        :param assembly: Assembly version.
        :return:
        """

        VCFmut.genomic_fasta = GenomicFasta(fasta_filename=fasta_filename, assembly=assembly)

    @staticmethod
    def set_reg_doms():
        """
        Set regulatory domains for VCFmut class.

        :return:
        """

        # Create a list of GeneTSS objects sorted by chromosome name,
        # TSS, strand and gene name from a TAB-separated file.
        genes_tss_list = create_regulatory_domains.GenesTSSList.load_genes_tss_file(
            genes_tss_filename=genes_tss_filename_dict[VCFmut.genomic_fasta.assembly]
        )

        # Calculate the regulatory domains for each gene.
        # See "create_basal_plus_extension_regdoms" for more information.
        reg_doms_list_per_chrom = create_regulatory_domains.create_basal_plus_extension_regdoms(
            genes_tss_list=genes_tss_list,
            maximum_extension=create_regulatory_domains.default_maximum_extension,
            basal_up=create_regulatory_domains.default_basal_up,
            basal_down=create_regulatory_domains.default_basal_down,
            chrom_sizes=create_regulatory_domains.ChromSizes(VCFmut.genomic_fasta.chrom_sizes_dict),
            assembly=VCFmut.genomic_fasta.assembly
        )

        # Store start and end position, tss and strand of each regulatory domain in a per chromosome numpy array.
        VCFmut.reg_doms_start_end_tss_strand_array_per_chrom = {
            chrom: numpy.array(
                [
                    [reg_dom.chrom_start, reg_dom.chrom_end, reg_dom.tss, 1 if reg_dom.strand == '+' else -1]
                    for reg_dom in reg_doms
                ]
            )
            for chrom, reg_doms in reg_doms_list_per_chrom.items()
        }

        # Store gene name for each regulatory domain in a per chromosome numpy array.
        VCFmut.reg_doms_genes_array_per_chrom = {
            chrom: numpy.array(
                [
                    reg_dom.name
                    for reg_dom in reg_doms
                ]
            )
            for chrom, reg_doms in reg_doms_list_per_chrom.items()
        }

    @staticmethod
    def set_tads():
        """
        Set TAD domains for human (hg19) VCFmut class.

        :return:
        """

        # Load TADs for primary tissues and cell types and store start and end position of each TAD in a per chromosome
        # numpy array.
        VCFmut.tads_start_end_array_per_chrom = dict()

        tads_filenames_glob = os.path.join(
            os.path.dirname(__file__),
            'data',
            'tads',
            VCFmut.genomic_fasta.assembly,
            'TADs_*.bed'
        )

        # Get all filenames with TADs for the current assembly.
        tads_filenames = sorted(glob.glob(tads_filenames_glob))

        for tads_filename in tads_filenames:
            # Extract TADs identifier from TADs filename.
            tads_id = os.path.basename(tads_filename)[:-4]

            VCFmut.tads_start_end_array_per_chrom[tads_id] \
                = TADs.load_tads_from_file(tads_filename=tads_filename)

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
    def from_mut_ids_file(mut_ids_filename, mut_ids_column=1, return_input_line=False):
        """
        Create VCFmut objects from a file which contains mutation IDs.

        :param mut_ids_filename: Filename which contains a list of mutation IDs.
        :param mut_ids_column: Column number which contains the mutation ID (default=1).
        :param return_input_line: If set to True, also yield the input line for the mutation ID.
        :return: yield VCFmut objects for each mutation ID in the mut_ids_filename, and input line for the mutation if
                 return_input_line is set to True.
        """

        with open(mut_ids_filename, 'r') as mut_ids_fh:
            for line in mut_ids_fh:
                line = line.rstrip('\r\n')

                if line.startswith('#'):
                    continue

                columns = line.split('\t')

                if mut_ids_column > len(columns):
                    raise ValueError(
                        'Mutation line "{0:s}" needs to contain at least {1:d} columns.'.format(
                            line, mut_ids_column
                        )
                    )

                mut_id = columns[mut_ids_column - 1]

                if return_input_line:
                    yield VCFmut.from_mut_id(mut_id), line
                else:
                    yield VCFmut.from_mut_id(mut_id)

    @staticmethod
    def from_zero_based_no_ref_specified(chrom, start, ref, mut):
        """
        Create a VCFmut object from a mutation with zero-based coordinate and without start reference base specified.

        Examples:
            VCFmut.from_zero_based_no_ref_specified('chr10', '100038800', 'TTTTG', '')
            VCFmut.from_zero_based_no_ref_specified('chr10', '100038800', 'TTTTG', '-')
            VCFmut.from_zero_based_no_ref_specified('chr10', '100038800', 'TTTTG', '----')
              ==> mut_id: chr10__100038800__tTTTTG__t__DEL

            VCFmut.from_zero_based_no_ref_specified('chr10', '10011659', '', 'AT')
            VCFmut.from_zero_based_no_ref_specified('chr10', '10011659', '-', 'AT')
            VCFmut.from_zero_based_no_ref_specified('chr10', '10011659', '--', 'AT')
              ==> mut_id: chr10__10011659__A__AAT__INS

            VCFmut.from_zero_based_no_ref_specified('chr10', '100061061', 'C', 'T')
              ==> mut_id: chr10__100061062__C__T__SNV

        :param chrom: Chromosome name on which the mutation is located.
        :param start: Start position of the mutation (zero-based coordinate).
        :param ref: Reference sequence for the mutation (no start reference base for deletions and insertions).
        :param mut: Mutation sequence for the mutation (no start reference base for deletions and insertions).
        :return: VCFmut object.
        """

        try:
            # Fix zero-based start coordinate position to one-based coordinate position.
            start = int(start) + 1
        except ValueError:
            raise ValueError(
                'Mutation position {0:s} is not an integer ({1:s} zero_based_no_ref_specified).'.format(
                    str(start),
                    chrom + '_' + str(start) + '_' + ref + '_' + mut
                )
            )

        # Remove dashes for insertions.
        ref = ref.replace('-', '')
        # Remove dashes for deletions.
        mut = mut.replace('-', '')

        return VCFmut(chrom, start, ref, mut, no_ref_specified=True)

    @staticmethod
    def from_bedlike_mut_id(bedlike_mut_id):
        """
        Create a VCFmut object from a BED-like mutation ID.

        Examples:
            chr10_100038800_100038801_TTTTG_-----_DEL
            chr10_10011659_10011660_--_AT_INS
            chr10_100142677_100142678_T_-_INDEL
            chr10_100061061_100061062_C_T_SNP
            chr10_100038800_100038801_TTTTG_-----
            chr10_100038800_100038801_TTTTG_-
            chr10_10011659_10011660_--_AT
            chr10_10011659_10011660_-_AT
            chr10_100142677_100142678_T_-
            chr10_100061061_100061062_C_T

        :param bedlike_mut_id: BED-like mutation ID.
        :return: VCFmut object.
        """

        if bedlike_mut_id.startswith('#'):
            return None

        try:
            chrom, start, _, ref, mut = bedlike_mut_id.split('_')[0:5]
        except ValueError:
            raise ValueError(
                'BED-like mutation ID "{0:s}" is not valid.'.format(
                    bedlike_mut_id
                )
            )

        return VCFmut.from_zero_based_no_ref_specified(chrom, start, ref, mut)

    @staticmethod
    def from_bedlike_mut_ids_file(bedlike_mut_ids_filename, bedlike_mut_ids_column=1, return_input_line=False):
        """
        Create VCFmut objects from a file which contains BED-like mutation IDs.

        :param bedlike_mut_ids_filename: Filename which contains a list of BED-like mutation IDs.
        :param bedlike_mut_ids_column: Column number which contains the BED-like mutation ID (default=1).
        :param return_input_line: If set to True, also yield the input line for the BED-like mutation ID.
        :return: yield VCFmut objects for each mutation ID in the bedlike_mut_ids_filename, and input line for the
                 mutation if return_input_line is set to True.
        """

        with open(bedlike_mut_ids_filename, 'r') as bedlike_mut_ids_fh:
            for line in bedlike_mut_ids_fh:
                line = line.rstrip('\r\n')

                if line.startswith('#'):
                    continue

                columns = line.split('\t')

                if bedlike_mut_ids_column > len(columns):
                    raise ValueError(
                        'BED-like mutation line "{0:s}" needs to contain at least {1:d} columns.'.format(
                            line, bedlike_mut_ids_column
                        )
                    )

                bedlike_mut_id = columns[bedlike_mut_ids_column - 1]

                if return_input_line:
                    yield VCFmut.from_bedlike_mut_id(bedlike_mut_id), line
                else:
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
    def from_vcf_file(vcf_filename, return_input_line=False):
        """
        Create VCFmut objects from a VCF file.

        :param vcf_filename: VCF filename.
        :param return_input_line: If set to True, also yield the input line for the VCF mutation.
        :return: yield a VCFmut object for each mutation in the VCF file, and input line for the mutation if
                 return_input_line is set to True.
        """

        with open(vcf_filename, 'r') as vcf_fh:
            for vcf_line in vcf_fh:
                if vcf_line.startswith('#'):
                    continue

                if return_input_line:
                    yield VCFmut.from_vcf_line(vcf_line), vcf_line
                else:
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

        :param fasta_seq_id:
            FASTA sequence ID in the following format: chrom__start__ref__mut__mut_type__bp_up_X__bp_down_Y__mut_or_wt
        :return:
            VCFmut object, bp_upstream, bp_downstream, is_wt
        """

        (chrom,
         start,
         ref,
         mut,
         mut_type,
         bp_upstream_str,
         bp_downstream_str,
         wt_or_mut_str) = fasta_seq_id.split('__')[0:8]

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

    def __init__(self, chrom, start, ref, mut, no_ref_specified=False):
        """
        Create a VCFmut object and do some checking to see if a valid mutation was provided.

        :param chrom: Chromosome name on which the mutation is located.
        :param start: Start position of the mutation (1-based coordinate).
        :param ref: Reference sequence for the mutation.
        :param mut: Mutation sequence for the mutation.
        :param no_ref_specified: Set to True, if deletions and insertions do not include
                                 the reference base like specified in VCF format.
        :return:
        """

        self.mut_line = '{0:s} {1:s} {2:s} {3:s}{4:s}'.format(chrom,
                                                              str(start),
                                                              ref,
                                                              mut,
                                                              ' no_ref_specified' if no_ref_specified else '')

        if not VCFmut.genomic_fasta:
            raise ValueError(
                'No genomic FASTA specified: run VCFmut.set_genomic_fasta() first.'
            )

        if not VCFmut.genomic_fasta.is_chromosome_name(chrom):
            raise ValueError(
                'Chromosome name "{0:s}" is not valid ({1:s}).'.format(chrom, self.mut_line)
            )

        self.chrom = chrom

        try:
            self.start = int(start)
        except ValueError:
            raise ValueError(
                'Mutation position {0:s} is not an integer ({1:s}).'.format(str(start), self.mut_line)
            )

        self.snv = False
        self.mnv = False
        self.deletion = False
        self.insertion = False

        ref_length = len(ref)
        mut_length = len(mut)

        if ref_length == mut_length:
            if ref_length == 1:
                # SNV: single nucleotide variant.
                self.snv = True
                self.mut_type = 'SNV'
            else:
                # MNV: multi nucleotide variant (two or more SNVs in succession).
                self.mnv = True
                self.mut_type = 'MNV'
        elif ref_length > mut_length:
            # Deletion.
            self.deletion = True
            self.mut_type = 'DEL'

            if no_ref_specified:
                # Fix start position, if no_ref_specified is True. The reference base will be added later.
                self.start -= 1
        elif ref_length < mut_length:
            # Insertion.
            self.insertion = True
            self.mut_type = 'INS'

            if no_ref_specified:
                # Fix start position, if no_ref_specified is True. The reference base will be added later.
                self.start -= 1

        if self.start <= 0:
            raise ValueError(
                'Mutation positions are one-based and cannot be zero or lower ({0:s}).'.format(self.mut_line)
            )
        elif self.start > VCFmut.genomic_fasta.chromosome_size(chrom):
            raise ValueError(
                'Mutation position {0:d} is higher than the chromosome length ({1:d}) '
                'for chromosome "{2:s}" ({3:s}).'.format(self.start,
                                                         VCFmut.genomic_fasta.chromosome_size(chrom),
                                                         self.chrom,
                                                         self.mut_line)
            )

        if no_ref_specified:
            if self.mut_type == 'DEL' or self.mut_type == 'INS':
                # Add the reference base to the reference and mutation sequence for deletions and insertions
                # if no_ref_specified is True. The start position is already corrected before.
                ref_at_first_pos = VCFmut(chrom, self.start, 'N', 'N').get_reference_sequence_at_vcfmut()
                ref = ref_at_first_pos + ref
                mut = ref_at_first_pos + mut

                # Update length of reference and mutation.
                ref_length = len(ref)
                mut_length = len(mut)

        if ref_length == 0:
            raise ValueError('Reference nucleotide can not be empty ({0:s}).'.format(self.mut_line))

        if mut_length == 0:
            raise ValueError('Mutation nucleotide can not be empty ({0:s}).'.format(self.mut_line))

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

        self.ref = ref
        self.mut = mut

        self.mut_id = '{0:s}__{1:d}__{2:s}__{3:s}__{4:s}'.format(
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
        Check if this mutation is a MNV (multi nucleotide variant: two or more SNVs in succession).

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

        :return: list with chromosome name, start position, end position and  mutation ID.
        """
        return self.chrom, self.start - 1, self.start, self.mut_id

    def get_associated_genes_and_distance_to_tss_and_tss(self):
        """
        Get gene names which regulatory domains contain the mutation and the distance of the mutation to the TSS and TSS.

        :return: dictionary of gene names as keys and distance to the TSS and TSS as values.
        """

        if not VCFmut.reg_doms_start_end_tss_strand_array_per_chrom:
            # If no regulatory domains were loaded, return an empty dictionary.
            return {}

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
        #       - Calculate distance of mutation to TSS of the associated genes:
        #           = (mutation_start_position - tss) * strand
        #         Where:
        #           - mutation_start_position: self.start
        #           - tss: VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2]
        #           - strand: VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 3]
        #       - TSS:
        #           - tss: VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2]
        return dict(
            zip(
                VCFmut.reg_doms_genes_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array].tolist(),
                numpy.append(
                    (self.start
                     - VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2])
                    * VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 3],
                    VCFmut.reg_doms_start_end_tss_strand_array_per_chrom[self.chrom][mutation_overlapping_with_reg_doms_array][:, 2],
                ).reshape([2, -1]).transpose().tolist()
            )
        )

    def tss_of_associated_gene_in_same_tad_as_mutation(self, tads_id, tss):
        """
        Check if mutations falls in TAD of a specific TADs identifier and that the TSS of the associated gene falls in
        the same TAD.

        :param tads_id: Use TADs identifier
        :param tss: TSS of associated gene.

        :return: True or False
        """

        if not VCFmut.tads_start_end_array_per_chrom:
            # If no TADs were loaded, return False.
            return False

        if tads_id not in VCFmut.tads_start_end_array_per_chrom:
            raise ValueError(
                'Unknown TADs identifier "{0:s}".'.format(tads_id)
            )

        if self.chrom not in VCFmut.tads_start_end_array_per_chrom[tads_id]:
            # If the chromosome on which the mutation is located, is not available in the TAD, return False.
            return False

        # Create a boolean array for the TADs located on the same chromosome as on which the mutation is
        # located ( = VCFmut.tads_start_end_array_per_chrom_for_primary_tissues_and_cell_types[primary_tissue_or_cell_type][self.chrom] ) and
        # for which which their interval contains the mutation start position:
        #   - The start position of a TAD
        #     ( = VCFmut.tads_start_end_array_per_chrom_for_primary_tissues_and_cell_types[primary_tissue_or_cell_type][self.chrom][:, 0] )
        #     is lower than or equal to the mutation position.
        #   - The end position of a TAD
        #     ( = VCFmut.tads_start_end_array_per_chrom_for_primary_tissues_and_cell_types[primary_tissue_or_cell_type][self.chrom][:, 1] )
        #     is greater than or equal to the mutation position.

        # If there was any overlap between a TAD and the mutation, return the TAD.
        mutation_in_which_tad = VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][
            (VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][:, 0] <= self.start)
            &
            (VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][:, 1] >= self.start)
        ]

        # If there was any overlap between a TAD and the TSS of the associated gene, return the TAD.
        tss_of_gene_in_which_tad = VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][
            (VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][:, 0] <= tss)
            &
            (VCFmut.tads_start_end_array_per_chrom[tads_id][self.chrom][:, 1] >= tss)
        ]

        if numpy.shape(mutation_in_which_tad) == (1, 2) and numpy.shape(tss_of_gene_in_which_tad) == (1, 2):
            # If both the mutation and the TSS of the associated gene are in a TAD, check if it is the same TAD or not.
            return numpy.all(mutation_in_which_tad == tss_of_gene_in_which_tad)
        elif numpy.shape(mutation_in_which_tad) == (0, 2) or numpy.shape(tss_of_gene_in_which_tad) == (0, 2):
            # If the mutation or the TSS of the associated gene is not in a TAD, return False.
            return False
        else:
            raise ValueError(
                'Mutation "{0:s}" or TSS ({1:d}) overlaps multiple TADs in cell line "{2:s}", which is unsupported.'.format(
                    self.mut_id,
                    tss,
                    tads_id
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
            VCFmut.genomic_fasta.fasta_sequences.sequence(
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
                            VCFmut.genomic_fasta.chromosome_size(self.chrom))

        # Reference sequence before mutation position.
        ref_seq_before_mut_pos = str(
            VCFmut.genomic_fasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': before_mut_start,
                 'stop': mut_start - 1,
                 'one_based': True}
            )
        )

        # Reference sequence at the mutation position.
        ref_seq_at_mut_pos = str(
            VCFmut.genomic_fasta.fasta_sequences.sequence(
                {'chr': self.chrom,
                 'start': mut_start,
                 'stop': mut_end,
                 'one_based': True}
            )
        )

        # Reference sequence before after mutation position.
        ref_seq_after_mut_pos = str(
            VCFmut.genomic_fasta.fasta_sequences.sequence(
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

        if (allow_first_reference_base_to_be_N
                and self.ref.upper()[0] == 'N'
                and self.mut.upper()[0] == 'N'
                and (self.is_deletion or self.is_insertion)):
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
