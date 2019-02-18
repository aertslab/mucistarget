#!/usr/bin/env python3

"""
Purpose:       Create BED file with regulatory domains from file with TSS information for each gene.

Copyright (C): 2016-2019 - Gert Hulselmans
"""

import argparse
import os.path
import sys


# Get genes with TSS for hg19 used by GREAT 3.0:
#   Website: http://bejerano.stanford.edu/help/display/GREAT/Genes
#   Download command:
#     wget -O data/hregulatory_domains/g19.great3.0.genes.txt 'http://bejerano.stanford.edu/help/download/attachments/2752609/hg19.great3.0.genes.txt?version=1&modificationDate=1443465966000&api=v2'
default_genes_tss_filename = os.path.join(os.path.dirname(__file__),
                                          'data',
                                          'regulatory_domains',
                                          'hg19.great3.0.genes.txt')

default_maximum_extension = 1000000
default_basal_up = 5000
default_basal_down = 1000


class ChromSizes:
    """
    Set and get chromosome lengths.
    """

    def __init__(self, chrom_sizes_dict=None):
        """
        Define chromosome sizes.

        :param chrom_sizes_dict: dictionary with chromosome names as key and chromosome lengths as values.

        :return: ChromSizes object
        """
        if chrom_sizes_dict:
            self.chrom_sizes_dict = chrom_sizes_dict
        else:
            # If no chromosome sizes where specified, use the following (chromosome sizes for hg19).
            self.chrom_sizes_dict = {
                'chr1': 249250621,
                'chr10': 135534747,
                'chr11': 135006516,
                'chr11_gl000202_random': 40103,
                'chr12': 133851895,
                'chr13': 115169878,
                'chr14': 107349540,
                'chr15': 102531392,
                'chr16': 90354753,
                'chr17': 81195210,
                'chr17_ctg5_hap1': 1680828,
                'chr17_gl000203_random': 37498,
                'chr17_gl000204_random': 81310,
                'chr17_gl000205_random': 174588,
                'chr17_gl000206_random': 41001,
                'chr18': 78077248,
                'chr18_gl000207_random': 4262,
                'chr19': 59128983,
                'chr19_gl000208_random': 92689,
                'chr19_gl000209_random': 159169,
                'chr1_gl000191_random': 106433,
                'chr1_gl000192_random': 547496,
                'chr2': 243199373,
                'chr20': 63025520,
                'chr21': 48129895,
                'chr21_gl000210_random': 27682,
                'chr22': 51304566,
                'chr3': 198022430,
                'chr4': 191154276,
                'chr4_ctg9_hap1': 590426,
                'chr4_gl000193_random': 189789,
                'chr4_gl000194_random': 191469,
                'chr5': 180915260,
                'chr6': 171115067,
                'chr6_apd_hap1': 4622290,
                'chr6_cox_hap2': 4795371,
                'chr6_dbb_hap3': 4610396,
                'chr6_mann_hap4': 4683263,
                'chr6_mcf_hap5': 4833398,
                'chr6_qbl_hap6': 4611984,
                'chr6_ssto_hap7': 4928567,
                'chr7': 159138663,
                'chr7_gl000195_random': 182896,
                'chr8': 146364022,
                'chr8_gl000196_random': 38914,
                'chr8_gl000197_random': 37175,
                'chr9': 141213431,
                'chr9_gl000198_random': 90085,
                'chr9_gl000199_random': 169874,
                'chr9_gl000200_random': 187035,
                'chr9_gl000201_random': 36148,
                'chrM': 16571,
                'chrUn_gl000211': 166566,
                'chrUn_gl000212': 186858,
                'chrUn_gl000213': 164239,
                'chrUn_gl000214': 137718,
                'chrUn_gl000215': 172545,
                'chrUn_gl000216': 172294,
                'chrUn_gl000217': 172149,
                'chrUn_gl000218': 161147,
                'chrUn_gl000219': 179198,
                'chrUn_gl000220': 161802,
                'chrUn_gl000221': 155397,
                'chrUn_gl000222': 186861,
                'chrUn_gl000223': 180455,
                'chrUn_gl000224': 179693,
                'chrUn_gl000225': 211173,
                'chrUn_gl000226': 15008,
                'chrUn_gl000227': 128374,
                'chrUn_gl000228': 129120,
                'chrUn_gl000229': 19913,
                'chrUn_gl000230': 43691,
                'chrUn_gl000231': 27386,
                'chrUn_gl000232': 40652,
                'chrUn_gl000233': 45941,
                'chrUn_gl000234': 40531,
                'chrUn_gl000235': 34474,
                'chrUn_gl000236': 41934,
                'chrUn_gl000237': 45867,
                'chrUn_gl000238': 39939,
                'chrUn_gl000239': 33824,
                'chrUn_gl000240': 41933,
                'chrUn_gl000241': 42152,
                'chrUn_gl000242': 43523,
                'chrUn_gl000243': 43341,
                'chrUn_gl000244': 39929,
                'chrUn_gl000245': 36651,
                'chrUn_gl000246': 38154,
                'chrUn_gl000247': 36422,
                'chrUn_gl000248': 39786,
                'chrUn_gl000249': 38502,
                'chrX': 155270560,
                'chrY': 59373566
            }

    @staticmethod
    def load_chrom_sizes_file(chrom_sizes_filename):
        """
        Create a ChromSizes object by reading chromosome sizes from a chromosome sizes file.

        Example of a chromosome sizes file:
            # Chromosome name    Chromosome length
            chr1    249250621
            chr2    243199373
            chr3    198022430
            chr4    191154276
            chr5    180915260
        """
        chrom_sizes_dict = dict()

        with open(chrom_sizes_filename, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                columns = line.rstrip('\n').split()
                if len(columns) == 2:
                    # Save chromosome name and chromsome size in dictionary.
                    chrom_sizes_dict[columns[0]] = int(columns[1])

        return ChromSizes(chrom_sizes_dict)

    def get_chrom_size_for_chrom(self, chrom):
        """
        Get chromsome length for specified chromosome name.

        :param chrom: chromosome name for which you want to know the length.

        :return: chromosome length (or None if chromosome name was not found).
        """
        return self.chrom_sizes_dict.get(chrom, None)


class GeneTSS:
    """
    Create a GeneTSS object.
    """

    def __init__(self, chrom, tss, strand, name):
        """
        Create a GeneTSS object.

        :param chrom: chromosome name
        :param tss: TSS
        :param strand: strand
        :param name: gene name

        :return: GeneTSS object
        """
        self.chrom = chrom
        self.tss = int(tss)
        self.strand = strand
        self.name = name

    def __lt__(self, other):
        """
        When sorting different GeneTSS objects, sort them by:
           - chromosome name
           - TSS
           - strand
           - gene name
        """
        if self.chrom == other.chrom:
            # It the chromosome is the same, check the TSS.
            if self.tss == other.tss:
                # If the TSS is the same check the strand.
                if self.strand == other.strand:
                    # If the stand is the same, check the gene name.
                    return self.name < other.name
                elif self.strand == '+':
                    # If the strand was different, take the one with the positive strand first.
                    return True
                else:
                    return False
            else:
                # If the TSS is different, take the smallest one first.
                return self.tss < other.tss
        else:
            # If the chromsome name is not the same, check which one is comes first.
            return self.chrom < other.chrom


class GenesTSSList:
    """
    Create a list of GeneTSS objects sorted by chromosome name, TSS, strand and gene name.
    """

    def __init__(self, genes_tss_list):
        """
        Create a list of GeneTSS objects sorted by chromosome name, TSS, strand and gene name.

        :param genes_tss_list: list of GeneTSS objects.

        :return: sorted (by chromosome, TSS, strand and gene name.) list of GeneTSS objects.
        """
        self.genes_tss_list = sorted(genes_tss_list)

    @staticmethod
    def load_genes_tss_file(genes_tss_filename):
        """
        Create a list of GeneTSS objects sorted by chromosome name, TSS, strand and gene name
        from a TAB-separated file in one of the following formats:
          - Chromosome name, TSS, strand, gene name.
          - ENSEMBL gene ID, chromosome name, TSS, strand, gene name.

        :param genes_tss_filename: Filename with TSS info for each gene.

        :return: GenesTSSList object which contains a sorted list of GeneTSS objects
                 by chromosome name, TSS, strand and gene name.
        """
        genes_tss_list = list()

        with sys.stdin if genes_tss_filename in ('-', 'stdin') else open(genes_tss_filename, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                columns = line.rstrip('\n').split('\t')
                if len(columns) == 4:
                    # Chromosome name, TSS, strand, gene name.
                    genes_tss_list.append(GeneTSS(columns[0], columns[1], columns[2], columns[3]))
                elif len(columns) == 5:
                    # ENSEMBL gene ID, chromosome name, TSS, strand, gene name.
                    genes_tss_list.append(GeneTSS(columns[1], columns[2], columns[3], columns[4]))

        return GenesTSSList(genes_tss_list)

    def __iter__(self):
        return self.genes_tss_list.__iter__()


class RegDom:
    """
    Create a regulatory domain.
    """

    def __init__(self, chrom, chrom_start, chrom_end, name, tss, strand, basal_up=5000, basal_down=1000,
                 chrom_sizes=ChromSizes()):
        """
        Create a regulatory domain.

        :param chrom: chromosome name
        :param chrom_start: chromosome start
        :param chrom_end: chromosome end
        :param name: gene name
        :param tss: TSS
        :param strand: strand
        :param basal_up: number of bases upstream of TSS to create a basal domain region.
        :param basal_down: number of bases downstream of TSS to create a basal domain region.
        :param chrom_sizes: ChromSizes object with chromosome lengths for each chromosome

        :return: regulatory domain object with:
                   - chrom: chromosome name on which the regulatory domain is located
                   - chrom_start: regulatory domain start
                   - chrom_end: regulatory domain end
                   - name: gene name
                   - tss: TSS
                   - strand: strand
                   - basal_start: start of basal domain region
                   - basal_end: end of basal domain region

        If basal_up and basal_down are both None:
          - basal_start = chrom_start
          - basal_end = chrom_end.
        This is useful when creating a RegDom for a curated regulatory domain.

        If chrom_start and chrom_end are both None:
          - chrom_start = basal_start
          - chrom_end = basal_end.
        This is useful when creating a RegDom for genes for which you want to
        create a basal domain based on TSS and basal_up and basal_down values.
        """
        self.chrom = chrom
        self.chrom_start = chrom_start
        self.chrom_end = chrom_end
        self.name = name
        self.tss = tss
        self.strand = strand

        if basal_up and basal_down:
            if self.strand == '+':
                self.basal_start = max(0, self.tss - basal_up)
                self.basal_end = min(chrom_sizes.get_chrom_size_for_chrom(self.chrom), self.tss + basal_down)
            elif self.strand == '-':
                self.basal_start = max(0, self.tss - basal_down)
                self.basal_end = min(chrom_sizes.get_chrom_size_for_chrom(self.chrom), self.tss + basal_up)

            if self.chrom_start is None:
                self.chrom_start = self.basal_start

            if self.chrom_end is None:
                self.chrom_end = self.basal_end
        elif self.chrom_start and self.chrom_end:
            self.basal_start = self.chrom_start
            self.basal_end = self.chrom_end

    def __str__(self):
        return '{0:s}\t{1:d}\t{2:d}\t{3:s}\t{4:d}\t{5:s}\t{6:d}\t{7:d}'.format(
            self.chrom,
            self.chrom_start,
            self.chrom_end,
            self.name,
            self.tss,
            self.strand,
            self.basal_start,
            self.basal_end)


class CuratedRegDoms:
    """
    List of curated regulatory domains for hg19.
    """
    def __init__(self, chrom_sizes=ChromSizes()):
        """
        Create a list of curated regulatory domains for hg19:
          http://bejerano.stanford.edu/help/display/GREAT/Association+Rules#AssociationRules-CuratedRegulatoryDomains

            Sonic hedgehog long-range enhancer:
                SHH: chr7:155438203-156584569

            HOXD global control region:
                KIAA1715 (= LNP): chr2:176714855-176947690
                EVX2: chr2:176714855-176953690
                HOXD13: chr2:176714855-176959529
                HOXD12: chr2:176714855-176967083
                HOXD11: chr2:176714855-176976491
                HOXD10: chr2:176714855-176982491

            Beta-globin locus control region:
                HBB: chr11:5226931-5314124
                HBD: chr11:5253302-5314124
                HBG1: chr11:5260859-5314124
                HBE1: chr11:5276088-5314124

        :param chrom_sizes: ChromSizes object with chromosome lengths.

        :return: CuratedRegDoms object
        """

        self.curated_reg_doms_dict = dict()

        # Sonic hedgehog long-range enhancer.
        self.curated_reg_doms_dict['SHH'] = RegDom(chrom='chr7',
                                                   chrom_start=155438203,
                                                   chrom_end=156584569,
                                                   name='SHH',
                                                   tss=155604967,
                                                   strand='-',
                                                   basal_up=None,
                                                   basal_down=None,
                                                   chrom_sizes=chrom_sizes)

        # HOXD global control region.
        self.curated_reg_doms_dict['KIAA1715'] = RegDom(chrom='chr2',
                                                        chrom_start=176714855,
                                                        chrom_end=176947690,
                                                        name='KIAA1715',
                                                        tss=176867073,
                                                        strand='-',
                                                        basal_up=None,
                                                        basal_down=None,
                                                        chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['EVX2'] = RegDom(chrom='chr2',
                                                    chrom_start=176714855,
                                                    chrom_end=176953690,
                                                    name='EVX2',
                                                    tss=176948641,
                                                    strand='-',
                                                    basal_up=None,
                                                    basal_down=None,
                                                    chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HOXD13'] = RegDom(chrom='chr2',
                                                      chrom_start=176714855,
                                                      chrom_end=176959529,
                                                      name='HOXD13',
                                                      tss=176957618,
                                                      strand='+',
                                                      basal_up=None,
                                                      basal_down=None,
                                                      chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HOXD12'] = RegDom(chrom='chr2',
                                                      chrom_start=176714855,
                                                      chrom_end=176967083,
                                                      name='HOXD12',
                                                      tss=176964457,
                                                      strand='+',
                                                      basal_up=None,
                                                      basal_down=None,
                                                      chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HOXD11'] = RegDom(chrom='chr2',
                                                      chrom_start=176714855,
                                                      chrom_end=176976491,
                                                      name='HOXD11',
                                                      tss=176972013,
                                                      strand='+',
                                                      basal_up=None,
                                                      basal_down=None,
                                                      chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HOXD10'] = RegDom(chrom='chr2',
                                                      chrom_start=176714855,
                                                      chrom_end=176982491,
                                                      name='HOXD10',
                                                      tss=176981306,
                                                      strand='+',
                                                      basal_up=None,
                                                      basal_down=None,
                                                      chrom_sizes=chrom_sizes)

        # Beta-globin locus control region.
        self.curated_reg_doms_dict['HBB'] = RegDom(chrom='chr11',
                                                   chrom_start=5226931,
                                                   chrom_end=5314124,
                                                   name='HBB',
                                                   tss=5248427,
                                                   strand='-',
                                                   basal_up=None,
                                                   basal_down=None,
                                                   chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HBD'] = RegDom(chrom='chr11',
                                                   chrom_start=5253302,
                                                   chrom_end=5314124,
                                                   name='HBD',
                                                   tss=5255878,
                                                   strand='-',
                                                   basal_up=None,
                                                   basal_down=None,
                                                   chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HBG1'] = RegDom(chrom='chr11',
                                                    chrom_start=5260859,
                                                    chrom_end=5314124,
                                                    name='HBG1',
                                                    tss=5271122,
                                                    strand='-',
                                                    basal_up=None,
                                                    basal_down=None,
                                                    chrom_sizes=chrom_sizes)
        self.curated_reg_doms_dict['HBE1'] = RegDom(chrom='chr11',
                                                    chrom_start=5276088,
                                                    chrom_end=5314124,
                                                    name='HBE1',
                                                    tss=5526834,
                                                    strand='-',
                                                    basal_up=None,
                                                    basal_down=None,
                                                    chrom_sizes=chrom_sizes)

    def get_curated_reg_doms_for_gene(self, gene):
        """
        Get the curated regulatory domain for a gene.

        :param gene: gene name

        :return: RegDom object if gene name has a curated regulatory domain or None if it does not have one.
        """
        return self.curated_reg_doms_dict.get(gene, None)


def create_basal_plus_extension_regdoms(genes_tss_list,
                                        maximum_extension=1000000,
                                        basal_up=5000,
                                        basal_down=1000,
                                        chrom_sizes=ChromSizes()):
    """
    Create regulatory domains from a sorted GenesTSSList object.

    The creation of regulatory domains is based on the Basal plus extension approach of GREAT:
      http://bejerano.stanford.edu/help/display/GREAT/Association+Rules#AssociationRules-Approach1Basalplusextension

    Each gene is assigned a basal regulatory domain of a minimum distance
    upstream and downstream of the TSS (regardless of other nearby genes).

    The gene regulatory domain is extended in both directions to the nearest
    gene's basal domain but no more than a maximum extension in one direction.

    When extending the regulatory domain of gene G beyond its basal domain,
    the extension to the "left" extends until it reaches the first basal domain
    of any gene whose transcription start site is "left" of G's transcription
    start site (and analogously for extending "right").

    When the TSS start site and strand are the same, both genes will have the
    same regulatory domain (regulatory domain size is then determined by the
    next closest gene).

    :param genes_tss_list:
                GenesTSSList object with all genes which you want to consider
                for making regulatory domains.
    :param maximum_extension:
                maximum extension size in base pairs a regulatory domain can be
                extended if it did not encounter a basal domain of the nearest gene.
    :param basal_up:
                # bp upstream of TSS for setting the basal start of the basal domain.
    :param basal_down:
                # bp downstream of TSS for setting the basal end of the basal domain.
    :param chrom_sizes:
                ChromSizes object with chromosome lengths.

    :return: List of RegDom objects per chromosome.
    """

    # Store all regulatory domains in a per chromosome list.
    reg_doms_list_per_chrom = {chrom: list() for chrom in chrom_sizes.chrom_sizes_dict}

    # Store basal regulatory domain starts and ends per chromosome.
    basal_starts_per_chrom_dict = {chrom: list() for chrom in chrom_sizes.chrom_sizes_dict}
    basal_ends_per_chrom_dict = {chrom: list() for chrom in chrom_sizes.chrom_sizes_dict}

    # Store gene index for each gene in basal_starts_per_chrom_dict and basal_ends_per_chrom_dict.
    gene_idx_in_per_chrom_dicts_dict = dict()

    # List of curated regulatory domains to add at the end of the reg_doms_list_per_chrom.
    curated_reg_doms_list_to_add_per_chrom = {chrom: list() for chrom in chrom_sizes.chrom_sizes_dict}

    # Get all curated regulatory domains.
    curated_reg_doms = CuratedRegDoms()

    # Loop over GenesTSSList.
    for curr_gene_tss in genes_tss_list:
        if curated_reg_doms.get_curated_reg_doms_for_gene(curr_gene_tss.name):
            # If the current GeneTSS has a curated regulatory domain, do not add them to the reg_doms_list_per_chrom
            # yet. They will be added after all genes with non-curated regulatory domains are added.
            curr_regdom = curated_reg_doms.get_curated_reg_doms_for_gene(curr_gene_tss.name)
            curated_reg_doms_list_to_add_per_chrom[curr_gene_tss.chrom].append(curr_regdom)
        else:
            # Create a regulatory domain based on the current GeneTSS.
            curr_regdom = RegDom(chrom=curr_gene_tss.chrom,
                                 chrom_start=None,
                                 chrom_end=None,
                                 name=curr_gene_tss.name,
                                 tss=curr_gene_tss.tss,
                                 strand=curr_gene_tss.strand,
                                 basal_up=basal_up,
                                 basal_down=basal_down,
                                 chrom_sizes=chrom_sizes)

            reg_doms_list_per_chrom[curr_gene_tss.chrom].append(curr_regdom)

        # Store basal regulatory domain start and end for the current regulatory domain per chromosome.
        basal_starts_per_chrom_dict[curr_gene_tss.chrom].append(curr_regdom.basal_start)
        basal_ends_per_chrom_dict[curr_gene_tss.chrom].append(curr_regdom.basal_end)

        # Store the gene index for the current GeneTSS in the per chromosome list
        # in basal_starts_per_chrom_dict and basal_ends_per_chrom_dict, so later
        # we will be able to get all basal start locations before a certain TSS of
        # a gene and all basal end locations of a certain TSS, so we do not extend
        # a regulatory domain to far.
        gene_idx_in_per_chrom_dicts_dict[
            (curr_gene_tss.name,
             curr_gene_tss.chrom,
             curr_gene_tss.tss,
             curr_gene_tss.strand)] = len(basal_starts_per_chrom_dict[curr_gene_tss.chrom]) - 1

    # Loop over all regulatory domains per chromosome.
    for chrom in reg_doms_list_per_chrom:
        for idx, curr_regdom in enumerate(reg_doms_list_per_chrom[chrom]):
            # Get the previous regulatory domain if this is not the first regulatory domain of the list.
            prev_regdom = (reg_doms_list_per_chrom[chrom][idx - 1]
                           if idx > 0
                           else None)

            # Get the next regulatory domain if this is not the last regulatory domain of the list.
            next_regdom = (reg_doms_list_per_chrom[chrom][idx + 1]
                           if idx < (len(reg_doms_list_per_chrom[chrom]) - 1)
                           else None)

            # Extend regulatory domain as much as possible to the left.
            tmp_start = min(curr_regdom.basal_start,
                            max(0, curr_regdom.tss - maximum_extension))

            if prev_regdom:
                if (prev_regdom.tss, prev_regdom.strand) == (curr_regdom.tss, curr_regdom.strand):
                    # The previous regulatory domain had the same TSS and strand.
                    tmp_start = min(curr_regdom.basal_start, prev_regdom.chrom_start)
                else:
                    # Get the gene index of the current gene in the per chromosome list of the basal_ends_per_chrom_dict
                    # and substract one.
                    prev_gene_idx_in_per_chrom_dicts = gene_idx_in_per_chrom_dicts_dict[(curr_regdom.name,
                                                                                         curr_regdom.chrom,
                                                                                         curr_regdom.tss,
                                                                                         curr_regdom.strand)] - 1

                    if prev_gene_idx_in_per_chrom_dicts >= 1:
                        # This is the third or higher gene on this chromosome.

                        # Adapt the start site of the current regulatory domain by taking the lowest genomic coordinate
                        # of:
                        #   - the current regulatory basal domain start.
                        #   - the highest genomic coordinate of:
                        #       - the previous regulatory basal domain end.
                        #       - the extend regulatory domain as much as possible to the left.
                        #       - the highest basal end of all previous genes on the current chromosome
                        #         (which is not always the previous gene, as strandness plays a role in the basal domain
                        #         creation).
                        tmp_start = min(curr_regdom.basal_start,
                                        max(prev_regdom.basal_end,
                                            tmp_start,
                                            max(basal_ends_per_chrom_dict[curr_regdom.chrom][
                                                0:prev_gene_idx_in_per_chrom_dicts
                                                ])
                                            )
                                        )
                    else:
                        # This is the second gene on this chromosome.
                        tmp_start = min(curr_regdom.basal_start,
                                        max(prev_regdom.basal_end, tmp_start))

            # Extend regulatory domain as much as possible to the right.
            tmp_end = max(curr_regdom.basal_end,
                          min(chrom_sizes.get_chrom_size_for_chrom(curr_regdom.chrom),
                              curr_regdom.tss + maximum_extension)
                          )

            if next_regdom:
                # The next regulatory domain was on the same chromosome as the current one.

                # Get the gene index of the current gene in the per chromosome list of the basal_starts_per_chrom_dict
                # and add one.
                next_gene_idx_in_per_chrom_dicts = gene_idx_in_per_chrom_dicts_dict[(curr_regdom.name,
                                                                                     curr_regdom.chrom,
                                                                                     curr_regdom.tss,
                                                                                     curr_regdom.strand)] + 1

                # Adapt the end site of the current regulatory domain by taking the highest genomic coordinate of:
                #   - the current regulatory basal domain end.
                #   - the lowest genomic coordinate of:
                #       - the next regulatory basal domain start.
                #       - the extend regulatory domain as much as possible to the right.
                #       - the lowest basal start of all next genes on the current chromosome
                #         (which is not always the next gene, as strandness plays a role in the basal domain creation).
                tmp_end = max(curr_regdom.basal_end,
                              min(next_regdom.basal_start,
                                  tmp_end,
                                  min(basal_starts_per_chrom_dict[curr_regdom.chrom][next_gene_idx_in_per_chrom_dicts:])
                                  )
                              )

            # Update the regulatory domain start and end in the reg_doms_list_per_chrom.
            reg_doms_list_per_chrom[chrom][idx].chrom_start = tmp_start
            reg_doms_list_per_chrom[chrom][idx].chrom_end = tmp_end

            if prev_regdom and (
                        (prev_regdom.tss, prev_regdom.strand) ==
                        (curr_regdom.tss, curr_regdom.strand)
            ):
                # The previous regulatory domain has the same chromosome name, TSS and strand.

                prev_idx = idx - 1

                # Fix previous regulatory domain starts and current regulatory domain start as long as
                # TSS and strand are the same as for the current regulatory domain.
                while prev_idx >= 0 and (
                            (reg_doms_list_per_chrom[chrom][prev_idx].tss,
                             reg_doms_list_per_chrom[chrom][prev_idx].strand) ==
                            (curr_regdom.tss,
                             curr_regdom.strand)
                ):
                    reg_doms_list_per_chrom[chrom][prev_idx].chrom_end = max(
                        reg_doms_list_per_chrom[chrom][prev_idx].chrom_end,
                        curr_regdom.chrom_end
                    )
                    reg_doms_list_per_chrom[chrom][idx].chrom_start = \
                        reg_doms_list_per_chrom[chrom][prev_idx].chrom_start

                    prev_idx -= 1

    # Add curated regulatory domains at the end.
    for chrom in curated_reg_doms_list_to_add_per_chrom:
        reg_doms_list_per_chrom[chrom].extend(curated_reg_doms_list_to_add_per_chrom[chrom])

    return reg_doms_list_per_chrom


def main():
    parser = argparse.ArgumentParser(
        description='Create BED file with regulatory domains from file with TSS information for each gene.'
    )

    parser.add_argument('--basal-up',
                        dest='basal_up',
                        action='store',
                        type=int,
                        required=False,
                        default=default_basal_up,
                        help='number of bases upstream of TSS to create a basal domain region '
                             '(default: "{0:d}").'.format(default_basal_up)
                        )
    parser.add_argument('--basal-down',
                        dest='basal_down',
                        action='store',
                        type=int,
                        required=False,
                        default=default_basal_down,
                        help='number of bases downstream of TSS to create a basal domain region '
                             '(default: "{0:d}").'.format(default_basal_down)
                        )
    parser.add_argument('--max-ext',
                        dest='maximum_extension',
                        action='store',
                        type=int,
                        required=False,
                        default=default_maximum_extension,
                        help='maximum extension size in base pairs a regulatory domain can be '
                             'extended if it did not encounter a basal domain of the nearest gene. '
                             '(default: "{0:d}").'.format(default_maximum_extension)
                        )
    parser.add_argument('--genes-tss',
                        dest='genes_tss_filename',
                        action='store',
                        type=str,
                        required=False,
                        default=default_genes_tss_filename,
                        help='TAB-separated file with 4 (chromosome name, TSS, strand, gene name) or '
                             '5 (ENSEMBL gene ID, chromosome name, TSS, strand, gene name) columns '
                             'which is used as input to create the regulatory domains. '
                             '(default: "{0:s}").'.format(default_genes_tss_filename)
                        )
    parser.add_argument('--regdoms',
                        dest='regulatory_domains_bed_filename',
                        action='store',
                        type=str,
                        required=False,
                        default='-',
                        help='Regulatory domains BED output file (default: "stdout").'
                        )

    args = parser.parse_args()

    # Create a list of GeneTSS objects sorted by chromosome name, TSS, strand and gene name from a TAB-separated file.
    genes_tss_list = GenesTSSList.load_genes_tss_file(
        genes_tss_filename=args.genes_tss_filename
    )

    # Calculate the regulatory domains for each gene.
    # See "create_basal_plus_extension_regdoms" for more information.
    reg_doms_list_per_chrom = create_basal_plus_extension_regdoms(genes_tss_list=genes_tss_list,
                                                                  maximum_extension=args.maximum_extension,
                                                                  basal_up=args.basal_up,
                                                                  basal_down=args.basal_down,
                                                                  chrom_sizes=ChromSizes())

    if args.regulatory_domains_bed_filename in ('-', 'stdout'):
        regulatory_domains_output_bed_fh = sys.stdout
    else:
        regulatory_domains_output_bed_fh = open(args.regulatory_domains_bed_filename, 'w')

    # Print the regulatory domain output to standard output.
    print('# chrom',
          'chrom_start',
          'chrom_end',
          'name',
          'tss',
          'strand',
          'basal_start',
          'basal_end',
          sep='\t',
          file=regulatory_domains_output_bed_fh)

    for chrom in sorted(reg_doms_list_per_chrom):
        for reg_dom in reg_doms_list_per_chrom[chrom]:
            print(reg_dom, file=regulatory_domains_output_bed_fh)

    regulatory_domains_output_bed_fh.close()


if __name__ == "__main__":
    main()
