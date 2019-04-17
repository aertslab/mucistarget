# Regulatory domains

## Download knownCanonical transcription start sites from UCSC

Get knownCanonical transcription start sites from [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables)
for usage with `create_regulatory_domains.py`, which reimplements *Basal plus extension* section of
*Associating genomic regions with genes* of [GREAT](http://bejerano.stanford.edu/great/public/html/).


### Human (hg19)

Go to [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables) and select the following:

|     Parameter     |                       Value                       |
| ----------------: | :-----------------------------------------------: |
|         **clade** | *Mammal*                                          |
|        **genome** | *Human*                                           |
|      **assembly** | *Feb. 2009 (GRCh37/hg19)*                         |
|         **group** | *Genes and Gene Predictions*                      |
|         **track** | *UCSC Genes*                                      |
|         **table** | *knownCanonical*                                  |
| **output format** | *selected fields from primary and related tables* |
|   **output file** | `hg19.knownCanonical.full_table.tsv`              |

Press **get output**.


Select Fields from **hg19.knownCanonical**:

|     Field     |                          Description                          |
| ------------: | ------------------------------------------------------------- |
| **clusterId** | Which cluster of transcripts this belongs to in knownIsoforms |


Select the following **Linked Tables**:

|    go    |   goaPart   |                   Description                    |
| -------: | :---------: | ------------------------------------------------ |
| **hg19** | *kgXref*    | Link together a Known Gene ID and a gene alias   |
| **hg19** | *knownGene* | Transcript from default gene set in UCSC browser |

Press **allow selection from checked tables**.


Select the following fields in **hg19.kgXref fields**:

|     Field      | Description |
| -------------: | :---------: |
| **geneSymbol** | Gene Symbol |
|     **refseq** | RefSeq ID   |


Select the following fields in **hg19.knownGene fields**:

|     Field      |                             Description                              |
| -------------: | -------------------------------------------------------------------- |
|       **name** | Name of gene                                                         |
|      **chrom** | Reference sequence chromosome or scaffold                            |
|     **strand** | + or - for strand                                                    |
|    **txStart** | Transcription start position (or end position for minus strand item) |
|      **txEnd** | Transcription end position (or start position for minus strand item) |
|   **cdsStart** | Coding region start (or end position if for minus strand item)       |
|     **cdsEnd** | Coding region end (or start position if for minus strand item)       |
|  **exonCount** | Number of exons                                                      |
| **exonStarts** | Exon start positions (or end positions for minus strand item)        |
|   **exonEnds** | Exon end positions (or start positions for minus strand item)        |
|  **proteinID** | UniProt display ID, UniProt accession, or RefSeq protein ID          |

Press **get output**.



### Human (hg38)

Go to [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables) and select the following:

|     Parameter     |                       Value                       |
| ----------------: | :-----------------------------------------------: |
|         **clade** | *Mammal*                                          |
|        **genome** | *Human*                                           |
|      **assembly** | *Dec. 2013 (GRCh38/hg38)*                         |
|         **group** | *Genes and Gene Predictions*                      |
|         **track** | *GENCODE v29*                                     |
|         **table** | *knownCanonical*                                  |
| **output format** | *selected fields from primary and related tables* |
|   **output file** | `hg38.knownCanonical.full_table.tsv`              |

Press **get output**.


Select Fields from **hg38.knownCanonical**:

|     Field     |                          Description                          |
| ------------: | ------------------------------------------------------------- |
| **clusterId** | Which cluster of transcripts this belongs to in knownIsoforms |


Select the following **Linked Tables**:

|    go    |   goaPart   |                   Description                    |
| -------: | :---------: | ------------------------------------------------ |
| **hg38** | *kgXref*    | Link together a Known Gene ID and a gene alias   |
| **hg38** | *knownGene* | Transcript from default gene set in UCSC browser |

Press **allow selection from checked tables**.


Select the following fields in **hg38.kgXref fields**:

|     Field      | Description |
| -------------: | :---------: |
| **geneSymbol** | Gene Symbol |
|     **refseq** | RefSeq ID   |


Select the following fields in **hg38.knownGene fields**:

|     Field      |                             Description                              |
| -------------: | -------------------------------------------------------------------- |
|       **name** | Name of gene                                                         |
|      **chrom** | Reference sequence chromosome or scaffold                            |
|     **strand** | + or - for strand                                                    |
|    **txStart** | Transcription start position (or end position for minus strand item) |
|      **txEnd** | Transcription end position (or start position for minus strand item) |
|   **cdsStart** | Coding region start (or end position if for minus strand item)       |
|     **cdsEnd** | Coding region end (or start position if for minus strand item)       |
|  **exonCount** | Number of exons                                                      |
| **exonStarts** | Exon start positions (or end positions for minus strand item)        |
|   **exonEnds** | Exon end positions (or start positions for minus strand item)        |
|  **proteinID** | UniProt display ID, UniProt accession, or RefSeq protein ID          |

Press **get output**.



### Mouse (mm9)

Go to [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables) and select the following:

|     Parameter     |                       Value                       |
| ----------------: | :-----------------------------------------------: |
|         **clade** | *Mammal*                                          |
|        **genome** | *Mouse*                                           |
|      **assembly** | *July 2007 (NCBI37/mm9)*                          |
|         **group** | *Genes and Gene Predictions*                      |
|         **track** | *UCSC Genes*                                      |
|         **table** | *knownCanonical*                                  |
| **output format** | *selected fields from primary and related tables* |
|   **output file** | `mm9.knownCanonical.full_table.tsv`               |

Press **get output**.


**Select Fields from mm9.knownCanonical**:

|     Field     |                          Description                          |
| ------------: | ------------------------------------------------------------- |
| **clusterId** | Which cluster of transcripts this belongs to in knownIsoforms |


Select the following **Linked Tables**:

|   go    |   goaPart   |                   Description                    |
| ------: | :---------: | ------------------------------------------------ |
| **mm9** | *kgXref*    | Link together a Known Gene ID and a gene alias   |
| **mm9** | *knownGene* | Transcript from default gene set in UCSC browser |

Press **allow selection from checked tables**.


Select the following fields in **mm9.kgXref fields**:

|     Field      | Description |
| -------------: | :---------: |
| **geneSymbol** | Gene Symbol |
|     **refseq** | RefSeq ID   |


Select the following fields in **mm9.knownGene fields**:

|     Field      |                             Description                              |
| -------------: | -------------------------------------------------------------------- |
|       **name** | Name of gene                                                         |
|      **chrom** | Reference sequence chromosome or scaffold                            |
|     **strand** | + or - for strand                                                    |
|    **txStart** | Transcription start position (or end position for minus strand item) |
|      **txEnd** | Transcription end position (or start position for minus strand item) |
|   **cdsStart** | Coding region start (or end position if for minus strand item)       |
|     **cdsEnd** | Coding region end (or start position if for minus strand item)       |
|  **exonCount** | Number of exons                                                      |
| **exonStarts** | Exon start positions (or end positions for minus strand item)        |
|   **exonEnds** | Exon end positions (or start positions for minus strand item)        |
|  **proteinID** | UniProt display ID, UniProt accession, or RefSeq protein ID          |

Press **get output**.



### Mouse (mm10)

Go to [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables) and select the following:

|     Parameter     |                       Value                       |
| ----------------: | :-----------------------------------------------: |
|         **clade** | *Mammal*                                          |
|        **genome** | *Mouse*                                           |
|      **assembly** | *Dec. 2011 (GRCm38/mm10)*                         |
|         **group** | *Genes and Gene Predictions*                      |
|         **track** | *GENCODE VM20*                                    |
|         **table** | *knownCanonical*                                  |
| **output format** | *selected fields from primary and related tables* |
|   **output file** | `mm10.knownCanonical.full_table.tsv`              |

Press **get output**.


**Select Fields from mm10.knownCanonical**:

|     Field     |                          Description                          |
| ------------: | ------------------------------------------------------------- |
| **clusterId** | Which cluster of transcripts this belongs to in knownIsoforms |


Select the following **Linked Tables**:

|    go    |   goaPart   |                   Description                    |
| -------: | :---------: | ------------------------------------------------ |
| **mm10** | *kgXref*    | Link together a Known Gene ID and a gene alias   |
| **mm10** | *knownGene* | Transcript from default gene set in UCSC browser |

Press **allow selection from checked tables**.


Select the following fields in **mm10.kgXref fields**:

|     Field      | Description |
| -------------: | :---------: |
| **geneSymbol** | Gene Symbol |
|     **refseq** | RefSeq ID   |


Select the following fields in **mm10.knownGene fields**:

|     Field      |                             Description                              |
| -------------: | -------------------------------------------------------------------- |
|       **name** | Name of gene                                                         |
|      **chrom** | Reference sequence chromosome or scaffold                            |
|     **strand** | + or - for strand                                                    |
|    **txStart** | Transcription start position (or end position for minus strand item) |
|      **txEnd** | Transcription end position (or start position for minus strand item) |
|   **cdsStart** | Coding region start (or end position if for minus strand item)       |
|     **cdsEnd** | Coding region end (or start position if for minus strand item)       |
|  **exonCount** | Number of exons                                                      |
| **exonStarts** | Exon start positions (or end positions for minus strand item)        |
|   **exonEnds** | Exon end positions (or start positions for minus strand item)        |
|  **proteinID** | UniProt display ID, UniProt accession, or RefSeq protein ID          |

Press **get output**.



## Convert knownCanonical transcription start sites from UCSC to TSS file for `create_regulatory_domains.py`

Convert knownCanonical transcription start sites fetched from
[UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables) to a TSV file that can be used by
`create_regulatory_domains.py`, which reimplements *Basal plus extension* section of
*Associating genomic regions with genes* of [GREAT](http://bejerano.stanford.edu/great/public/html/) and which is used
internally by `mucistarget.py`.


```bash
convert_knownCanonical_full_table_to_GREAT_TSS_file () {
    local knownCanonical_full_table_filename="${1}"
    local assembly="${2}"
    
    awk -F '\t' -v 'OFS=\t' -v assembly="${assembly}" '
        {
            if (NR == 1) {
                # Store header names to column indexes.

                # Strip of leading "#" for first column of the header.
                header_name_to_column_idx[substr($1, 2)] = 1;

                for (idx = 2; idx <= NF; idx ++) {
                    header_name_to_column_idx[$(idx)] = idx;
                }

                # Get column indexes for columns we will use.
                gene_symbol_column_idx = header_name_to_column_idx[assembly ".kgXref.geneSymbol"];
                ensembl_transcript_id_column_idx = header_name_to_column_idx[assembly ".knownGene.name"];
                chrom_column_idx = header_name_to_column_idx[assembly ".knownGene.chrom"];
                strand_column_idx = header_name_to_column_idx[assembly ".knownGene.strand"];
                tx_start_column_idx = header_name_to_column_idx[assembly ".knownGene.txStart"];
                tx_end_column_idx = header_name_to_column_idx[assembly ".knownGene.txEnd"];
                protein_id_column_idx = header_name_to_column_idx[assembly ".knownGene.proteinID"];

                # If one of the column names was not found the index will be zero.
                if (gene_symbol_column_idx * ensembl_transcript_id_column_idx * chrom_column_idx * strand_column_idx * tx_start_column_idx * tx_end_column_idx * protein_id_column_idx == 0) {
                    print "Error: One or more header columns could not be found." > "/dev/stderr";
                    exit(1);
                }
            } else {
                gene_symbol = $(gene_symbol_column_idx);
                ensembl_transcript_id = $(ensembl_transcript_id_column_idx);
                chrom = $(chrom_column_idx);
                strand = $(strand_column_idx);
                tx_start = $( tx_start_column_idx);
                tx_end = $(tx_end_column_idx);
                protein_id = $(protein_id_column_idx);

                # Delete arrays so they can be repopulated from scratch.
                delete chroms;
                delete tx_starts;
                delete tx_ends;

                if (chrom == "chrX,chrY,") {
                    # Remove trailing comma for all field for genes that
                    # are located at the same position in chrX and chrY.
                    sub(/,$/, "", gene_symbol);
                    sub(/,$/, "", ensembl_transcript_id);
                    sub(/,$/, "", strand);
                    sub(/,$/, "", tx_start);
                    sub(/,$/, "", tx_end);
                    sub(/,$/, "",  protein_id);

                    chroms[1] = "chrX";
                    chroms[2] = "chrY";

                    nbr_tx_starts = split(tx_start, tx_starts, ",");
                    nbr_tx_ends = split(tx_end, tx_ends, ",");

                    # Fill in the transcription start or end coordinate for chrY as chrX
                    # if no specific coordinate is provided for chrY.
                    if (nbr_tx_starts == 1) {
                        tx_starts[2] = tx_starts[1];
                    }

                    if (nbr_tx_ends == 1) {
                        tx_ends[2] = tx_ends[1];
                    }
                } else {
                    chroms[1] = chrom;
                    tx_starts[1] = tx_start;
                    tx_ends[1] = tx_end;
                }

                # Only keep genes that will be transcribed to a protein.
                if (protein_id != "") {
                    # Loop over all chromosomes:
                    #    - One chromosome in most cases.
                    #    - "chrX" and chrY" for genes located at same location in chrX and chrY.

                    for (chrom_idx in chroms) {
                        if (strand == "-") {
                            # Take transcription end minus 1 when gene is on negative strand (0-based TSS coordinate).
                            print ensembl_transcript_id, chroms[chrom_idx], tx_ends[chrom_idx] - 1, strand, gene_symbol;
                        } else {
                            # Take transcription start when gene is on positive strand (0-based TSS coordinate).
                            print ensembl_transcript_id, chroms[chrom_idx], tx_starts[chrom_idx], strand, gene_symbol;
                        }
                    }
                }
            }
        }' \
        "${knownCanonical_full_table_filename}" \
      | sort -k 2V -k 3n
}
```


Convert `${assembly}.knownCanonical.full_table.tsv` to `${assembly}.tss.tsv` for assemblies of interest.

```bash
# Specify assemblies of interest.
assemblies='hg19 hg38 mm9 mm10';

for assembly in ${assemblies} ; do
    # Convert knownCanonical full table to TSS TSV file.
    convert_knownCanonical_full_table_to_GREAT_TSS_file "${assembly}.knownCanonical.full_table.tsv" "${assembly}" > "${assembly}.tss.tsv";
    
    # Create BED TSS file for loading in UCSC.
    awk -F '\t' -v 'OFS=\t' '{print $2, $3, $3 + 1, $5, "1000", $4}' "${assembly}.tss.tsv" > "${assembly}.tss.bed";
done
```

If you want to check the regulatory domains that will be used by `mucistarget.py`, run the following from the root of
this git repo:

```bash
# Specify assemblies of interest.
assemblies='hg19 hg38 mm9 mm10';

# Go to root of mucistarget git repo or set the following variable correctly.
mucistarget_git_root_dir='.'

for assembly in ${assemblies} ; do
    # Create regulatory domains BED file for assembly: "${mucistarget_git_root_dir}/data/regulatory_domains/${assembly}.regdoms.bed"
    ./create_regulatory_domains.py \
        --basal-up 5000 \
        --basal-down 1000 \
        --max-ext 1000000 \
        --chrom-sizes "${mucistarget_git_root_dir}/data/genomic_fasta/${assembly}.chrom.sizes" \
        --assembly "${assembly}" \
        --genes-tss "${mucistarget_git_root_dir}/data/regulatory_domains/${assembly}.tss.tsv" \
        --regdoms "${mucistarget_git_root_dir}/data/regulatory_domains/${assembly}.regdoms.bed";
done
```