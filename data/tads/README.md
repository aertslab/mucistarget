# TADs: topologically associating domains

## Download TADs

Download TADs in hg19, hg38, mm9 and mm10 predicted using pipeline described in
[Topological domains in mammalian genomes identified by analysis of chromatin interactions.
 Dixon JR, Selvaraj S, Yue F, Kim A, Li Y, Shen Y, Hu M, Liu JS, Ren B.
 Nature. 2012 Apr 11;485(7398):376-80. doi: 10.1038/nature11082.
](https://www.ncbi.nlm.nih.gov/pubmed/22495300)

See [http://promoter.bx.psu.edu/hi-c/publications.html](http://promoter.bx.psu.edu/hi-c/publications.html)
for included datasets.

```bash
# Download TAD info.
wget 'http://promoter.bx.psu.edu/hi-c/downloads/TAD-info.txt'

# Download TADs.
wget 'http://promoter.bx.psu.edu/hi-c/downloads/hg19.TADs.zip'
wget 'http://promoter.bx.psu.edu/hi-c/downloads/hg38.TADs.zip'
wget 'http://promoter.bx.psu.edu/hi-c/downloads/mm9.TADs.zip'
wget 'http://promoter.bx.psu.edu/hi-c/downloads/mm10.TADs.zip'
```



## Create TAD BED files from raw files.

Create TAD BED files from raw files:
  - Cleanup filename a little bit.
  - Add "chr" to chromosome name if necessary.


```bash
# Extract TAD files in assembly specific directory.
unzip hg19.TADs.zip -d hg19
unzip mm9.TADs.zip -d mm9

unzip hg38.TADs.zip
unzip mm10.TADs.zip

for assembly in hg19 hg38 mm9 mm10 ; do
    # Loop over all raw TAD files.
    for tad_file in ${assembly}/*.txt ${assembly}/*.domains ; do
        # Skip ${tad_file} which contain "*.txt" or "*.domains" literally (no files matching glob pattern).
        if [ "${tad_file/\*./}" = "${tad_file}" ] ; then
            # Remove "-raw*" "_raw*" from TAD filename.
            tad_bed_file="${tad_file%?raw*}";

            # Remove directory part of TAD BED filename.
            TAD_name="${tad_bed_file##*/}";

            # Add "TADs_" in the beginning of the TAD BED file basename.
            tad_bed_file="${tad_bed_file/\///TADs_}";

            # Add ".bed" extension to TAD BED filename.
            tad_bed_file="${tad_bed_file}.bed";


            echo "Creating \"${tad_bed_file}\" with name \"${TAD_name}\".";

            # Write a 4 column BED file:
            #   - Add "chr" to chromosome name in necessary.
            #   - Add TAD name as column 4.
            awk -F '\t' -v 'OFS=\t' -v TAD_name="${TAD_name}" \
                '{ if ($1 !~ /^chr/) { $1 = "chr" $1; } print $1, $2, $3, TAD_name; }' \
                "${tad_file}" \
              > "${tad_bed_file}";
        fi
    done
done
```



## Remove raw TAD files.

Raw TAD files can be removed.

```bash
# Remove raw TAD files as they are converted to BED files.
for assembly in hg19 hg38 mm9 mm10 ; do
    # Loop over all raw TAD files.
    rm ${assembly}/*.txt ${assembly}/*.domains
done

# Remove downloaded TAD archive files.
rm hg19.TADs.zip hg38.TADs.zip mm9.TADs.zip mm10.TADs.zip

```
