# Genomic FASTA files


## Download genomic 2bit and chromosome sizes files.

Download genomic 2bit and chromosome sizes files from [UCSC goldenPath](http://hgdownload.cse.ucsc.edu/goldenPath/):

```bash
# Specify assemblies of interest.
assemblies='hg19 hg38 mm9 mm10';

for assembly in ${assemblies} ; do
    printf "\nDownloading 2bit and chromosome file for %s ...\n\n" "${assembly}";

    # Download 2bit file for assembly.
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.2bit
    rsync -aP rsync://hgdownload.soe.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.2bit .

    # Download chromosome sizes file for assembly (only needed when manually running "create_regulatory_domains.py").
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.chrom.sizes
    rsync -aP rsync://hgdownload.soe.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.chrom.sizes .
done
```



## Convert genomic 2bit file to FASTA file.

Download `twoBitToFa` or add `twoBitToFa` to your `PATH` if you already have it:

```bash
# Download precompiled version of twhBitToFa.
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/twoBitToFa .

# Make twoBitToFa binary executable.
chmod a+x twoBitToFa

export PATH="${PWD}:${PATH}"
```

Convert genomic 2bit file to FASTA file with `twoBitToFa`.

```bash
# Specify assemblies of interest.
assemblies='hg19 hg38 mm9 mm10';

for assembly in ${assemblies} ; do
    # Convert 2bit file to FASTA file for each assembly.
    twoBitToFa "${assembly}.2bit" "${assembly}.fa";
done

# 2bit files are not needed by mucistarget and can be removed.
for assembly in ${assemblies} ; do
    # Convert 2bit file to FASTA file for each assembly.
    rm "${assembly}.2bit";
done
```



## Create flattened FASTA file.

Load each FASTA file with `pyfasta` so it creates:
  - flattened FASTA file: `${assemlby}.fa.flat`
  - pickle of the start, stop locations of each header in the FASTA file: `${assemlby}.fa.gdx`


Start python3 interpreter and run:
```python
import pyfasta

pyfasta.Fasta('hg19.fa')
pyfasta.Fasta('hg38.fa')
pyfasta.Fasta('mm9.fa')
pyfasta.Fasta('mm10.fa')
```
