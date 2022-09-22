# bs-ATLAS-seq analysis
## Background
Script to call L1 insertions and their methylation levels in bs-ATLAS-seq experiments.

## Installation
### Dependencies
- [cutadapt](https://github.com/marcelm/cutadapt) (tested with v3.1)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.4.1)
- [bismark](https://github.com/FelixKrueger/Bismark) (version >= 0.22.0 which supports soft-clipping, tested with v0.22.1)
- [samtools](http://www.htslib.org) (tested with v1.3)
- [Picard tools](http://broadinstitute.github.io/picard/) (tested with v2.19)
- [bedtools](https://github.com/arq5x/bedtools2) (tested with v2.25.0)
- [methpat](https://bjpop.github.io/methpat/) (tested with v2.1.0)
- [seqtk](https://github.com/lh3/seqtk) (tested with v1.3)
- [GNU parallel](https://www.gnu.org/software/parallel/) (tested v20200922)

**Make sure these programs are included in your $PATH directory.**

### Procedure
Note that `<PATH>` should be replaced by the name of the folder in which you have installed the `bs-ATLAS-seq/` folder.

1. Download the whole project folder from GitHub:
```
git clone https://github.com/retrogenomics/bs-ATLAS-seq.git
```
2. Download a human reference genome (e.g. hg38):
````
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O <PATH>/bs-ATLAS-seq/references/hg38/hg38.fa.gz
gunzip <PATH>/bs-ATLAS-seq/references/hg38/hg38.fa.gz
```
3. Prepare Bismark genomes and indexes:
```
bismark_genome_preparation --verbose <PATH>/bs-ATLAS-seq/references/hg38/
bismark_genome_preparation --verbose <PATH>/bs-ATLAS-seq/references/L1/
```
4. Edit the `config_dir.sh` file in the `bs-ATLAS-seq/` folder according to your configuration

## How to use?
### To process sequencing data and get L1 insertions
```bash
<PATH>/bs-ATLAS-seq/scripts/bs-atlas-seq_calling.sh \
   -1 Sample1_R1.fastq.gz \
   -2 Sample2_R2.fastq.gz \
   -o <OUTPUT_DIRECTORY> \
   -p "${k}" -n 10 -s 2 -t 12
```
