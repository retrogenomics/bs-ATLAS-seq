# bs-ATLAS-seq analysis
## Background
Script to call L1 insertions and their methylation levels in bs-ATLAS-seq experiments.

## Installation
### Dependencies
Make sure to include the path to these programs in your `$PATH` variable by editing the appropriate file depending on your system (it could be ~/.bash_profile, ~/.bashrc, or ~/.profile)
- [cutadapt](https://github.com/marcelm/cutadapt) (tested with v3.1)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested with v2.4.1)
- [bismark](https://github.com/FelixKrueger/Bismark) (version >= 0.22.0 which supports soft-clipping, tested with v0.22.1)
- [samtools](http://www.htslib.org) (tested with v1.3)
- [Picard tools](http://broadinstitute.github.io/picard/) (tested with v2.19)
- [bedtools](https://github.com/arq5x/bedtools2) (tested with v2.25.0)
- [methpat](https://bjpop.github.io/methpat/) (tested with v2.1.0)
- [seqtk](https://github.com/lh3/seqtk) (tested with v1.3)
- [GNU parallel](https://www.gnu.org/software/parallel/) (tested v20200922)

### Procedure
Note that `<DOWNLOAD>` is the name of the folder in which you have downloaded the repository.

1. Download the project repository from GitHub:
```
wget 'https://github.com/retrogenomics/bs-ATLAS-seq/archive/refs/heads/main.zip' -O <DOWNLOAD>/main.zip
cd <DOWNLOAD>
unzip main.zip
mv bs-ATLAS-seq-main/ bs/
```
2. Download a human reference genome (e.g. hg38):
```
mkdir -p bs/references/hg38
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O bs/references/hg38/hg38.fa.gz
gunzip bs/references/hg38/hg38.fa.gz
```
3. Prepare Bismark genomes and indexes:
Note that human genome preparation can take several hours.
```
bismark_genome_preparation --verbose bs/references/L1/
bismark_genome_preparation --verbose bs/references/hg38/
```
4. Edit the `config_dir.sh` file in the `bs-ATLAS-seq/` folder according to your configuration

## How to use?
### To process sequencing data and get L1 insertions
```bash
<DOWNLOAD>/bs/scripts/bs-atlas-seq_calling.sh \
   -1 Sample1_R1.fastq.gz \
   -2 Sample2_R2.fastq.gz \
   -o ./results \
   -p Sample1 -n 10 -s 2 -t 8
```
