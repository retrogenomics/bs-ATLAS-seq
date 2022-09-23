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
4. Prepare repeatmasker file:
```
gunzip bs/annotations/hg38_rmsk_L1_repPercent.bed.gz
```
5. Prepare Bismark human reference genome:
Note that this can take several hours.
```
bismark_genome_preparation --verbose bs/references/hg38/
```
6. Edit the `config_dir.sh` file in the `bs-ATLAS-seq/` folder according to your configuration

## How to use?
```
usage:	bs-atlas-seq_calling.sh [options] -1 read1.fastq.gz -2 read2.fastq.gz -p prefix
options:
	-h Print this help menu.
	-v What version of bs-atlas-seq_calling are you using?
	-o Output directory.
	-n Total read number threshold to call an insertion [default=3]
	-s Split read number threshold to call an insertion [default=2]
	-u Subsampling of input fastq file (no=1; or indicate fraction of reads to consider, e.g. 0.01) [default=1]
	-t Number of threads used for mapping [default=4]
	-c Activate cleanup (deletion of temporary files) [default=off]
	-d Configuration file [default=./config_dir.sh]
```
Example:
The output directory will be `~/Downloads/test_220923` and output files will all be prefixed with `test_IMR90.`. Here we use the subsampling `-u` option to use only 1% of the reads for testing purposes.
```
./bs/scripts/bs-atlas-seq_calling.sh \
	-1 ~/Lab/bioinfo/data/bs-atlas-seq/180525_NB501350_0027_AH37TJBGX7/IMR90_S6_R1.fastq.gz \
	-2 ~/Lab/bioinfo/data/bs-atlas-seq/180525_NB501350_0027_AH37TJBGX7/IMR90_S6_R2.fastq.gz \
	-o ~/Downloads/test_220923 \
	-p test_IMR90 \
	-n 3 -s 1 -t 6 -u 0.01 \
	-d ~/Downloads/bs/config_dir.sh
```
## Output
The pipeline generates a number of files which are described below:
| File name | Description |
|---|---|
|`<PREFIX>.bs-atlas-seq_calling.sh.script.sh`|A copy of the script used to generate the data.|
|`<PREFIX>.stats.txt`|A summary file with the number of reads, processed at each step, and the number of L1 detected.|
|`<PREFIX>.log`|A log file with the information displayed on the screen while the script is running.|
|`<PREFIX>.hg38.bam`|Alignment of the reads to the reference genome provided.|
|`<PREFIX>.hg38.bai`|Index file of the `.bam` file.|
|`<PREFIX>.all_L1.flanks_and_internal.hg38.bedGraph.gz`|BedGraph file with mCG levels.|
|`<PREFIX>.all_L1.hg38.bed`|Bed file containing all called L1 as well as their average mCG level as score (5th column).|
|`<PREFIX>.L1HS.hg38.bed`|Similar to `<PREFIX>.all_L1.hg38.bed` but filtered to keep only L1HS elements.|
|`<PREFIX>.all_L1.hg38.html`|An `.html` file displaying the single-molecule profiles of methylation for each L1 locus identified. Can be opened with a web browser (tested with Chrome).|
|`<PREFIX>.L1HS.hg38.html`|Similar to `<PREFIX>.all_L1.hg38.html` but filtered to keep only L1HS elements.|

 




