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

### Other requirements
- A human reference genome sequence (ex: [`hg38.fa`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/))
- and its bismark index (ex: `hg19.fa.amb, .ann, .bwt, .fai, .pac, .sa`)

### Procedure

1. Download from GitHub:
```
git clone https://github.com/retrogenomics/bs-ATLAS-seq.git
```
2. Edit the `CONFIG` file in the `iss/scripts/` folder according to your configuration
3. Edit the experiment configuration file (sequencing metadata) according to your sequencing data.
   - This file is stored in the `iss/experiments/` folder and its name.
   - Its location should be included in the `CONFIG`file as `EXP_CONFIG="${EXPERIMENTS}/<your file name>"`.
   - `atlas-neo-R01-to-R09.txt` is the default example.

## How to use?
### To process sequencing data and get L1 insertions
```bash
cd iss/scripts
./atlas-neo_all_samples.sh
```
