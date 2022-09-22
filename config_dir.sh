#!/usr/bin/env bash

# path to root directory
path_to_bioinfo="${HOME}/bs-ATLAS-seq"

# path to software dependencies
path_to_bwt2="/opt/local/bin"
picard="/usr/local/picard/picard.jar"

# reference sequences and bismark indexes
reference="hg38"
path_to_reference="${path_to_bioinfo}/datasets/references/hg38"
path_to_L1HS="${path_to_bioinfo}/datasets/references/L1/bismark-bwt2"

# path to annotation files
rmsk="${path_to_bioinfo}/datasets/annotations/hg38_rmsk_L1_repPercent.bed"
genes="${path_to_bioinfo}/datasets/annotations/hg38_NCBIRefSeq_refGene.bed"
blacklist="${path_to_bioinfo}/datasets/annotations/hg38.encode4_unified_blacklist.ENCFF356LFX.bed"

echo -e "Default directories loaded."
