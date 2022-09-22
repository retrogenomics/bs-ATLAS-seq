#!/usr/bin/env bash

################################################################################
# L1 and methylation calling for bs-atlas-seq experiments                      #
################################################################################

# dependencies -----------------------------------------------------------------
#   - cutadapt 3.1 (python 3.8)
#   - bowtie2 2.4.1
#   - bismark 0.22.1
#   - samtools 1.3
#   - picard tools
#   - bedtools 2.29.2
#   - methpat 2.1.0 (python 2.7)
#   - seqtk 1.3-r115-dirty
#   - GNU parallel 20200922

###### load preferences and define global settings #############################
script_name="bs-atlas-seq_calling"
script_version='1.2.0'

start_time=$( date +%s )
day=$(date +"[%d-%m-%Y] [%T]")
starline=$( printf "*%.0s" {1..105} ) # separator for output
step=1	# store progress through the pipeline

LC_NUMERIC_OLD=$LC_NUMERIC
export LC_NUMERIC="en_US.UTF-8"

shopt -s extglob

###### script argument parsing and processing ##################################

# set defaults argument values -------------------------------------------------
distance=100 # window (bp) for merging reads into clusters
threads=4
subsampling=1
read_threshold=3
splitread_threshold=2
output_dir="$( pwd )"
cleanup="off"
config_file="./config_dir.sh"

# store usage explanations -----------------------------------------------------
USAGE="\
$script_name v$script_version:\tanalysis of bs-ATLAS-seq runs from demultiplexed fastq.gz files to L1 and CG methylation calling. \n\n\
usage:\t$( basename $0 ) [options] -1 read1.fastq.gz -2 read2.fastq.gz -p prefix \n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
\t-o Output directory. \n\
\t-n Total read number threshold to call an insertion [default=$read_threshold] \n\
\t-s Split read number threshold to call an insertion [default=$splitread_threshold] \n\
\t-u Subsampling of input fastq file (no=1; or indicate fraction of reads to consider, e.g. 0.01) [default=$subsampling] \n\
\t-t Number of threads used for mapping [default=$threads]\n\
\t-c Activate cleanup (deletion of temporary files) [default=off]\n\
\t-d Configuration file [default=${config_dir}\n\
"

# parse script arguments -------------------------------------------------------
while getopts ':hcv1:2:p:o:t:s:n:u:d:' opt ; do
	case $opt in
    1) input_file_R1=$OPTARG ;;
    2) input_file_R2=$OPTARG ;;
    o) output_dir=$OPTARG ;;
		n) read_threshold=$OPTARG ;;
		s) splitread_threshold=$OPTARG ;;
    p) prefix=$OPTARG ;;
		t) threads=$OPTARG ;;
		u) subsampling=$OPTARG ;;
		d) config_dir=$OPTARG ;;
		c) cleanup="on" ;;
		h) echo -e "\n$USAGE"; exit 1 ;;
		v) echo -e "${script_name} v${script_version}" ; exit 1 ;;
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e "${USAGE}"; exit 1 ;;
	esac
done

# check for mandatory positional parameters
if [[ -z "${input_file_R1}" || ! -f "${input_file_R1}" ]];
then
	echo -e "\nInput file read1.fastq.gz not specified or not existing.\n";
	echo -e "${USAGE}" ; exit 1
fi

if [[ -z "${input_file_R2}" || ! -f "${input_file_R2}" ]];
then
	echo -e "\nInput file read2.fastq.gz not specified or not existing.\n";
	echo -e "${USAGE}" ; exit 1
fi

if [[ -z "${prefix}" ]];
then
	echo -e "\nPrefix not specified.\n";
	echo -e "${USAGE}" ; exit 1
fi

# set an upper limit for the number of parallel instances of bismark to be run concurrently
# this is necessary, since each bismark instance already use 4 bowtie2 process.
if [[ "${threads}" -ge 3 ]];
then
	bs_parallel=3;
else
	bs_parallel="${threads}"
fi

# define environmental variables and paths -------------------------------------
if [[ -z "${config_file}" || ! -f "${config_file}" || ! -s "${config_file}" ]];
then
	echo -e "\nConfiguration file not specified, not existing or empty.\n";
	echo -e "${USAGE}" ; exit 1
else
	source "${config_file}"
	echo -e "\nConfiguration loaded from: ${config_file}.\n";
fi
wd="${output_dir}/temp"
mkdir -p "${output_dir}" "${wd}"

# print header
output="\n\
${starline}\n\
$day bs-ATLAS-seq analysis - Script v$script_version \n\
$starline\n\
Command line:\n$( basename $0 ) $@ \n\
${starline}\n\
PE-sequencing data: \n\
\t- R1: ${input_file_R1} \n\
\t- R2: ${input_file_R2} \n\
Ref genome:\n\
\t- ${reference}: ${path_to_reference}/${reference}.fa \n\
Sample:\n\
\t- ${prefix}\n\
${starline}\n\
Steps:\n\
"
log+="${output}"
echo -ne "${output}" | fold -w 105


####### test if subsampling required and if yes, generate subsampled input file

output=$( printf "%-3s- Data loading " "${step}" )
log+="${output}"
echo -ne "${output}"

if [ $( echo "${subsampling}<1" | bc -l ) -eq 1 ];
then
	printf "with subsampling to %s fraction of sequences..." "${subsampling}"
	R1="${wd}/sample_${subsampling}_$( basename "${input_file_R1}" )"
	R2="${wd}/sample_${subsampling}_$( basename "${input_file_R2}" )"
	seqtk sample -s 100 "${input_file_R1}" "${subsampling}" \
	| gzip \
	> "${R1}"
	seqtk sample -s 100 "${input_file_R2}" "${subsampling}" \
	| gzip \
	> "${R2}"
else
	printf "without subsampling..."
	R1="${input_file_R1}"
	R2="${input_file_R2}"
fi

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))


##### pre-processing ###########################################################

# remove adapters and low quality reads ----------------------------------------
# discard pairs not containing L1 primer in R2, and too short
output=$( printf "%-3s- Read trimming..." "${step}" )
log+="${output}"
echo -ne "${output}"

cutadapt \
  --cores="${threads}" \
  --nextseq-trim=20 \
  --minimum-length=36:189 \
  --trim-n \
  --discard-untrimmed \
  --action=none \
  --quiet \
  -G L1_target_rc="^AAAAACTCCCTAACCCCTTA;e=0.15;o=5" \
  -o "${wd}/${prefix}_R1.trimmed.fastq.gz" \
  -p "${wd}/${prefix}_R2.trimmed.fastq.gz" \
  "${R1}" \
  "${R2}"

# remove trailing adapter sequences from mates ---------------------------------
# observed when sequencing of R1 reached end of R2 and vice versa (should be rare)
cutadapt \
  --cores="${threads}" \
  --nextseq-trim=20 \
  --minimum-length=36:189 \
  --overlap=5 \
  --error-rate=0.15 \
  --trim-n \
  --quiet \
  -a L1_target="GTAAGGGGTTAGGGAGTTTTT" \
  -A Rd1SP_rc="AAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
  -o "${wd}/${prefix}_R1.retrimmed.fastq.gz" \
  -p "${wd}/${prefix}_R2.retrimmed.fastq.gz" \
  "${wd}/${prefix}_R1.trimmed.fastq.gz" \
  "${wd}/${prefix}_R2.trimmed.fastq.gz"

# remove intermediate files ----------------------------------------------------
rm "${wd}/${prefix}_R"*.trimmed.fastq.gz
mv "${wd}/${prefix}_R1.retrimmed.fastq.gz" "${wd}/${prefix}_R1.fastq.gz"
mv "${wd}/${prefix}_R2.retrimmed.fastq.gz" "${wd}/${prefix}_R2.fastq.gz"

# read count for report --------------------------------------------------------
count_raw=$( pigz -cdk "${R1}" | wc -l | awk '{print $1/4}' )
count_trimmed=$( pigz -cdk "${wd}/${prefix}_R1.fastq.gz" | wc -l | awk '{print $1/4}' )

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

###### identify PE reads with R2 mapping to L1 #################################

output=$( printf "%-3s- Mapping R2 to L1..." "${step}" )
log+="${output}"
echo -ne "${output}"

# build bismark index for L1HS 5end --------------------------------------------
bismark_genome_preparation \
  --bowtie2 \
  --path_to_aligner "${path_to_bwt2}" \
  "${path_to_L1HS}" #&> /dev/null

# map R2 mates to L1HS ---------------------------------------------------------
bismark \
  --genome "${path_to_L1HS}" \
  "${wd}/${prefix}_R2.fastq.gz" \
  --output_dir "${wd}" \
  --unmapped \
  --bowtie2 \
  --path_to_bowtie2 "${path_to_bwt2}" \
  --local \
  --pbat \
  --multicore "${bs_parallel}" \
	#&> /dev/null

mv "${wd}/${prefix}_R2_bismark_bt2.bam" "${wd}/${prefix}.R2toL1HS.bam"
mv "${wd}/${prefix}_R2.fastq.gz_unmapped_reads.fq.gz" "${wd}/${prefix}.R2toL1HS.unmapped.fastq.gz"

# list fragments with R2 mate mapping to L1HS and correctly oriented -----------
samtools view -F 16 "${wd}/${prefix}.R2toL1HS.bam" \
| awk '{split($1,read,"_"); print read[1]}' \
> "${wd}/${prefix}.R2toL1HS.list"

# extract raw sequence for pairs with R2 mate mapping to L1HS ------------------
for mate in R1 R2;
do
  seqtk subseq "${wd}/${prefix}_${mate}.fastq.gz" "${wd}/${prefix}.R2toL1HS.list" \
  | awk -F " " '{if (NR%4==1) {print $1} else {print $0}}' \
  | gzip \
  > "${wd}/${prefix}.L1HS-filtered.${mate}.fastq.gz"
done

# read count for report --------------------------------------------------------
count_mapped_1=$( pigz -cdk "${wd}/${prefix}.L1HS-filtered.R1.fastq.gz" | wc -l | awk '{print $1/4}' )

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

###### identify properly paired and discordant reads mapped to ref genome ######

output=$( printf "%-3s- Mapping PE reads to %s..." "${step}" "${reference}" )
log+="${output}"
echo -ne "${output}"

# map to ref genome PE fragments with R2 mate matching L1HS --------------------
# end-to-end mode by default, with score_min adjusted as in bs-seeker2 ---------

bismark \
  --genome "${path_to_reference}" \
  -1 "${wd}/${prefix}.L1HS-filtered.R1.fastq.gz" \
  -2 "${wd}/${prefix}.L1HS-filtered.R2.fastq.gz" \
  --output_dir "${wd}" \
  --unmapped \
  --bowtie2 \
  --path_to_bowtie2 "${path_to_bwt2}" \
  --minins 250 \
  --maxins 1250 \
	--score_min L,-0.6,-0.6 \
	--multicore "${bs_parallel}"	\
  &> /dev/null

mv "${wd}/${prefix}.L1HS-filtered.R1.fastq.gz_unmapped_reads_1.fq.gz" "${wd}/${prefix}.L1HS-filtered.unmapped_R1.fastq.gz"
mv "${wd}/${prefix}.L1HS-filtered.R2.fastq.gz_unmapped_reads_2.fq.gz" "${wd}/${prefix}.L1HS-filtered.unmapped_R2.fastq.gz"

# recalibrate MAPQ scores ------------------------------------------------------
# set to 50, allowing proper visualization in IGV
samtools view -h "${wd}/${prefix}.L1HS-filtered.R1_bismark_bt2_pe.bam" \
| awk -v OFS="\t" '
	($1~/^@/) {print}
	($1!~/^@/) {
		printf $1 "\t" $2 "\t" $3 "\t" $4 "\t50\t";
		for (i=6;i<NF;i++) {printf $i "\t"};
		printf $NF "\n";
	}
 ' \
| samtools view -b - \
> "${wd}/${prefix}.L1HS-filtered.PEto${reference}.bam"

# sort reads by coordinates ----------------------------------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} SortSam \
  --INPUT "${wd}/${prefix}.L1HS-filtered.PEto${reference}.bam" \
  --OUTPUT "${wd}/${prefix}.ref-L1.${reference}.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
	--SORT_ORDER coordinate \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --CREATE_INDEX TRUE

# list fragments not mapped as pairs to hg38 (discordant pairs) ----------------
zcat "${wd}/${prefix}.L1HS-filtered.unmapped_R1.fastq.gz" \
| sed -n '1~4p' \
| sed 's/^@//g' \
> "${wd}/${prefix}.discordant_PE_${reference}.list"

# extract R1 raw sequence for discordant pairs ---------------------------------
seqtk subseq "${wd}/${prefix}_R1.fastq.gz" "${wd}/${prefix}.discordant_PE_${reference}.list" \
| awk -F " " '{if (NR%4==1) {print $1} else {print $0}}' \
| gzip \
> "${wd}/${prefix}.discordant_PE_${reference}.R1.fastq.gz"

# calculate covered regions for ref-L1 insertions ------------------------------
samtools view -b -f 3 "${wd}/${prefix}.ref-L1.${reference}.bam" \
| samtools sort -@ "${threads}" -n \
| bedtools bamtobed -bedpe -mate1 -i - \
| awk -v OFS="\t" '
    ($NF=="-") {print $1,$2,$6,$7,$6-$2,$9}
    ($NF=="+") {print $1,$5,$3,$7,$3-$5,$9}
  ' \
| sort -k1,1 -k2,2n \
| bedtools merge -s -c 4,5,6 -o distinct,count,distinct -i - \
> "${wd}/${prefix}.refL1.PE-clusters.bed"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

###### remap discordant reads to ref genome (R1 only as SE read) ###############

output=$( printf "%-3s- Mapping discordant R1 reads to %s..." "${step}" "${reference}" )
log+="${output}"
echo -ne "${output}"

# map R1 from discordant pairs to ref genome -----------------------------------
# local mode
bismark \
  --genome "${path_to_reference}" \
  "${wd}/${prefix}.discordant_PE_${reference}.R1.fastq.gz" \
  --output_dir "${wd}" \
  --unmapped \
  --bowtie2 \
  --path_to_bowtie2 "${path_to_bwt2}" \
  --local \
  --multicore "${bs_parallel}" \
	&> /dev/null

mv "${wd}/${prefix}.discordant_PE_${reference}.R1.fastq.gz_unmapped_reads.fq.gz" \
	"${wd}/${prefix}.discordant_PE_${reference}.unmapped_R1.fastq.gz"

# filter out R1 that are: ------------------------------------------------------
# 	- soft-clipped at their 5' end or
# 	- mapping to clusters of PE reads
bedtools intersect -v \
	-abam "${wd}/${prefix}.discordant_PE_${reference}.R1_bismark_bt2.bam" \
	-b "${wd}/${prefix}.refL1.PE-clusters.bed" \
| samtools view -h  \
| awk -v OFS="\t" '($1~/^@/) || ($1!~/^@/ && $2==0 && $6!~/^[0-9]+S/) || ($1!~/^@/ && $2==16 && $6!~/[0-9]+S$/)' \
| samtools view -b -h \
> "${wd}/${prefix}.discordant_PE.R1to${reference}.bam"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

##### SPLIT-READS ##############################################################

output=$( printf "%-3s- Identify split R1 reads..." "${step}" )
log+="${output}"
echo -ne "${output}"

# minimal length of seed for split-read search ---------------------------------
seed=15

# extract 3' soft-clipped fragment of R1 reads to identify split reads ---------
samtools view "${wd}/${prefix}.discordant_PE.R1to${reference}.bam" \
| awk -v OFS="\n" -v seed="${seed}" '
	($2==0 && $6~/[0-9]+S$/) {
		n=split($6,cigar,"[S|M|I|D]");
		softclip_seq=substr($10,length($10)-cigar[n-1]+1)
		softclip_qual=substr($11,length($11)-cigar[n-1]+1)
		if (length(softclip_seq)>=seed) {print "@"$1,softclip_seq,"+",softclip_qual}
	}
	($2==16 && $6~/^[0-9]+S/){
		n=split($6,cigar,"[S|M|I|D]");
		softclip_seq=substr($10,1,cigar[1])
		softclip_qual=substr($11,1,cigar[1])
		if (length(softclip_seq)>=seed) {print "@"$1,softclip_seq,"+",softclip_qual}
	}' \
| gzip \
> "${wd}/${prefix}.split_reads_R1.fastq.gz"

# map soft-clipped region of R1 reads to L1HS ----------------------------------
bismark \
  --genome "${path_to_L1HS}" \
  "${wd}/${prefix}.split_reads_R1.fastq.gz" \
  --output_dir "${wd}" \
  --bowtie2 \
  --path_to_bowtie2 "${path_to_bwt2}" \
  --local \
	--non_directional \
	-L "${seed}" \
	--score_min G,12,0.8 \
  --multicore "${bs_parallel}" \
	&> /dev/null

mv "${wd}/${prefix}.split_reads_R1_bismark_bt2.bam" \
	"${wd}/${prefix}.split_reads.R1toL1HS.bam"

# list split reads mapping both on ref genome and L1HS -------------------------
samtools view -F 16 "${wd}/${prefix}.split_reads.R1toL1HS.bam" \
| awk '{print $1}' \
> "${wd}/${prefix}.split_reads.R1toL1HS.reads"

# extract split reads mapped on ref genome -------------------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar "${picard}" FilterSamReads \
	--INPUT "${wd}/${prefix}.discordant_PE.R1to${reference}.bam" \
	--OUTPUT "${wd}/${prefix}.split_reads.R1to${reference}.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
	--READ_LIST_FILE "${wd}/${prefix}.split_reads.R1toL1HS.reads" \
	--FILTER includeReadList \
	--WRITE_READS_FILES FALSE \
	--USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
	--SORT_ORDER coordinate \
	--CREATE_INDEX TRUE

# dedup split reads mapped on ref genome ---------------------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} MarkDuplicates \
  --INPUT "${wd}/${prefix}.split_reads.R1to${reference}.bam" \
  --OUTPUT "${wd}/${prefix}.split_reads.R1to${reference}.dedup.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
  --METRICS_FILE "${wd}/${prefix}.split_reads.R1to${reference}.dedup.log" \
  --REMOVE_DUPLICATES TRUE \
  --ASSUME_SORT_ORDER coordinate \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES \
	--CREATE_INDEX TRUE

# cluster split reads mapped on ref genome -------------------------------------
bedtools bamtobed -i "${wd}/${prefix}.split_reads.R1to${reference}.dedup.bam" \
| sort -k1,1 -k2,2n \
| bedtools merge -s -d "${distance}" -c 4,5,6 -o distinct,count,distinct -i - \
> "${wd}/${prefix}.split_reads.R1to${reference}.clusters.bed"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

##### DISCORDANT READS #########################################################

output=$( printf "%-3s- Merge all reads..." "${step}" )
log+="${output}"
echo -ne "${output}"

# list fragments mapped as singletons to ref genome (discordant pairs) ---------
samtools view "${wd}/${prefix}.discordant_PE.R1to${reference}.bam" \
| awk '{print $1}' \
> "${wd}/${prefix}.singletons.R1to${reference}.reads"

# extract R2 raw sequence for discordant pairs ---------------------------------
seqtk subseq "${wd}/${prefix}_R2.fastq.gz" "${wd}/${prefix}.singletons.R1to${reference}.reads" \
| awk -F " " '{if (NR%4==1) {print $1} else {print $0}}' \
| gzip \
> "${wd}/${prefix}.discordant_R2.fastq.gz"

# identify HQ R2 reads (most likely corresponding to non-ref L1HS - to exclude
# chimeric reads with older L1PA families -> biggest difference with version 1.0
# of this script)

cutadapt \
  --cores="${threads}" \
	--minimum-length=189 \
  --overlap 100 \
	--trim-n \
	--action=none \
  --discard-untrimmed \
	--error-rate=0.045 \
  --quiet \
	-a L1_target="AAAAACTCCCTAACCCCTTACRCTTCCCARRTRARRCAATRCCTCRCCCTRCTTCRRCTCRCRCACRRTRCRCACACACACTRRCCTRCRCCCACTRTCTRRCACTCCCTARTRARATRAACCC" \
	-o "${wd}/${prefix}.discordant_R2_HQ.fastq.gz" \
	"${wd}/${prefix}.discordant_R2.fastq.gz"

# extract read names for discordant HQ R2 reads
zcat "${wd}/${prefix}.discordant_R2_HQ.fastq.gz" \
| awk 'NR%4==1 {print substr($1,2)}' \
> "${wd}/${prefix}.discordant_R2_HQ.reads"

# keep only R1 singletons corresponding to these reads
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar "${picard}" FilterSamReads \
	--INPUT "${wd}/${prefix}.discordant_PE.R1to${reference}.bam" \
	--OUTPUT "${wd}/${prefix}.discordant_PE_HQ.R1to${reference}.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
	--READ_LIST_FILE "${wd}/${prefix}.discordant_R2_HQ.reads" \
	--FILTER includeReadList \
	--WRITE_READS_FILES FALSE \
	--USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
	--SORT_ORDER coordinate \
	--CREATE_INDEX TRUE

# transform HQ R2 fastq into unmapped bam file ------------------------------------
	java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} FastqToSam \
  -FASTQ "${wd}/${prefix}.discordant_R2_HQ.fastq.gz" \
  --OUTPUT "${wd}/${prefix}.discordant_R2_HQ.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
  --SAMPLE_NAME "${wd}/${prefix}.discordant_R2_HQ.bam" \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --SORT_ORDER queryname

# merge HQ R1 and R2 bam files for discordant PE ----------------------------------
samtools merge - \
  "${wd}/${prefix}.discordant_PE_HQ.R1to${reference}.bam" \
  "${wd}/${prefix}.discordant_R2_HQ.bam" \
| samtools sort -n \
| samtools view \
| awk -v OFS="\t" '
    (NR%2==1) {
      r1=$2";"$3";"$4
      if ($2==0) {
        printf $1 "\t73\t" $3 "\t" $4 "\t" "50" "\t" $6 "\t=\t" $4 "\t";
        for (i=9; i<NF; i++){printf $i "\t"};
        printf $NF "\n";
      }
      else if ($2==16){
        printf $1 "\t89\t" $3 "\t" $4 "\t" "50" "\t" $6 "\t=\t" $4 "\t";
        for (i=9; i<NF; i++){printf $i "\t"};
        printf $NF "\n"
      }
    }
    (NR%2==0) {
      split(r1, a, ";");
      if (a[1]==0) {
        printf $1 "\t133\t" a[2] "\t" a[3] "\t" $5 "\t" $6 "\t=\t" a[3] "\t";
        for (i=9; i<NF; i++){printf $i "\t"};
        printf $NF "\n"
      }
      else if (a[1]==16) {
        printf $1 "\t165\t" a[2] "\t" a[3] "\t" $5 "\t" $6 "\t=\t" a[3] "\t";
        for (i=9; i<NF; i++){printf $i "\t"};
        printf $NF "\n"
      }
    }
  ' \
| sed "s/\tRG:Z:[^\t]*//" \
> "${wd}/${prefix}.discordant_PE.${reference}.sam"

# rebuild header for merged sam file -------------------------------------------
samtools merge - \
  "${wd}/${prefix}.discordant_PE_HQ.R1to${reference}.bam" \
  "${wd}/${prefix}.discordant_R2_HQ.bam" \
| samtools view -H - \
| grep -v "^@RG" \
| grep -v "^@PG" \
> "${wd}/header.txt"

# prepend header to merged sam file and convert to bam -------------------------
cat "${wd}/header.txt" "${wd}/${prefix}.discordant_PE.${reference}.sam" \
| samtools view -b - \
| samtools sort -@ "${threads}" - \
> "${wd}/${prefix}.discordant_PE.${reference}.bam"

##### pool L1 reads mapped as proper or discordant pairs and deduplicate #######

# merge ref and non-ref bs-atlas-seq reads in bam file -------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} MergeSamFiles \
  --INPUT "${wd}/${prefix}.ref-L1.${reference}.bam" \
  --INPUT "${wd}/${prefix}.discordant_PE.${reference}.bam" \
  --OUTPUT "${wd}/${prefix}.merged.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --SORT_ORDER queryname

# remove duplicate reads in merged file ----------------------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} MarkDuplicates \
  --INPUT "${wd}/${prefix}.merged.bam" \
  --OUTPUT "${wd}/${prefix}.merged.dedup.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
  --METRICS_FILE "${wd}/${prefix}.merged.dedup.log" \
  --REMOVE_DUPLICATES TRUE \
  --ASSUME_SORT_ORDER queryname \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES

# sort deduplicated merged file ------------------------------------------------
java -Xmx8g -Dpicard.useLegacyParser=FALSE \
	-jar ${picard} SortSam \
  --INPUT "${wd}/${prefix}.merged.dedup.bam" \
  --OUTPUT "${wd}/${prefix}.merged.dedup.sorted.bam" \
	--VERBOSITY ERROR \
	--QUIET TRUE \
  --SORT_ORDER coordinate \
  --USE_JDK_DEFLATER TRUE \
  --USE_JDK_INFLATER TRUE \
  --CREATE_INDEX TRUE

# read count for report --------------------------------------------------------
count_mapped_proper=$( samtools view -@ "${threads}" -f 67 "${wd}/${prefix}.merged.bam" | wc -l )
count_mapped_discordant=$( samtools view -@ "${threads}" -f 73 "${wd}/${prefix}.merged.bam" | wc -l )
count_mapped_total=$(( "${count_mapped_proper}" + "${count_mapped_discordant}"))
count_mapped_dedup=$( samtools view -@ "${threads}" -f 65 "${wd}/${prefix}.merged.dedup.sorted.bam" | wc -l )

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

##### call L1 insertions #######################################################

output=$( printf "%-3s- Call reference L1 elements..." "${step}" )
log+="${output}"
echo -ne "${output}"

# ref-L1 elements (using only R2) ----------------------------------------------
# name=<chr>:<start-end 0-based>:<strand>:<Family>:<REF>:<% integrity>:<score>
# score=deduplicated read count
samtools view -b -f 131 "${wd}/${prefix}.merged.dedup.bam" \
| bedtools bamtobed -i - \
| sort -k1,1 -k2,2n \
| bedtools merge -s -c 4,5,6 -o distinct,count,distinct -i - \
| bedtools intersect -S -F 0.7 -wa -wb -a "${rmsk}" -b - \
| awk -v OFS="\t" '{
		split($4,family,"|");
		print $1,$2,$3,$1":"$2"-"$3":"$6":"family[3]":REF:"$5":"$(NF-1),$(NF-1),$6
	}' \
> "${wd}/${prefix}.called-refL1.bed"

# ref-L1 filtering based on user defined threshold, encode blacklisted regions,
# and extra-chromosomes --------------------------------------------------------
awk -v t="${read_threshold}" '$5>=t' "${wd}/${prefix}.called-refL1.bed" \
| awk '$1~/^chr[0-9]+|^chrX|^chrY/' \
| bedtools intersect -v -a - -b "${blacklist}" \
> "${wd}/${prefix}.called-refL1.filtered.bed"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

# non-ref L1 elements (using only R1) ------------------------------------------

output=$( printf "%-3s- Call non-reference L1 elements..." "${step}" )
log+="${output}"
echo -ne "${output}"

# calculate covered regions for non-refL1 insertions (using only R1) -----------
# score=deduplicated read count
# filters to call non-refL1:
#		- requires 5% of softclipped reads at 3' end (or at least 1 softclipped read
#			if count of reads <= 20)
#		- requires at least 60% of non-soft clipped reads

samtools view -b -f 73 "${wd}/${prefix}.merged.dedup.bam" \
| bedtools bamtobed -i - -cigar \
| sort -k1,1 -k2,2n \
| bedtools merge -s -d "${distance}" -c 4,5,6,7 -o collapse,count,distinct,collapse -i - \
| bedtools intersect -v -s -a - -b "${wd}/${prefix}.refL1.PE-clusters.bed" \
| bedtools intersect -s -wa -wb -a - -b "${wd}/${prefix}.split_reads.R1to${reference}.clusters.bed" \
| awk -v OFS="\t" '
		{
			softclip=0;
			n=split($7, cigar, ",");
			for (i=1; i<=n; i++) {
				if (($6=="+" && cigar[i]~/[0-9]+S$/)||($6=="-" && cigar[i]~/^[0-9]+S/)) {
					softclip++
				}
			}
			if ((($5<20 && softclip>=1) || (softclip>=5*n/100)) && softclip/n<=0.4) {
				print $1,$2,$3,$4,$5,$6,softclip,$12
			}
		}' \
> "${wd}/${prefix}.non-refL1.SE-clusters.bed"

# call non-ref L1 (precise integration site based on clusters) -----------------
# name=<chr>:<start-end 0-based>:<strand>:<Family>:<NONREF>
# 	:<assumed 100% integrity>:<total read count, softclipped read count, split read count>
# read counts are all related to deduplicated alignments
# score=deduplicated read count

awk -v OFS="\t" '
    ($6=="+") {print $1,$3-1,$3+1,$1":"$3-1"-"$3+1":"$6":L1HS:NONREF:100:"$5","$7","$8,$5,$6}
    ($6=="-") {print $1,$2-1,$2+1,$1":"$2-1"-"$2+1":"$6":L1HS:NONREF:100:"$5","$7","$8,$5,$6}
  ' "${wd}/${prefix}.non-refL1.SE-clusters.bed" \
> "${wd}/${prefix}.called-nonrefL1.bed"

# filter non-ref L1 insertions based on -n and -s, and exclude blacklisted regions,
# as well as extra-chromosomes -------------------------------------------------
awk -v n="${read_threshold}" -v s="${splitread_threshold}" -v OFS="\t" '
		($5>=n) {
			split($4, name, ":");
			split(name[7],count,",");
			if (count[3]>=s) {print};
		}' "${wd}/${prefix}.called-nonrefL1.bed" \
| awk '$1~/^chr[0-9]+|^chrX|^chrY/' \
| bedtools intersect -v -a - -b "${blacklist}" \
> "${wd}/${prefix}.called-nonrefL1.filtered.bed"

# insertion count for report ---------------------------------------------------
count_refL1_raw=$( grep -c L1 "${wd}/${prefix}.called-refL1.bed" )
count_refL1_filtered=$( grep -c L1 "${wd}/${prefix}.called-refL1.filtered.bed" )
count_refL1HS_raw=$( grep -c L1HS "${wd}/${prefix}.called-refL1.bed" )
count_refL1HS_filtered=$( grep -c L1HS "${wd}/${prefix}.called-refL1.filtered.bed" )
count_nonrefL1HS_raw=$( wc -l "${wd}/${prefix}.called-nonrefL1.bed" | awk '{print $1}' )
count_nonrefL1HS_filtered=$( wc -l "${wd}/${prefix}.called-nonrefL1.filtered.bed" | awk '{print $1}' )
count_L1_raw="$(( ${count_refL1_raw} + ${count_nonrefL1HS_raw} ))"
count_L1_filtered="$(( ${count_refL1_filtered} + ${count_nonrefL1HS_filtered} ))"
count_L1HS_raw="$(( ${count_refL1HS_raw}+${count_nonrefL1HS_raw} ))"
count_L1HS_filtered="$(( ${count_refL1HS_filtered}+${count_nonrefL1HS_filtered} ))"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

##### DNA methylation calling ##################################################

# ref L1 -----------------------------------------------------------------------

output=$( printf "%-3s- Call CpG methylation for reference L1 elements..." "${step}" )
log+="${output}"
echo -ne "${output}"

# list ref-L1 read pairs from merged dedup bam file ----------------------------
awk -v OFS="\t" '
		($6=="+"){print $1,$2-1,$2+1,$4,$5,$6}
		($6=="-"){print $1,$3-1,$3+1,$4,$5,$6}
	' "${wd}/${prefix}.called-refL1.filtered.bed" \
| bedtools intersect -s -wa -a "${wd}/${prefix}.refL1.PE-clusters.bed" -b - \
| samtools view -b "${wd}/${prefix}.merged.dedup.bam" -L - \
| samtools sort -n \
| samtools fixmate -O bam - - \
| samtools view -b -f 3 - \
> "${wd}/${prefix}.refL1.dedup.name_sorted.bam"

# extract methylation for ref L1 -----------------------------------------------
bismark_methylation_extractor \
--paired-end \
--comprehensive \
--merge_non_CpG \
--no_overlap \
--multicore "${bs_parallel}" \
--mbias_off \
--no_header \
--output "${wd}" \
--bedGraph \
--zero_based \
--cutoff "${read_threshold}" \
--buffer_size 50% \
"${wd}/${prefix}.refL1.dedup.name_sorted.bam" \
&> /dev/null


# prepare methpat .tsv file for ref-L1 called by bs-atlas-seq ------------------
bedtools intersect -s \
	-a "${wd}/${prefix}.called-refL1.filtered.bed" \
	-b "${wd}/${prefix}.refL1.PE-clusters.bed" \
| awk -v OFS="\t" '
    ($6=="+") {print $1,$2,$3,$4,$3-$2,0,24}
    ($6=="-") {print $1,$2,$3,$4,$3-$2,24,0}
  ' \
> "${wd}/${prefix}.called-refL1.filtered.tsv"

grep L1HS "${wd}/${prefix}.called-refL1.filtered.tsv" \
> "${wd}/${prefix}.called-refL1HS.filtered.tsv"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

# non-ref L1 -------------------------------------------------------------------

output=$( printf "%-3s- Call CpG methylation for non-reference L1 elements..." "${step}" )
log+="${output}"
echo -ne "${output}"

# create temporary folder to store read cluster individually -------------------
rm -f -r "${wd}/tmp/"
mkdir -p "${wd}/tmp"

# extract read names of each cluster (dedup) -----------------------------------
bedtools intersect -s -wa -wb \
	-a "${wd}/${prefix}.non-refL1.SE-clusters.bed" \
	-b "${wd}/${prefix}.called-nonrefL1.filtered.bed" \
| cut -f 1-8,12 \
> "${wd}/${prefix}.non-refL1.SE-clusters.filtered.bed"

while read aline;
do
	cluster=$( echo -ne "${aline}" | awk '{print $9}' ); \
	echo -e "${aline}" \
  | awk '{print $4}' \
  | awk -F "," '{gsub("/1",""); OFS="\n"; $NF=$NF; print $0}' \
  > "${wd}/tmp/${cluster}.reads" ; \
done < "${wd}/${prefix}.non-refL1.SE-clusters.filtered.bed"

# reformat read names in initial bam file (mapping to L1HS) to keep only
# the read base name -----------------------------------------------------------
samtools view -h "${wd}/${prefix}.R2toL1HS.bam" \
| awk -v OFS="\t" '
  ($1~/^@/) {print}
  ($1!~/^@/) {
    split($1,l,"_");
    printf l[1];
    for (i=2;i<NF;i++){printf "\t" $i};
    printf "\t" $NF "\n";
  }
' \
| samtools view -b - \
| samtools sort -@ "${threads}" - \
> "${wd}/${prefix}.R2toL1HS.reformated.bam"

# create one bam file for each cluster to call methylation separately ----------
read_extraction(){
	file="$1"
	picard="$2"
	wd="$3"
	prefix="$4"
	threads="$5"

	file_name="$( basename "${file}" .reads )"

  # extract aligned reads on L1HS corresponding to each cluster ----------------
  java -Dpicard.useLegacyParser=FALSE \
		-jar "${picard}" FilterSamReads \
    --INPUT "${wd}/${prefix}.R2toL1HS.reformated.bam" \
    --OUTPUT "${wd}/tmp/${file_name}.bam" \
		--VERBOSITY ERROR \
		--QUIET TRUE \
    --READ_LIST_FILE "${wd}/tmp/${file_name}.reads" \
    --FILTER includeReadList \
    --WRITE_READS_FILES FALSE \
		--USE_JDK_DEFLATER TRUE \
	  --USE_JDK_INFLATER TRUE \
    --SORT_ORDER coordinate

  # reformat bam file (append cluster coordinates to read names) ---------------
  samtools view -h "${wd}/tmp/${file_name}.bam" \
  | awk -v OFS="\t" -v f="${file_name}" '
    ($1~/^@/) {print}
    ($1!~/^@/) {
      printf $1 "|" f;
      for (i=2;i<NF;i++){printf "\t" $i};
      printf "\t" $NF "\n";
    }
  ' \
  | samtools view -@ "${threads}" -b - \
  | samtools sort -@ "${threads}" - \
  > "${wd}/tmp/${file_name}.reformated.bam"
}

export -f read_extraction

find "${wd}/tmp" -name "chr*.reads" | parallel --jobs "${threads}" read_extraction \
	{} \
	"${picard}" \
	"${wd}" \
	"${prefix}" \
	"${threads}" ;

# merge all cluster bam files into a single file
# (cluster coordinates appended to read names) ---------------------------------
samtools merge -@ "${threads}" "${wd}/${prefix}.R2toL1HS.nonref-L1HS.filtered.bam" "${wd}/tmp/"*.reformated.bam

# extract methylation of each mCG in each read of each cluster -----------------
bismark_methylation_extractor \
  --single-end \
  --comprehensive \
  --merge_non_CpG \
  --multicore "${bs_parallel}" \
  --mbias_off \
  --no_header \
  --output "${wd}" \
	--bedGraph \
	--zero_based \
	--cutoff "${read_threshold}" \
	--buffer_size 50% \
  "${wd}/${prefix}.R2toL1HS.nonref-L1HS.filtered.bam" \
	&> /dev/null

# reformat non-ref L1 bismark_methylation_extractor output to get insertion
# coordinates ------------------------------------------------------------------
awk -v OFS="\t" '
  {
    split($1,n,"|");
    print n[1],$2,n[2],$4,$5
  }' "${wd}/CpG_context_${prefix}.R2toL1HS.nonref-L1HS.filtered.txt" \
| sort -k3,3V -k1,1 -k4,4n \
> "${wd}/CpG_context_${prefix}.R2toL1HS.nonref-L1HS.filtered.renamed.txt"

# prepare methpat .tsv file for non-ref L1 called by bs-atlas-seq --------------
echo -ne "" > "${wd}/${prefix}.called-non-refL1HS.filtered.tsv"
awk -v OFS="\t" '{
    print $3, 1, 228, $3, 228, 0, 24
  }' "${wd}/CpG_context_${prefix}.R2toL1HS.nonref-L1HS.filtered.renamed.txt" \
| uniq \
>> "${wd}/${prefix}.called-non-refL1HS.filtered.tsv"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

####### final outputs for both ref- and non-ref L1 #############################

output=$( printf "%-3s- Calculate methylation patterns and entropy for each L1 locus..." "${step}" )
log+="${output}"
echo -ne "${output}"

# pool ref and non-ref L1 methylation calls and methpat .tsv files -------------
header="Bismark methylation extractor"

awk '($1!~/^Bismark/) {print}' \
  "${wd}/CpG_context_${prefix}.refL1.dedup.name_sorted.txt" \
  "${wd}/CpG_context_${prefix}.R2toL1HS.nonref-L1HS.filtered.renamed.txt" \
| sort -k3,3V -k1,1 -k4,4n \
| sed "1i ${header}" \
> "${wd}/CpG_context_${prefix}_pooled-L1.txt"

cat "${wd}/${prefix}.called-refL1HS.filtered.tsv" "${wd}/${prefix}.called-non-refL1HS.filtered.tsv" \
| sort -k4,4V \
> "${wd}/${prefix}.pooled-L1HS.tsv"

cat "${wd}/${prefix}.called-refL1.filtered.tsv" "${wd}/${prefix}.called-non-refL1HS.filtered.tsv" \
| sort -k4,4V \
> "${wd}/${prefix}.pooled-L1.tsv"

# calculate stat and organize display for pooled L1HS (no threshold for calculation)
methpat \
  --amplicons "${wd}/${prefix}.pooled-L1HS.tsv" \
  --logfile "${wd}/${prefix}.methpat_pooled-L1HS_CpG.log" \
  --html "${wd}/${prefix}.methpat_pooled-L1HS_CpG.html" \
  --webassets online \
  --title "${prefix} - L1HS family only" \
	--min_cpg_percent 10 \
  "${wd}/CpG_context_${prefix}_pooled-L1.txt" \
> "${wd}/${prefix}_methpat_pooled-L1HS_CpG.output"

# calculate stat and organize display for pooled L1 (no threshold for calculation)
methpat \
  --amplicons "${wd}/${prefix}.pooled-L1.tsv" \
  --logfile "${wd}/${prefix}.methpat_pooled-L1_CpG.log" \
  --html "${wd}/${prefix}.methpat_pooled-L1_CpG.html" \
  --webassets online \
  --title "${prefix} - all called L1 families" \
	--min_cpg_percent 10 \
  "${wd}/CpG_context_${prefix}_pooled-L1.txt" \
> "${wd}/${prefix}_methpat_pooled-L1_CpG.output"

# calculate average methylation and other stats for each L1 instance -----------
# ME=methylation entropy (from Xie et al. NAR 2011 and correction in 2013)
#   ME=0 -> homogenous pattern of call_methylation
#   ME=0.25 -> 50% for a given pattern and 50% for a second one
#   ME=1 -> all possible patterns observed

awk -v OFS="\t" 'BEGIN {print "#chr","start","end","name","mCG_frac","strand","family","ref","read_count","mCG_count","CG_count","ME"}' \
> "${wd}/${prefix}_methpat_pooled-L1_perCpG.txt"

awk '
  BEGIN {OFS="\t"}
  ($1!~/^</) {
    gsub("-","",$5);
    n=split($5,m,"");
    for (i=1;i<=n;i++){mCG[$1]+=m[i]*$6};
    pattern[$1","$5]=$6;
    read_count[$1]+=$6;
    CG[$1]+=n*$6;
  }
  END {
    for (id in CG) {
      split(id,a,":");
      split(a[2],coord,"-");
      for (i in pattern) {
        if (index(i,id)) {
          k[id]+=(-pattern[i]/read_count[id])*log(pattern[i]/read_count[id])/log(10)
        }
      };
      b[id]=CG[id]/read_count[id]
      ME[id]=(log(10)/(log(2)*b[id]))*k[id];
      print a[1],coord[1],coord[2],id,mCG[id]/CG[id],a[3],a[4],a[5],read_count[id],mCG[id],CG[id],ME[id];
    }
  }' "${wd}/${prefix}_methpat_pooled-L1_CpG.output" \
| sort -k1,1 -k2,2n \
>> "${wd}/${prefix}_methpat_pooled-L1_perCpG.txt"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

# extract methylation of flanking regions --------------------------------------

output=$( printf "%-3s- Generate bedGraph files for CpG upstream of called L1 loci..." "${step}" )
log+="${output}"
echo -ne "${output}"

# 	- for refL1
zcat "${wd}/${prefix}.refL1.dedup.name_sorted.bedGraph.gz" \
| bedtools intersect -v -a - -b "${wd}/${prefix}.called-refL1.filtered.bed" \
| gzip \
> "${wd}/${prefix}.refL1.flank.bedGraph.gz"

# 	- for non-refL1
samtools view -b -f 73 "${wd}/${prefix}.merged.dedup.bam" \
| bedtools intersect -s -a - -b "${wd}/${prefix}.non-refL1.SE-clusters.filtered.bed" \
| samtools view -h \
| awk -v OFS="\t" '
	($1~/^@/) {print}
	($1!~/@/) {
		printf $1"\t";
		if ($2==73) {printf "0\t"} else if ($2==89) {printf "16\t"};
		for (i=3;i<NF;i++){printf $i"\t"};
		printf $NF "\n";
	}' \
| samtools view -b \
> "${wd}/${prefix}.non-refL1.flank.bam"

bismark_methylation_extractor \
  --single-end \
  --comprehensive \
  --merge_non_CpG \
  --multicore "${bs_parallel}" \
  --mbias_off \
  --no_header \
  --output "${wd}" \
	--bedGraph \
	--zero_based \
	--cutoff "${read_threshold}" \
	--buffer_size 50% \
"${wd}/${prefix}.non-refL1.flank.bam" \
&> /dev/null

if [[ ! -f "${wd}/${prefix}.non-refL1.flank.bedGraph.gz" ]];
then
	echo -e "track type=bedGraph" \
	| gzip \
	> "${wd}/${prefix}.non-refL1.flank.bedGraph.gz"
fi

#		- merge flanking met of flanking sequence for ref and non-ref L1
zcat "${wd}/${prefix}.refL1.flank.bedGraph.gz" "${wd}/${prefix}.non-refL1.flank.bedGraph.gz" \
| sort -k1,1 -k2,2n \
| awk 'BEGIN {print "track type=bedGraph"} ($1!~/^track/) {print}' \
| gzip \
> "${wd}/${prefix}.pooled-L1.flank.bedGraph.gz"

#		- merge all bedgraph files for both ref and non-ref L1, internal and flanks
zcat "${wd}/${prefix}.refL1.dedup.name_sorted.bedGraph.gz" \
| bedtools intersect -a - -b "${wd}/${prefix}.called-refL1.filtered.bed" \
| sort -k1,1 -k2,2n \
| awk 'BEGIN {print "track type=bedGraph"} ($1!~/^track/) {print}' \
| gzip \
> "${wd}/${prefix}.refL1.bedGraph.gz"

zcat "${wd}/${prefix}.refL1.bedGraph.gz" "${wd}/${prefix}.refL1.flank.bedGraph.gz" "${wd}/${prefix}.non-refL1.flank.bedGraph.gz" \
| sort -k1,1 -k2,2n \
| awk 'BEGIN {print "track type=bedGraph"} ($1!~/^track/) {print}' \
| gzip \
> "${wd}/${prefix}.pooled-CG.bedGraph.gz"

output=$( printf "Done" )"\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

# Summarize and annotate called L1 elements ------------------------------------
output=$( printf "%-3s- Annotate identified L1 elements..." "${step}" )
log+="${output}"
echo -ne "${output}"

# find closest gene (+/- 10 kb) ------------------------------------------------
awk -v OFS="\t" 'BEGIN {print "#chr","start","end","name","mCG_frac","strand","family","ref","read_count","mCG_count","CG_count","ME","closest_gene"}' \
> "${wd}/${prefix}_final_summary_L1.bed"

sort -k1,1 -k2,2n "${genes}" \
| bedtools closest -d -t first \
  -a "${wd}/${prefix}_methpat_pooled-L1_perCpG.txt" \
  -b - \
| awk -F "\t" -v OFS="\t" '
  ($1!~/^#/) {
    if ($NF<=10000 && $NF>=0) {
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16
    } else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"."}
  }' \
>> "${wd}/${prefix}_final_summary_L1.bed"

# subset output for L1HS family only -------------------------------------------
awk '$1~/^#/ || $4~/:L1HS:/' "${wd}/${prefix}_final_summary_L1.bed" \
> "${wd}/${prefix}_final_summary_L1HS.bed"

# count of insertions with methylation call ------------------------------------
count_L1_met=$( awk '$1!~/^#/' "${wd}/${prefix}_final_summary_L1.bed" | wc -l )
count_L1HS_met=$( grep -c L1HS "${wd}/${prefix}_final_summary_L1HS.bed" )
count_nonrefL1HS_met=$( grep -c NONREF "${wd}/${prefix}_final_summary_L1HS.bed" )

output=$( printf "Done" )"\n${starline}\n"
log+="${output}"
echo -ne "${output}"
(( step++ ))

##### create report ############################################################
output="\
Sample\t${prefix}\n\
Unprocessed reads\t${count_raw}\n\
Trimmed reads\t${count_trimmed}\n\
Mapped to L1\t${count_mapped_1}\n\
Mapped to ${reference} as PE\t${count_mapped_proper}\n\
Mapped to ${reference} as SE\t${count_mapped_discordant}\n\
Total mapped to ${reference}\t${count_mapped_total}\n\
Deduplicated mapped reads\t${count_mapped_dedup}\n\
Called L1 unfiltered\t${count_L1_raw}\n\
- including L1HS\t${count_L1HS_raw}\n\
- including nonref L1HS\t${count_nonrefL1HS_raw}\n\
Called L1 filtered\t${count_L1_filtered}\n\
- including L1HS\t${count_L1HS_filtered}\n\
-- including nonref L1HS\t${count_nonrefL1HS_filtered}\n\
L1 with methylation call\t${count_L1_met}\n\
- including L1HS\t${count_L1HS_met}\n\
-- including nonref L1HS\t${count_nonrefL1HS_met}\n\
"
log+="${output}"
echo -ne "${output}" | tee "${wd}/${prefix}.count.txt"

output="${starline} \n"
log+="${output}"
echo -ne "${output}"

# calculate runtime for the whole pipeline -------------------------------------
end_time=$( date +%s )
runtime=$( date -u -d @$(( end_time - start_time )) +"%T" )
day=$(date +"[%d-%m-%Y] [%T]")

output="\
${day} \tRunning time: $runtime (hh:mm:ss) \n\
${starline} \n\
"
log+="${output}"
echo -ne "${output}"

##### cleanup ##################################################################

# move important files to output directory -------------------------------------
mv "${wd}/${prefix}.merged.dedup.sorted.bam" "${output_dir}/${prefix}.${reference}.bam"
mv "${wd}/${prefix}.merged.dedup.sorted.bai" "${output_dir}/${prefix}.${reference}.bai"
mv "${wd}/${prefix}.methpat_pooled-L1_CpG.html" "${output_dir}/${prefix}.all_L1.${reference}.html"
mv "${wd}/${prefix}.methpat_pooled-L1HS_CpG.html" "${output_dir}/${prefix}.L1HS.${reference}.html"
mv "${wd}/${prefix}_final_summary_L1.bed" "${output_dir}/${prefix}.all_L1.${reference}.bed"
mv "${wd}/${prefix}_final_summary_L1HS.bed" "${output_dir}/${prefix}.L1HS.${reference}.bed"
mv "${wd}/${prefix}.pooled-CG.bedGraph.gz" "${output_dir}/${prefix}.all_L1.flanks_and_internal.${reference}.bedGraph.gz"
mv "${wd}/${prefix}.pooled-L1.flank.bedGraph.gz" "${output_dir}/${prefix}.all_L1.flanks_only.${reference}.bedGraph.gz"
mv "${wd}/${prefix}.refL1.bedGraph.gz" "${output_dir}/${prefix}.refL1.internal.${reference}.bedGraph.gz"
mv "${wd}/${prefix}.count.txt" "${output_dir}/${prefix}.stats.txt"
echo -ne "${log}" > "${output_dir}/${prefix}.log"
cat "$0" > "${output_dir}/${prefix}.$( basename "$0" ).script.sh"

# remove temporary files -------------------------------------------------------
if [[ "${cleanup}" == "on" ]];
then
	rm -r "${wd}"
fi

# restore system variables -----------------------------------------------------
export LC_NUMERIC=$LC_NUMERIC_OLD
