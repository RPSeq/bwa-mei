#!/bin/sh

#Check if correct command line input

if [ $# -ne 6 ]; then
    echo "Usage: bash pipeline.sh <results prefix> <BAM file path> <MEI reference path> <output directory name> <ref_mapq_cutoff> <mei_mapq_cutoff>"
    exit 1
fi

#Get command line input
PREFIX=$1
INPUT_BAM=$2
MEI_REF=$3
OUTPUT_DIR_NAME=$4
REF_CUTOFF=$5
MEI_CUTOFF=$6


#Set relative dirs
WORKING_DIR=$(readlink -f $PWD)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FILE_LABEL=$(basename $INPUT_BAM ".bam")
FILE_DIR=$(dirname $INPUT_BAM)
SPLITTERS=${FILE_DIR}/${FILE_LABEL}.splitters.bam

#Set output directory
RESULTS_DIR=${WORKING_DIR}/${OUTPUT_DIR_NAME}

##Check if output dir exists and ask to overwrite
if [ -d ${RESULTS_DIR} ]; then

    read -p "Output directory already exists, do you want to overwrite? (yes/no): " overwrite

    if [ "$overwrite" == "y" ] || [ "$overwrite" == "yes" ] || \
    [ "$overwrite" == "Y" ] || [ "$overwrite" == "YES" ]; then
        #rm -f ${RESULTS_DIR}/*
        echo "OVERWRITTEN"
    else
        echo "Cancelled."
        exit 1
    fi

else
    mkdir ${RESULTS_DIR}
fi


# 1. get MEI candidates
sambamba view -t 2 -f bam -l 0 -F "\
paired and [MQ] != 0 and (mapping_quality == 0 or unmapped) \
and not (proper_pair or secondary_alignment \
or duplicate or mate_is_unmapped)" \
${INPUT_BAM} | bamToFastq -i /dev/stdin -fq /dev/stdout | bwa mem -t 8 -k 11 -M ${MEI_REF} - | \
sambamba view -t 4 --sam-input -f bam -F "not (secondary_alignment or duplicate or unmapped)" /dev/stdin | \
sambamba sort -t 4 -m 2G -n /dev/stdin -o ${RESULTS_DIR}/${PREFIX}.disc.fixed.bam;

# 2. get anchors for MEI candidates
sambamba view -t 4 -f bam -F "\
paired and ([MQ] == 0 or mate_is_unmapped) and mapping_quality > 0 \
and not (proper_pair or secondary_alignment \
or duplicate or unmapped)" \
${INPUT_BAM} > ${RESULTS_DIR}/${PREFIX}.disc.anchors.bam;

#print read counts to log
OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.disc.anchors.bam | wc -l)";
echo -e "Paired-end candidates:\t${OUTPUT}" > ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.disc.fixed.bam | wc -l)";
echo -e "Paired-end BWA MEM hits:\t${OUTPUT}" >> ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

#   a. Get IDS of mapped reads
samtools view ${RESULTS_DIR}/${PREFIX}.disc.fixed.bam | awk '{print $1}' | uniq > ${RESULTS_DIR}/${PREFIX}.disc.hits.id;

#   b. Pull hits from anchored bam file
python ${SCRIPTS_DIR}/remove_nonpairs.py -a ${RESULTS_DIR}/${PREFIX}.disc.anchors.bam \
-b ${RESULTS_DIR}/${PREFIX}.disc.hits.id | sambamba sort -t 4 -m 2G -n /dev/stdin -o ${RESULTS_DIR}/${PREFIX}.disc.anchored.mates.bam;

#   c. Merge bams from anchored and unanchored file.
sambamba merge -t 4 ${RESULTS_DIR}/${PREFIX}.disc.merged.bam ${RESULTS_DIR}/${PREFIX}.disc.anchored.mates.bam ${RESULTS_DIR}/${PREFIX}.disc.fixed.bam;

# 6. Fix the flags and mate coordinates for read pairs (currently excludes secondary alignments)
python $SCRIPTS_DIR/fix_pair_flags.py -i ${RESULTS_DIR}/${PREFIX}.disc.merged.bam | \
python ${SCRIPTS_DIR}/mapq_filter.py -S -a $REF_CUTOFF -u $MEI_CUTOFF > ${RESULTS_DIR}/${PREFIX}.disc.repaired.merged.bam;

# 7. Coordinate sort before running LUMPY
samtools sort ${RESULTS_DIR}/${PREFIX}.disc.repaired.merged.bam ${RESULTS_DIR}/${PREFIX}.disc.sorted.repaired.merged;

OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.disc.sorted.repaired.merged.bam | wc -l)";
OUTPUT=$((OUTPUT/2));
echo -e "Total read-pairs fed to LUMPY:\t${OUTPUT}" >> ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

#	b. Get insert size distribution file
samtools view ${INPUT_BAM} \
| head -n 1000000 \
| tail -n 100000 \
| python $SCRIPTS_DIR/pairend_distro.py -r 101 -X 4 -N 100000 -o ${RESULTS_DIR}/${PREFIX}.histo.out;

#extract and align splitters
python ${SCRIPTS_DIR}/extract_splitters.py -c 20 -i ${SPLITTERS} -a ${RESULTS_DIR}/${PREFIX}.split.anchors.bam -u /dev/stdout | bwa mem -t 8 -k 11 -M -C ${MEI_REF} - | \
sambamba view --sam-input -t 4 -h -F "not (unmapped or secondary_alignment)" /dev/stdin | \
sambamba view -t 1 --sam-input -f bam -l 0 /dev/stdin | sambamba sort -t 2 -m 2G -n /dev/stdin -o ${RESULTS_DIR}/${PREFIX}.split.hit.unanchors.bam;

OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.split.anchors.bam | wc -l)";
echo -e "Split-read candidates:\t${OUTPUT}" >> ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.split.hit.unanchors.bam | wc -l)";
echo -e "Split-read BWA MEM hits:\t${OUTPUT}" >> ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

#make list of successfully realigned unanchors
samtools view ${RESULTS_DIR}/${PREFIX}.split.hit.unanchors.bam | awk '{print $1}' > ${RESULTS_DIR}/${PREFIX}.split.hit.ids;

#pull corresponding anchors from anchors bam file
python ${SCRIPTS_DIR}/remove_nonpairs.py -a ${RESULTS_DIR}/${PREFIX}.split.anchors.bam -b ${RESULTS_DIR}/${PREFIX}.split.hit.ids | sambamba sort -t 4 -m 2G -n /dev/stdin -o ${RESULTS_DIR}/${PREFIX}.split.hit.anchors.bam;

#merge anchors and realigned unanchors
sambamba merge -t 4 ${RESULTS_DIR}/${PREFIX}.split.merged.bam ${RESULTS_DIR}/${PREFIX}.split.hit.anchors.bam ${RESULTS_DIR}/${PREFIX}.split.hit.unanchors.bam;

#call repair script
python ${SCRIPTS_DIR}/repair_cigars.py -c 20 -i ${RESULTS_DIR}/${PREFIX}.split.merged.bam | \
python ${SCRIPTS_DIR}/mapq_filter.py -S -a $REF_CUTOFF -u $MEI_CUTOFF > ${RESULTS_DIR}/${PREFIX}.split.repaired.bam;

OUTPUT="$(samtools view ${RESULTS_DIR}/${PREFIX}.split.repaired.bam | wc -l)";
OUTPUT=$((OUTPUT/2));
echo -e "Total split-read pairs fed to LUMPY:\t${OUTPUT}" >> ${RESULTS_DIR}/${PREFIX}_reads_flow.log;

sambamba sort -t 2 ${RESULTS_DIR}/${PREFIX}.split.repaired.bam -o ${RESULTS_DIR}/${PREFIX}.split.repaired.sorted.bam

#call lumpy
/gscmnt/gc2719/halllab/users/rsmith/git/lumpy-sv/bin/lumpy -mw 1 -tt 0 -pe bam_file:${RESULTS_DIR}/${PREFIX}.disc.sorted.repaired.merged.bam,histo_file:${RESULTS_DIR}/${PREFIX}.histo.out,mean:319.551326228,stdev:74.2952533362,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:0 -sr bam_file:${RESULTS_DIR}/${PREFIX}.split.repaired.sorted.bam,back_distance:10,min_mapping_threshold:0,weight:1,id:10,min_clip:20 > ${RESULTS_DIR}/${PREFIX}.vcf 2> ${RESULTS_DIR}/${PREFIX}.lumpyerr;

awk '{OFS="\t"; FS="\t"} {if($0 ~ /^E/ || $0 ~ /^A/ || $0 ~ /^C/) print $0}' ${RESULTS_DIR}/${PREFIX}.lumpyerr > ${RESULTS_DIR}/${PREFIX}.errsplits;