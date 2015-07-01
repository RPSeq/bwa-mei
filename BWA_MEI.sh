#!/bin/sh

#Check if correct command line input
if [ $# -ne 6 ]; then
    echo "Usage: bash pipeline.sh <results prefix> <BAM file path> <MEI reference path> <output directory name> <REF MAPQ cutoff> <MEI MAPQ cutoff>"
    exit 1
fi

#Get command line input
PREFIX=$1
INPUT_BAM=$2
MEI_REF=$3
OUTPUT_DIR_NAME=$4
REF_MAPQ=$5
MEI_MAPQ=$6

#Set relative dirs
WORKING_DIR=$(readlink -f $PWD)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FILE_LABEL=$(basename $INPUT_BAM ".bam")
FILE_DIR=$(dirname $INPUT_BAM)

#Set output directory
RESULTS_DIR=${WORKING_DIR}/${OUTPUT_DIR_NAME}
OUTPUT=$RESULTS_DIR/$PREFIX

##Check if output dir exists and ask to overwrite
if [ -d ${RESULTS_DIR} ]; then
    echo "Error: output directory already exists."
    exit
else
    mkdir ${RESULTS_DIR}
fi

#Pipes for streaming candidate alignments
mkfifo ${OUTPUT}.pairs_fq
mkfifo ${OUTPUT}.clips_fq

#	Note- These three commands run simultaneously. 
#	1. Extract MEI realignment candidates from input bam file
sambamba view -t 4 -f bam -l 0 $INPUT_BAM \
	| python ${SCRIPTS_DIR}/extract_candidates.py -pa ${OUTPUT}.pair_anchors.sam -pf ${OUTPUT}.pairs_fq -ca ${OUTPUT}.clip_anchors.sam -cf ${OUTPUT}.clips_fq -c 20 &

# 		a. Realign the paired candidates
bwa mem -t 8 -k 11 -M $MEI_REF ${OUTPUT}.pairs_fq \
	| samtools view -F 260 -b -u - \
	| sambamba sort -t 2 -m 2G -n -o ${OUTPUT}.pairs.bam /dev/stdin &

# 		b. Realign the clipped candidates
bwa mem -t 8 -k 11 -M $MEI_REF ${OUTPUT}.clips_fq \
	| samtools view -F 260 -b -u - \
	| sambamba sort -t 2 -m 2G -n -o ${OUTPUT}.clips.bam /dev/stdin &

#Wait for extraction and realignment to finish and delete pipes
wait
rm -f ${OUTPUT}.pairs_fq
rm -f ${OUTPUT}.clips_fq

#	2. Use zjoin to get clip anchors that match with clip unanchors that hit MEIs.
#		Then, namesort and merge for repair_cigars script, followed by MAPQ filtering.
#		Coordinate sort the resulting BAM. This is the SR file for LUMPY.
zjoin -a ${OUTPUT}.clip_anchors.sam -b <(samtools view ${OUTPUT}.clips.bam) -wa \
	| cat <(grep "^@" ${OUTPUT}.clip_anchors.sam) - \
	| samtools view -b -u - \
	| sambamba sort -t 2 -m 3G -n -o /dev/stdout /dev/stdin \
	| sambamba merge /dev/stdout ${OUTPUT}.clips.bam /dev/stdin \
	| python ${SCRIPTS_DIR}/repair_cigars.py -i /dev/stdin -c 20 \
	| python ${SCRIPTS_DIR}/mapq_filter.py -S -a $REF_MAPQ -u $MEI_MAPQ \
	| sambamba sort -m 3G -t 2 -o ${OUTPUT}.clips_merged.bam /dev/stdin

#Pipes for streaming paired read hits
mkfifo ${OUTPUT}.pair_anchor_hits
mkfifo ${OUTPUT}.pair_mei_hits

#	3. Filter the realigned pairs, namesort and merge the filtered MEI hits and anchored mates,
#		and fix the pair flags. Then filter by MAPQ and coordinate sort.
python ${SCRIPTS_DIR}/filter_pair_realigns.py -i ${OUTPUT}.pairs.bam -a ${OUTPUT}.pair_anchors.sam -o ${OUTPUT}.pair_mei_hits -ao ${OUTPUT}.pair_anchor_hits &
sambamba merge /dev/stdout <(sambamba sort -m 3G -t 2 -n -o /dev/stdout ${OUTPUT}.pair_mei_hits) <(sambamba sort -m 3G -t 2 -n -o /dev/stdout ${OUTPUT}.pair_anchor_hits) \
	| python ${SCRIPTS_DIR}/fix_pair_flags.py -i /dev/stdin \
	| python ${SCRIPTS_DIR}/mapq_filter.py -S -a $REF_MAPQ -u $MEI_MAPQ \
	| sambamba sort -m 3G -t 2 -o ${OUTPUT}.pairs_merged.bam /dev/stdin

#Wait for filtering script to finish then delete the pipes
wait
rm -f ${OUTPUT}.pair_anchor_hits
rm -f ${OUTPUT}.pair_mei_hits

# 	4. Get insert size distribution file
HIST_OUT=$(samtools view ${INPUT_BAM} | head -n 1000000 | tail -n 100000 | python $SCRIPTS_DIR/pairend_distro.py -r 101 -X 4 -N 100000 -o ${RESULTS_DIR}/${PREFIX}.histo.out)
MEAN=$(echo $HIST_OUT | tr ' ' '\n'  | cut -f 2 -d ":" | tr '\n' '\t' | cut -f 1)
STDEV=$(echo $HIST_OUT | tr ' ' '\n'  | cut -f 2 -d ":" | tr '\n' '\t' | cut -f 2)

#	5. Call lumpy
/gscmnt/gc2719/halllab/users/rsmith/git/lumpy-sv/bin/lumpy -mw 1 -tt 0 -pe bam_file:${OUTPUT}.pairs_merged.bam,histo_file:${OUTPUT}.histo.out,mean:${MEAN},stdev:${STDEV},read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:0 -sr bam_file:${OUTPUT}.clips_merged.bam,back_distance:10,min_mapping_threshold:0,weight:1,id:10,min_clip:20 > ${OUTPUT}.vcf 2> ${OUTPUT}.lumpyerr;

#	6. Cluster LUMPY breakpoint calls and collapse into indivual MEI events.
python ${SCRIPTS_DIR}/cluster.py -i <(vcf-sort ${OUTPUT}.vcf) > ${OUTPUT}.clusters.bed