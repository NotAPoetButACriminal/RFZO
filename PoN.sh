#!/bin/bash
#
#SBATCH -J PoN
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 512G
#SBATCH --time 3-00:00:00

module load samtools
module load gatk

set -ex

### VARIABLES ###

# RUN VARIABLES
readarray -t SAMPLES < $1
PON=$2
WDIR="/lustre/imgge/RFZO"

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
INTERVALS="${WDIR}/refs/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_Combined_Mito.bed"

mkdir -p \
  output/${PON}/bams/metrics \
  output/${PON}/vcfs \
  output/${PON}/counts

### START ###

for SAMPLE in "${SAMPLES[@]}"
do
  gatk Mutect2 \
    -R ${REF} \
    -L ${INTERVALS} \
    -ip 10 \
    -I ${WDIR}/input/${SAMPLE}.bam \
    -O ${WDIR}/output/${PON}/vcfs/${SAMPLE}_normal.vcf.gz \
    --max-mnp-distance 0 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --native-pair-hmm-threads 1 &
done

wait

VCFS=$(ls ${WDIR}/output/${PON}/vcfs/*_normal.norm.vcf.gz | sed -e 's/^/ -V /g')

gatk GenomicsDBImport \
  -R ${REF} \
  -L ${INTERVALS} \
  --merge-input-intervals \
  $VCFS \
  --genomicsdb-workspace-path ${WDIR}/output/${PON}/gdb \
  --genomicsdb-shared-posixfs-optimizations true \
  --verbosity DEBUG

gatk CreateSomaticPanelOfNormals \
  -R ${REF} \
  --germline-resource ${WDIR}/refs/gatk_af-only-gnomad.hg38.vcf.gz \
  -V gendb://${WDIR}/output/${PON}/gdb \
  -O ${WDIR}/output/${PON}/PoN_SNV_${PON}.vcf.gz

# Remove multi-allelic sites because CreateSomaticPanelOfNormals falsely includes all of them with FRACTION=1.00
# See github issue: https://github.com/broadinstitute/gatk/issues/8916
# Ended up not removing them after all...

# zcat ${WDIR}/output/${PON}/PoN_SNV_${PON}.vcf.gz \
#  | awk '$5 !~ /,/' \
#  | bgzip -k \
#  > ${WDIR}/output/${PON}/PoN_SNV_${PON}_fixed.vcf.gz \
#  && tabix ${WDIR}/output/${PON}/PoN_SNV_${PON}_fixed.vcf.gz

echo "Done creating PoN for SNVs"

gatk PreprocessIntervals \
  -R ${REF} \
  -L ${INTERVALS} \
  -O ${WDIR}/output/${PON}/cnv_bins.interval_list \
  --bin-length 0 \
  -imr OVERLAPPING_ONLY \
  --padding 250

sed -i '$d' ${WDIR}/output/${PON}/cnv_bins.interval_list

for SAMPLE in "${SAMPLES[@]}"
do
  gatk CollectReadCounts \
    -R ${REF} \
    -L ${WDIR}/output/${PON}/cnv_bins.interval_list \
    -imr OVERLAPPING_ONLY \
    -I ${WDIR}/input/${SAMPLE}.bam \
    -O ${WDIR}/output/${PON}/counts/${SAMPLE}.hdf5 &
done

wait

HDFS=$(ls ${WDIR}/output/${PON}/counts/*.hdf5 | sed -e 's/^/ -I /g')

gatk CreateReadCountPanelOfNormals \
  $HDFS \
  -O ${WDIR}/output/${PON}/PoN_CNV_${PON}.hdf5

echo "Done creating PoN for CNVs"