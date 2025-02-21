#!/bin/bash
#
#SBATCH -J varwolf_ces
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 3-00:00:00

module load fastp
module load bwa
module load samtools
module load gatk

set -ex

### VARIABLES ###

# RUN VARIABLES
SAMPLE=$1
WDIR="/lustre/imgge/RFZO"
THREADS=${SLURM_CPUS_PER_TASK}

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
INTERVALS="${WDIR}/refs/TruSight_One_TargetedRegions_v1.1.hg38.bed"

mkdir -p \
  output/${SAMPLE}/bams/metrics \
  output/${SAMPLE}/vcfs \
  output/${SAMPLE}/counts

### START ###

for LANE in {1..4}
do
  FLOWCELL=$(zcat "${WDIR}"/input/${SAMPLE}_*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 3)
  LIBRARY=$(zcat ${WDIR}/input/${SAMPLE}_*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
  fastp \
    -w 4 \
    -i ${WDIR}/input/${SAMPLE}_*L00${LANE}_R1*.fastq.gz \
    -I ${WDIR}/input/${SAMPLE}_*L00${LANE}_R2*.fastq.gz \
    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --stdout \
    -j ${WDIR}/output/${SAMPLE}/fastp.json \
    -h ${WDIR}/output/${SAMPLE}/fastp.html \
  | bwa mem \
    -t ${THREADS} \
    -M -p \
    -R "@RG\tID:${FLOWCELL}.LANE${LANE}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE}" \
    ${REF} - \
  | samtools sort \
    -@ 4 -n\
  > ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_L${LANE}.bam
done

echo "Finished mapping reads!"

BAMSHARDS=$(ls ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_L* | sed -e 's/^/ -I /g')
gatk MarkDuplicatesSpark \
  -R ${REF} \
  ${BAMSHARDS} \
  -O ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -M ${WDIR}/output/${SAMPLE}/bams/metrics/${SAMPLE}_mdmetrics.txt \
  --spark-runner LOCAL \
  --spark-master local[${THREADS}]

echo "Finished marking duplicates!"

gatk BQSRPipelineSpark \
  -R ${REF} \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -O ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  --known-sites ${DBSNP} \
  --spark-runner LOCAL \
  --spark-master local[${THREADS}]

echo "Finished recalibrating bases!"

gatk Mutect2 \
  -R ${REF} \
  --panel-of-normals ${WDIR}/refs/gatk_1000g_pon.hg38.vcf.gz \
  --germline-resource ${WDIR}/refs/gatk_af-only-gnomad.hg38.vcf.gz \
  --af-of-alleles-not-in-resource 0.0000025 \
  -L ${WDIR}/refs/sample3_depth_intervals.bed \
  -ip 50 \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz \
  --f1r2-tar-gz ${WDIR}/output/${SAMPLE}/vcfs/metrics/f1r2.tar.gz \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --native-pair-hmm-threads ${THREADS}

gatk FilterMutectCalls \
  -R ${REF} \
  -V ${WDIR}/output/SAMPLE5/vcfs/SAMPLE5_raw.vcf.gz \
  --contamination-table ${WDIR}/output/SAMPLE5/vcfs/metrics/SAMPLE5_contamination.table \
  --tumor-segmentation ${WDIR}/output/SAMPLE5/vcfs/metrics/SAMPLE5_segments.table \
  -O ${WDIR}/output/SAMPLE5/vcfs/SAMPLE5_filtered.vcf.gz

rm ${WDIR}/output/${SAMPLE}/bams/*_*
