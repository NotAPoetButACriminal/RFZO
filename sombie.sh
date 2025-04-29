#!/bin/bash
#
#SBATCH -J sombie
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
INTERVALS="${WDIR}/refs/avenio_hg38.bed"

mkdir -p \
  ${WDIR}/output/${SAMPLE}/bams/metrics \
  ${WDIR}/output/${SAMPLE}/vcfs/metrics \
  ${WDIR}/output/${SAMPLE}/counts

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
    -@ 4 -n \
  > ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_L${LANE}.bam
done

echo "Finished mapping reads!"

BAMSHARDS=$(ls ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_L* | sed -e 's/^/ -I /g')
gatk MarkDuplicatesSpark \
  -R ${REF} \
  ${BAMSHARDS} \
  -O ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -M ${WDIR}/output/${SAMPLE}/bams/metrics/mdmetrics.txt \
  --spark-runner LOCAL \
  --spark-master local[${THREADS}]

echo "Finished marking duplicates!"

gatk BQSRPipelineSpark \
  -R ${REF} \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -O ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  --known-sites /lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz \
  --spark-runner LOCAL \
  --spark-master local[${THREADS}]

echo "Finished recalibrating bases!"

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_avenio.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5

gatk GetPileupSummaries \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/vcfs/metrics/pileup.table \
  -V ${WDIR}/refs/gatk_af-only-gnomad.hg38.vcf.gz \
  -L ${WDIR}/refs/gatk_af-only-gnomad-common-biallelic.vcf.gz

gatk CalculateContamination \
  -I ${WDIR}/output/${SAMPLE}/vcfs/metrics/pileup.table \
  -O ${WDIR}/output/${SAMPLE}/vcfs/metrics/contamination.table \
  -segments ${WDIR}/output/${SAMPLE}/vcfs/metrics/segments.table

gatk Mutect2 \
  -R ${REF} \
  --panel-of-normals ${WDIR}/refs/gatk_1000g_pon.hg38.vcf.gz \
  --germline-resource ${WDIR}/refs/gatk_af-only-gnomad.hg38.vcf.gz \
  --af-of-alleles-not-in-resource 0.0000025 \
  -L ${INTERVALS} \
  -ip 100 \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz \
  --f1r2-tar-gz ${WDIR}/output/${SAMPLE}/vcfs/metrics/f1r2.tar.gz \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --native-pair-hmm-threads ${THREADS}

gatk LearnReadOrientationModel \
  -I ${WDIR}/output/${SAMPLE}/vcfs/metrics/f1r2.tar.gz \
  -O ${WDIR}/output/${SAMPLE}/vcfs/metrics/orientation_model.tar.gz

gatk FilterMutectCalls \
  -R ${REF} \
  -V ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}.vcf.gz \
  --contamination-table ${WDIR}/output/${SAMPLE}/vcfs/metrics/contamination.table \
  --tumor-segmentation ${WDIR}/output/${SAMPLE}/vcfs/metrics/segments.table \
  --ob-priors ${WDIR}/output/${SAMPLE}/vcfs/metrics/orientation_model.tar.gz

rm ${WDIR}/output/${SAMPLE}/bams/*_*