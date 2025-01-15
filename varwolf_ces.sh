#!/bin/bash
#
#SBATCH -J varwolf_ces
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 100
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

module load fastp
module load bwa
module load samtools
module load gatk
module load miniconda3

set -eux

### VARIABLES ###

# RUN VARIABLES
readarray -t SAMPLE < $1
COHORT=$2
WDIR="/lustre/imgge/RFZO"
THREADS=${SLURM_CPUS_PER_TASK}

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
INTERVALS="${WDIR}/refs/TSOne_Expanded_Final_TargetedRegions_v2_hg38.bed"

mkdir -p output/${COHORT}/bams output/${COHORT}/vcfs

### START ###

for i in $(seq 0 $((${#SAMPLE[@]} -1)))
do
    for LANE in {1..4}
    do
        FLOWCELL=$(zcat "${WDIR}"/fastq/${SAMPLE[${i}]}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 3)
        LIBRARY=$(zcat ${WDIR}/fastq/${SAMPLE[${i}]}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
        fastp \
		    -w 1 \
	        -i ${WDIR}/fastq/${SAMPLE[${i}]}*L00${LANE}_R1*.fastq.gz \
            -I ${WDIR}/fastq/${SAMPLE[${i}]}*L00${LANE}_R2*.fastq.gz \
            --stdout \
        | bwa mem \
		    -t 2 \
		    -M -p -R "@RG\tID:${FLOWCELL}.LANE${LANE}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE[${i}]}" \
            ${REF} - \
        | samtools sort \
		    -@ 1 -n \
            > ${WDIR}/output/${COHORT}/bams/${SAMPLE[${i}]}_L${LANE}.bam
    done &
done

wait

echo "Finished aligning samples"

gatk MarkDuplicatesSpark \
	-R ${REF} \
	${BAMSHARDS} \
	-O ${WDIR}/bams/${SAMPLE}_dd.bam \
	-M ${WDIR}/bams/metrics/${SAMPLE}_mdmetrics.txt \
	--spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished writing ${SAMPLE}"

gatk BQSRPipelineSpark \
	-R ${REF} \
	-I ${WDIR}/bams/${SAMPLE}_dd.bam \
	-O ${WDIR}/bams/${SAMPLE}.bam \
	--known-sites ${DBSNP} \
	--spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished recalibrating ${SAMPLE}"

gatk CollectReadCounts \
	-R ${REF} \
    -L ${WDIR}/refs/read_counts_ces.interval_list \
    -imr OVERLAPPING_ONLY \
    -I ${WDIR}/bams/${SAMPLE}.bam \
    -O ${WDIR}/counts/${SAMPLE}.hdf5

gatk HaplotypeCaller \
    -R ${REF} \
    -L ${INTERVALS} \
    -ip 10 \
    -I ${WDIR}/bams/${SAMPLE}.bam \
    -O ${WDIR}/vcfs/${SAMPLE}_raw.vcf.gz

gatk VariantFiltration \
        -V ${WDIR}/vcfs/${SAMPLE}_raw.vcf.gz \
        -filter "DP < 5.0" --filter-name "DP5" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${WDIR}/vcfs/${SAMPLE}.vcf

sed -i 's/##source=HaplotypeCaller/##source=HaplotypeCaller\n##reference=hg38.fasta/g' ${WDIR}/vcfs/${SAMPLE}.vcf

bgzip -@ 16 ${WDIR}/vcfs/${SAMPLE}.vcf && tabix ${WDIR}/vcfs/${SAMPLE}.vcf.gz

rm ${WDIR}/bams/*${SAMPLE}*_* ${WDIR}/vcfs/${SAMPLE}_* ${WDIR}/vcfs/${SAMPLE}*idx
