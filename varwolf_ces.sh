#!/bin/bash
#
#SBATCH -J varwolf_ces
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 100
#SBATCH --mem 256Gsq
#SBATCH --time 3-00:00:00

module load fastp
module load bwa
module load samtools
module load gatk
module load miniconda3

set -ex

### VARIABLES ###

# RUN VARIABLES
readarray -t SAMPLES < $1
COHORT=$2
WDIR="/lustre/imgge/RFZO"
THREADS=${SLURM_CPUS_PER_TASK}

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
INTERVALS="${WDIR}/refs/TruSight_One_TargetedRegions_v1.1.hg38.bed"

mkdir -p \
    output/${COHORT}/bams/metrics \
    output/${COHORT}/vcfs \
    output/${COHORT}/counts

### START ###

# for SAMPLE in "${SAMPLES[@]}"
# do
#     for LANE in {1..4}
#     do
#         FLOWCELL=$(zcat "${WDIR}"/fastq/${SAMPLE}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 3)
#         LIBRARY=$(zcat ${WDIR}/fastq/${SAMPLE}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
#         fastp \
# 		    -w 1 \
# 	        -i ${WDIR}/fastq/${SAMPLE}*L00${LANE}_R1*.fastq.gz \
#             -I ${WDIR}/fastq/${SAMPLE}*L00${LANE}_R2*.fastq.gz \
#             --stdout \
#             -j ${WDIR}/output/${COHORT}/fastp.json \
#             -h ${WDIR}/output/${COHORT}/fastp.html \
#         | bwa mem \
# 		    -t 2 \
# 		    -M -p \
#             -R "@RG\tID:${FLOWCELL}.LANE${LANE}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE}" \
#             ${REF} - \
#         | samtools sort \
# 		    -@ 1 \
#             > ${WDIR}/output/${COHORT}/bams/${SAMPLE}_L${LANE}.bam
#     done &
# done

# wait

# echo "Finished mapping reads!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     BAMSHARDS=$(ls ${WDIR}/output/${COHORT}/bams/${SAMPLE}_L* | sed -e 's/^/ -I /g')
#     gatk MarkDuplicates \
# 	    -R ${REF} \
# 	    ${BAMSHARDS} \
# 	    -O ${WDIR}/output/${COHORT}/bams/${SAMPLE}_dd.bam \
# 	    -M ${WDIR}/output/${COHORT}/bams/metrics/${SAMPLE}_mdmetrics.txt &
# done

# wait

# echo "Finished marking duplicates!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     gatk BaseRecalibrator \
# 	    -R ${REF} \
# 	    -I ${WDIR}/output/${COHORT}/bams/${SAMPLE}_dd.bam \
# 	    -O ${WDIR}/output/${COHORT}/bams/metrics/${SAMPLE}_bqsr.report \
# 	    --known-sites ${DBSNP} &
# done

# wait

# for SAMPLE in "${SAMPLES[@]}"
# do
#     gatk ApplyBQSR \
# 	    -R ${REF} \
# 	    -I ${WDIR}/output/${COHORT}/bams/${SAMPLE}_dd.bam \
# 	    -O ${WDIR}/output/${COHORT}/bams/${SAMPLE}.bam \
#         --bqsr ${WDIR}/output/${COHORT}/bams/metrics/${SAMPLE}_bqsr.report &
# done

# wait

# echo "Finished recalibrating bases!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     gatk HaplotypeCaller \
# 	    -R ${REF} \
#         -L ${INTERVALS} \
#         -ip 10 \
# 	    -I ${WDIR}/output/${COHORT}/bams/${SAMPLE}.bam \
#         -O ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_raw.vcf.gz &
# done

# wait

# echo "Finished calling variants!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     gatk VariantFiltration \
# 	    -V ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_raw.vcf.gz \
#         -filter "DP < 5.0" --filter-name "DP5" \
#         -filter "QD < 2.0" --filter-name "QD2" \
#         -filter "QUAL < 30.0" --filter-name "QUAL30" \
#         -filter "SOR > 3.0" --filter-name "SOR3" \
#         -filter "FS > 60.0" --filter-name "FS60" \
#         -filter "MQ < 40.0" --filter-name "MQ40" \
#         -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#         -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#         -O ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_filtered.vcf.gz &
# done

# wait

# echo "Finished filtering variants!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     zcat ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_filtered.vcf.gz \
#     | sed 's/##source=HaplotypeCaller/##source=HaplotypeCaller\n##reference=hg38.fasta/g' \
#     | bgzip -@ 2 -o ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.vcf.gz
#     tabix ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.vcf.gz
# done

# echo "Done with SNVs!"

# for SAMPLE in "${SAMPLES[@]}"
# do
#     gatk CollectReadCounts \
#         -R ${REF} \
#         -L ${WDIR}/refs/read_counts_ces.interval_list \
#         -imr OVERLAPPING_ONLY \
#         -I ${WDIR}/output/${COHORT}/bams/${SAMPLE}.bam \
#         -O ${WDIR}/output/${COHORT}/counts/${SAMPLE}.hdf5 &
# done

# wait

# echo "Finished counting reads!"

eval "$(conda shell.bash hook)"
conda activate gatk

HDF5S=$(ls ${WDIR}/output/${COHORT}/counts/*.hdf5 | sed -e 's/^/ -I /g')

gatk FilterIntervals \
  -L ${WDIR}/refs/read_counts_ces.interval_list \
  --annotated-intervals ${WDIR}/refs/read_counts_ces_annotated.interval_list \
  -imr OVERLAPPING_ONLY \
  $HDF5S \
  -O ${WDIR}/output/${COHORT}/filtered.interval_list \
  --low-count-filter-percentage-of-samples 65

gatk IntervalListTools \
  -I ${WDIR}/output/${COHORT}/filtered.interval_list  \
  -O ${WDIR}/output/${COHORT}/interval_scatters \
  --SUBDIVISION_MODE INTERVAL_COUNT \
  --SCATTER_CONTENT 3000

SCATTERS=$(basename ${WDIR}/output/${COHORT}/interval_scatters/temp_0001_of_* \
  | cut -d "_" -f 4)

echo "Finished scattering intervals!"

gatk DetermineGermlineContigPloidy \
  -L ${WDIR}/output/${COHORT}/filtered.interval_list \
  -imr OVERLAPPING_ONLY \
  $HDF5S \
  -O ${WDIR}/output/${COHORT}/ \
  --output-prefix ploidy \
  --contig-ploidy-priors ${WDIR}/refs/contig_ploidy_priors.tsv

echo "Finished determining ploidy!"

for SCATTER in $(seq -w 0001 00${SCATTERS})
do
  gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L ${WDIR}/output/${COHORT}/interval_scatters/temp_${SCATTER}_of_${SCATTERS}/scattered.interval_list \
    --annotated-intervals ${WDIR}/refs/read_counts_ces_annotated.interval_list \
    -imr OVERLAPPING_ONLY \
    $HDF5S \
    -O ${WDIR}/output/${COHORT}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/counts/${COHORT}/ploidy-calls \
    --verbosity DEBUG &
done

wait

echo "Finished calling CNVs per scatter"

MODELS=$(ls -p ${WDIR}/output/${COHORT}/gcnvcaller_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/output/${COHORT}/gcnvcaller_scatters/#g")
CALLS=$(ls -p ${WDIR}/output/${COHORT}/gcnvcaller_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/output/${COHORT}/gcnvcaller_scatters/#g")

for i in $(seq 0 $((${#SAMPLES[@]} -1)))
do 
  gatk PostprocessGermlineCNVCalls \
    $MODELS \
    $CALLS \
    --sample-index ${i} \
    --output-genotyped-intervals ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_intervals.cnv.vcf.gz \
    --output-genotyped-segments ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_raw.cnv.vcf.gz \
    --output-denoised-copy-ratios ${WDIR}/output/${COHORT}/${SAMPLES[${i}]}_denoised_copy_ratios.tsv \
    --contig-ploidy-calls ${WDIR}/output/${COHORT}/ploidy-calls/ \
    --allosomal-contig chrX --allosomal-contig chrY \
    --sequence-dictionary ${WDIR}/refs/hg38.dict &
done

wait

echo "Finished calling CNVs per sample"

for SAMPLE in "${SAMPLES[@]}"
do
  gatk VariantFiltration \
    -V ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_raw.cnv.vcf.gz \
    -filter "QUAL < 30.0" \
    --filter-name "CNVQUAL" \
    -O ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_filtered.cnv.vcf.gz &
done

for SAMPLE in "${SAMPLES[@]}"
do
    zgrep -P -v "CNVQUAL|N\t\." ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
    | sed -e 's/##source=VariantFiltration/##source=VariantFiltration\n##reference=hg38.fasta/g' \
        -e 's/\tEND/\tSVTYPE=CNV;END/g' \
    | bgzip -o ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
    tabix ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
done

echo "Done with CNVs!"

rm ${WDIR}/output/${COHORT}/bams/*_* ${WDIR}/output/${COHORT}/vcfs/*_*

echo "All Done!"