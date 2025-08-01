#!/bin/bash
#
#SBATCH -J varwolf_ces_case
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 3-00:00:00

module load fastp
module load bwa
module load samtools
module load gatk
module load miniconda3
module load manta

set -ex

### VARIABLES ###

# RUN VARIABLES
SAMPLE=$1
COHORT=$2
WDIR="/lustre/imgge/RFZO"
THREADS=${SLURM_CPUS_PER_TASK}

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
INTERVALS="${WDIR}/refs/TruSight_One_TargetedRegions_v1.1.hg38.bed"

mkdir -p \
  output/${COHORT}/${SAMPLE}/bams/metrics \
  output/${COHORT}/${SAMPLE}/vcfs \
  output/${COHORT}/${SAMPLE}/counts

### START ###

for LANE in {1..4}
do
  (FLOWCELL=$(zcat "${WDIR}"/input/{SAMPLE}_*L001_R1*.fastq.gz \
    | head -1 \
    | cut -d ":" -f 3)
  LIBRARY=$(zcat ${WDIR}/input/${SAMPLE}_*L001_R1*.fastq.gz \
    | head -1 \
    | cut -d ":" -f 2 \
    | sed 's/^/Lib/g')
  fastp \
    -w 4 \
    -i ${WDIR}/input/${SAMPLE}_*L00${LANE}_R1*.fastq.gz \
    -I ${WDIR}/input/${SAMPLE}_*L00${LANE}_R2*.fastq.gz \
    --stdout \
    -j ${WDIR}/output/${COHORT}/${SAMPLE}/fastp.json \
    -h ${WDIR}/output/${COHORT}/${SAMPLE}/fastp.html \
  | bwa mem \
	  -t ${THREADS} \
	  -M -p \
    -R "@RG\tID:${FLOWCELL}.LANE${LANE}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE}" \
    ${REF} - \
  | samtools sort \
	  -@ 4 -n\
  > ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}_L${LANE}.bam)
done

echo "Finished mapping reads!"

BAMSHARDS=$(ls ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}_L* \
  | sed -e 's/^/ -I /g')
gatk MarkDuplicatesSpark \
  -R ${REF} \
  ${BAMSHARDS} \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -M ${WDIR}/output/${COHORT}/${SAMPLE}/bams/metrics/${SAMPLE}_mdmetrics.txt \
  --spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished marking duplicates!"

gatk BQSRPipelineSpark \
  -R ${REF} \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}_dd.bam \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}.bam \
  --known-sites ${DBSNP} \
  --spark-runner LOCAL \
	--spark-master local[${THREADS}]

echo "Finished recalibrating bases!"

gatk HaplotypeCaller \
  -R ${REF} \
  -L ${INTERVALS} \
  -ip 30 \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz
echo "Finished calling variants!"

gatk VariantFiltration \
  -V ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz \
  -filter "DP < 5.0" --filter-name "DP5" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}.vcf.gz

echo "Finished filtering variants!"

echo "Done with SNVs!"

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_ces.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5

echo "Finished counting reads!"

eval "$(conda shell.bash hook)"
conda activate gcnv

gatk DetermineGermlineContigPloidy \
  --model ${WDIR}/output/TSO250214/ploidy-model/ \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/ \
  --output-prefix ploidy

echo "Finished determining ploidy!"

gatk GermlineCNVCaller \
  --run-mode CASE \
  --model ${WDIR}/output/TSO250214/gcnvcaller/TSO250214-model \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller \
  --output-prefix ${SAMPLE} \
  --contig-ploidy-calls ${WDIR}/output/${COHORT}/${SAMPLE}/ploidy-calls \
  --verbosity DEBUG

echo "Finished calling CNVs"

gatk PostprocessGermlineCNVCalls \
  --model-shard-path ${WDIR}/output/TSO250214/gcnvcaller/TSO250214-model \
  --calls-shard-path ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller/${SAMPLE}-calls \
  --sample-index 0 \
  --output-genotyped-intervals ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz \
  --output-genotyped-segments ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  --output-denoised-copy-ratios ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller/${SAMPLE}_denoised_copy_ratios.tsv \
  --contig-ploidy-calls ${WDIR}/output/${COHORT}/${SAMPLE}/ploidy-calls/ \
  --allosomal-contig chrX --allosomal-contig chrY \
  -XL ${WDIR}/refs/par.bed \
  --sequence-dictionary ${WDIR}/refs/hg38.dict

echo "Finished processing CNVs"

gatk VariantFiltration \
  -V ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  -filter "QUAL < 100.0" --filter-name "CNVQUAL" \
  -filter "QUAL < 30.0" --filter-name "CNVRMV" \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz

echo "Finished filtering CNV calls"

zgrep -P -v "CNVRMV|N\t\." ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
  | sed 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}.cnv.vcf.gz
tabix ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}.cnv.vcf.gz

echo "Done with CNVs!"

configManta.py \
  --referenceFasta ${REF} \
  --callRegions refs/ces_manta.bed.gz \
  --bam ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}.bam \
  --runDir ${WDIR}/output/${COHORT}/${SAMPLE}/manta/ \
  --exome
  
${WDIR}/output/${COHORT}/${SAMPLE}/manta/runWorkflow.py -j ${THREADS}

mv ${WDIR}/output/${COHORT}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
  ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}.sv.vcf.gz
mv ${WDIR}/output/${COHORT}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz.tbi \
  ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}.sv.vcf.gz.tbi

echo "Done with SVs!"

rm ${WDIR}/output/${COHORT}/${SAMPLE}/bams/*_* ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/*_*

echo "All Done!"