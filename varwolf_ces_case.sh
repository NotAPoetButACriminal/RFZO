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
module load miniconda3

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
  FLOWCELL=$(zcat "${WDIR}"/input/${SAMPLE}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 3)
  LIBRARY=$(zcat ${WDIR}/input/${SAMPLE}*L001_R1*.fastq.gz | head -1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
  fastp \
    -w 4 \
    -i ${WDIR}/input/${SAMPLE}*L00${LANE}_R1*.fastq.gz \
    -I ${WDIR}/input/${SAMPLE}*L00${LANE}_R2*.fastq.gz \
    --stdout \
    -j ${WDIR}/output/${SAMPLE}/fastp.json \
    -h ${WDIR}/output/${SAMPLE}/fastp.html \
  | bwa mem \
	  -t ${THREADS} \
	  -M -p \
    -R "@RG\tID:${FLOWCELL}.LANE${LANE}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE}" \
    ${REF} - \
  | samtools sort \
	  -@ 4 \
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

gatk HaplotypeCaller \
  -R ${REF} \
  -L ${INTERVALS} \
  -ip 10 \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz
echo "Finished calling variants!"

gatk VariantFiltration \
  -V ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.vcf.gz \
  -filter "DP < 5.0" --filter-name "DP5" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_filtered.vcf.gz

echo "Finished filtering variants!"


zcat ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_filtered.vcf.gz \
  | sed 's/##source=HaplotypeCaller/##source=HaplotypeCaller\n##reference=hg38.fasta/g' \
  | bgzip -@ 2 -o ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}.vcf.gz
tabix ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}.vcf.gz

echo "Done with SNVs!"

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_ces.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5

echo "Finished counting reads!"

eval "$(conda shell.bash hook)"
conda activate gatk

SCATTERS=$(basename ${WDIR}/output/TSO250117/interval_scatters/temp_0001_of_* \
  | cut -d "_" -f 4)

gatk DetermineGermlineContigPloidy \
  --model ${WDIR}/output/TSO250117/ploidy-model/ \
  -I ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/output/${SAMPLE}/ \
  --output-prefix ploidy

echo "Finished determining ploidy!"

for SCATTER in $(seq -w 0001 00${SCATTERS})
do
  gatk GermlineCNVCaller \
    --run-mode CASE \
    --model ${WDIR}/output/TSO250117/gcnvcaller_scatters/scatter_${SCATTER}-model \
    -I ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5 \
    -O ${WDIR}/output/${SAMPLE}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/output/${SAMPLE}/ploidy-calls \
    --verbosity DEBUG
done

echo "Finished calling CNVs per scatter"

MODELS=$(ls -p ${WDIR}/output/TSO250117/gcnvcaller_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/output/TSO250117/gcnvcaller_scatters/#g")
CALLS=$(ls -p ${WDIR}/output/${SAMPLE}/gcnvcaller_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/output/${SAMPLE}/gcnvcaller_scatters/#g")


gatk PostprocessGermlineCNVCalls \
  $MODELS \
  $CALLS \
  --sample-index 0 \
  --output-genotyped-intervals ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz \
  --output-genotyped-segments ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  --output-denoised-copy-ratios ${WDIR}/output/${SAMPLE}/gcnvcaller_scatters/${SAMPLE}_denoised_copy_ratios.tsv \
  --contig-ploidy-calls ${WDIR}/output/${SAMPLE}/ploidy-calls/ \
  --allosomal-contig chrX --allosomal-contig chrY \
  --sequence-dictionary ${WDIR}/refs/hg38.dict

echo "Finished calling CNVs per sample"

gatk VariantFiltration \
  -V ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  -filter "QUAL < 30.0" \
  --filter-name "CNVQUAL" \
  -O ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz

echo "Finished filtering CNV calls"

zgrep -P -v "CNVQUAL|N\t\." ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
  | sed -e 's/##source=VariantFiltration/##source=VariantFiltration\n##reference=hg38.fasta/g' \
    -e 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}.cnv.vcf.gz
tabix ${WDIR}/output/${SAMPLE}/vcfs/${SAMPLE}.cnv.vcf.gz

echo "Done with CNVs!"

rm ${WDIR}/output/${SAMPLE}/bams/*_* ${WDIR}/output/${SAMPLE}/vcfs/*_*

echo "All Done!"