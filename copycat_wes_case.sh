#!/bin/bash
#
#SBATCH -J copycat_wes_case
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

module load samtools
module load gatk
module load miniconda3

set -ex

# RUN VARIABLES
SAMPLE=$1
COHORT=$2
WDIR="/lustre/imgge/RFZO"
THREADS=${SLURM_CPUS_PER_TASK}

# REFERENCE VARIABLES
REF="${WDIR}/refs/hg38.fasta"
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
INTERVALS="${WDIR}/refs/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_Combined_Mito.bed"

mkdir -p \
  output/${COHORT}/${SAMPLE}/bams/metrics \
  output/${COHORT}/${SAMPLE}/vcfs \
  output/${COHORT}/${SAMPLE}/counts

eval "$(conda shell.bash hook)"
conda activate gatk

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_wes.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5

echo "Finished counting reads!"

SCATTERS=$(basename ${WDIR}/output/WES250529/interval_scatters/temp_0001_of_* \
  | cut -d "_" -f 4)

gatk DetermineGermlineContigPloidy \
  --model ${WDIR}/output/WES250529/ploidy-model/ \
  -I ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/output/${COHORT}/${SAMPLE}/ \
  --output-prefix ploidy

echo "Finished determining ploidy!"

for SCATTER in $(seq -w 0001 00${SCATTERS})
do
  gatk GermlineCNVCaller \
    --run-mode CASE \
    --model ${WDIR}/output/WES250529/gcnvcaller_scatters/scatter_${SCATTER}-model \
    -I ${WDIR}/output/${COHORT}/${SAMPLE}/counts/${SAMPLE}.hdf5 \
    -O ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/output/${COHORT}/${SAMPLE}/ploidy-calls \
    --verbosity DEBUG &
done
wait

echo "Finished calling CNVs per scatter"

MODELS=$(ls -p ${WDIR}/output/WES250529/gcnvcaller_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/output/WES250529/gcnvcaller_scatters/#g")
CALLS=$(ls -p ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller_scatters/#g")

gatk PostprocessGermlineCNVCalls \
  $MODELS \
  $CALLS \
  --sample-index 0 \
  --output-genotyped-intervals ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz \
  --output-genotyped-segments ${WDIR}/output/${COHORT}/${SAMPLE}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  --output-denoised-copy-ratios ${WDIR}/output/${COHORT}/${SAMPLE}/gcnvcaller_scatters/${SAMPLE}_denoised_copy_ratios.tsv \
  --contig-ploidy-calls ${WDIR}/output/${COHORT}/${SAMPLE}/ploidy-calls/ \
  --allosomal-contig chrX --allosomal-contig chrY \
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
  --callRegions refs/wes_manta.bed.gz \
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
