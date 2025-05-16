#!/bin/bash
#
#SBATCH -J copycat_wes
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

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
INTERVALS="${WDIR}/refs/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_Combined_Mito.bed"

mkdir -p \
  output/${COHORT}/bams/metrics \
  output/${COHORT}/vcfs \
  output/${COHORT}/counts

### START ###

eval "$(conda shell.bash hook)"
conda activate gatk

for SAMPLE in "${SAMPLES[@]}"
do
  gatk CollectReadCounts \
    -R ${REF} \
    -L ${WDIR}/refs/read_counts_wes.interval_list \
    -imr OVERLAPPING_ONLY \
    -I ${WDIR}/input/${SAMPLE}.bam \
    -O ${WDIR}/output/${COHORT}/counts/${SAMPLE}.hdf5 &
done

wait

echo "Finished counting reads!"

HDF5S=$(ls ${WDIR}/output/${COHORT}/counts/*.hdf5 | sed -e 's/^/ -I /g')

gatk FilterIntervals \
  -L ${WDIR}/refs/read_counts_wes.interval_list \
  --annotated-intervals ${WDIR}/refs/read_counts_wes_annotated.interval_list \
  -imr OVERLAPPING_ONLY \
  $HDF5S \
  -O ${WDIR}/output/${COHORT}/filtered.interval_list \
  --low-count-filter-percentage-of-samples 65

gatk IntervalListTools \
  -I ${WDIR}/output/${COHORT}/filtered.interval_list  \
  -O ${WDIR}/output/${COHORT}/interval_scatters \
  --SUBDIVISION_MODE INTERVAL_COUNT \
  --SCATTER_CONTENT 26000

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
    --annotated-intervals ${WDIR}/refs/read_counts_wes_annotated.interval_list \
    -imr OVERLAPPING_ONLY \
    $HDF5S \
    -O ${WDIR}/output/${COHORT}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/output/${COHORT}/ploidy-calls \
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
  (SAMPLE=$(cat ${WDIR}/output/${COHORT}/gcnvcaller_scatters/scatter_0001-calls/SAMPLE_${i}/sample_name.txt)
  gatk PostprocessGermlineCNVCalls \
    $MODELS \
    $CALLS \
    --sample-index ${i} \
    --output-genotyped-intervals ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz \
    --output-genotyped-segments ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
    --output-denoised-copy-ratios ${WDIR}/output/${COHORT}/gcnvcaller_scatters/${SAMPLE}_denoised_copy_ratios.tsv \
    --contig-ploidy-calls ${WDIR}/output/${COHORT}/ploidy-calls/ \
    --allosomal-contig chrX --allosomal-contig chrY \
    --sequence-dictionary ${WDIR}/refs/hg38.dict) &
done
wait

echo "Finished calling CNVs per sample"

for SAMPLE in "${SAMPLES[@]}"
do
  gatk VariantFiltration \
    -V ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
    -filter "QUAL < 30.0" \
    --filter-name "CNVQUAL" \
    -O ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz &
done
wait

echo "Finished filtering CNV calls"

for SAMPLE in "${SAMPLES[@]}"
do
  zgrep -P -v "CNVQUAL|N\t\." ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
  | sed -e 's/##source=VariantFiltration/##source=VariantFiltration\n##reference=hg38.fasta/g' \
    -e 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
  tabix ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
done

echo "Done with CNVs!"

for SAMPLE in "${SAMPLES[@]}"
do
  (configManta.py \
    --referenceFasta ${REF} \
    --callRegions refs/wes_manta.bed.gz \
    --bam ${WDIR}/input/${SAMPLE}.bam \
    --runDir ${WDIR}/output/${COHORT}/manta/${SAMPLE}/ \
    --exome
    ${WDIR}/output/${COHORT}/manta/${SAMPLE}/runWorkflow.py -j 2
    mv ${WDIR}/output/${COHORT}/manta/${SAMPLE}/results/variants/diploidSV.vcf.gz \
      ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.sv.vcf.gz
    mv ${WDIR}/output/${COHORT}/manta/${SAMPLE}/results/variants/diploidSV.vcf.gz.tbi \
      ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.sv.vcf.gz.tbi) &
done
wait

echo "Done with SVs!"

rm ${WDIR}/output/${COHORT}/vcfs/*_*

echo "All Done!"