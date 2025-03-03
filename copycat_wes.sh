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

gatk DetermineGermlineContigPloidy \
  -L ${WDIR}/output/${COHORT}/filtered.interval_list \
  -imr OVERLAPPING_ONLY \
  $HDF5S \
  -O ${WDIR}/output/${COHORT}/ \
  --output-prefix ploidy \
  --contig-ploidy-priors ${WDIR}/refs/contig_ploidy_priors.tsv

echo "Finished determining ploidy!"

gatk GermlineCNVCaller \
  --run-mode COHORT \
  -L ${WDIR}/output/${COHORT}/filtered.interval_list \
  --annotated-intervals ${WDIR}/refs/read_counts_wes_annotated.interval_list \
  -imr OVERLAPPING_ONLY \
  $HDF5S \
  -O ${WDIR}/output/${COHORT}/gcnvcaller \
  --output-prefix ${COHORT} \
  --contig-ploidy-calls ${WDIR}/output/${COHORT}/ploidy-calls \
  --verbosity DEBUG

echo "Finished calling CNVs per scatter"

for i in $(seq 0 $((${#SAMPLES[@]} -1)))
do 
  gatk PostprocessGermlineCNVCalls \
    --model-shard-path ${WDIR}/output/${COHORT}/gcnvcaller/${COHORT}-model/ \
    --calls-shard-path ${WDIR}/output/${COHORT}/gcnvcaller/${COHORT}-calls/ \
    --sample-index ${i} \
    --output-genotyped-intervals ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_intervals.cnv.vcf.gz \
    --output-genotyped-segments ${WDIR}/output/${COHORT}/vcfs/${SAMPLES[${i}]}_raw.cnv.vcf.gz \
    --output-denoised-copy-ratios ${WDIR}/output/${COHORT}/gcnvcaller_scatters/${SAMPLES[${i}]}_denoised_copy_ratios.tsv \
    --contig-ploidy-calls ${WDIR}/output/${COHORT}/ploidy-calls/ \
    --allosomal-contig chrX --allosomal-contig chrY \
    --sequence-dictionary ${WDIR}/refs/hg38.dict &
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
  | sed -e 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
  tabix ${WDIR}/output/${COHORT}/vcfs/${SAMPLE}.cnv.vcf.gz
done

echo "Done with CNVs!"

rm ${WDIR}/output/${COHORT}/bams/*_* ${WDIR}/output/${COHORT}/vcfs/*_*

echo "All Done!"