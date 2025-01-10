#!/bin/bash
#
#SBATCH -J copycat_ces
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

module load samtools
module load gatk
module load miniconda3

set -ex

SAMPLES=$(cat $1)
COHORT=$2
WDIR="/lustre/imgge/RFZO"
REF=$WDIR/refs/hg38.fasta
DICT=$WDIR/refs/hg38.dict

eval "$(conda shell.bash hook)"
conda activate gatk

INPUTHDF5S=$(
  for i in $SAMPLES
  do
    echo -n "-I ${WDIR}/counts/${i}.hdf5 "
  done
)

mkdir -p counts/${COHORT}/ploidy

gatk FilterIntervals \
  -L ${WDIR}/refs/read_counts_ces.interval_list \
  --annotated-intervals ${WDIR}/refs/read_counts_ces_annotated.interval_list \
  -imr OVERLAPPING_ONLY \
  $INPUTHDF5S \
  -O ${WDIR}/counts/${COHORT}/filtered.interval_list \
  --low-count-filter-percentage-of-samples 65

gatk IntervalListTools \
  -I ${WDIR}/counts/${COHORT}/filtered.interval_list \
  -O ${WDIR}/counts/${COHORT}/interval_scatters \
  --SUBDIVISION_MODE INTERVAL_COUNT \
  --SCATTER_CONTENT 3000

SCATTERS=$(basename ${WDIR}/counts/${COHORT}/interval_scatters/temp_0001_of_* \
  | cut -d "_" -f 4)

gatk DetermineGermlineContigPloidy \
  -L ${WDIR}/counts/${COHORT}/filtered.interval_list \
  -imr OVERLAPPING_ONLY \
  $INPUTHDF5S \
  -O ${WDIR}/counts/${COHORT}/ploidy \
  --output-prefix ploidy \
  --contig-ploidy-priors ${WDIR}/refs/contig_ploidy_priors.tsv

for SCATTER in $(seq -w 0001 00${SCATTERS})
do
  gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L ${WDIR}/counts/${COHORT}/interval_scatters/temp_${SCATTER}_of_${SCATTERS}/scattered.interval_list \
    --annotated-intervals ${WDIR}/refs/read_counts_ces_annotated.interval_list \
    -imr OVERLAPPING_ONLY \
    $INPUTHDF5S \
    -O ${WDIR}/counts/${COHORT}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/counts/${COHORT}/ploidy/ploidy-calls \
    --verbosity DEBUG &
done

wait

MODELS=$(ls -p ${WDIR}/counts/${COHORT}/gcnvcaller_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/counts/${COHORT}/gcnvcaller_scatters/#g")
CALLS=$(ls -p ${WDIR}/counts/${COHORT}/gcnvcaller_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/counts/${COHORT}/gcnvcaller_scatters/#g")

for i in $(seq 0 $(( $(echo $SAMPLES | wc -w ) - 1 )))
do 
  SAMPLE_PREFIX=$(cat ${WDIR}/counts/${COHORT}/gcnvcaller_scatters/scatter_0001-calls/SAMPLE_${i}/sample_name.txt)
  gatk PostprocessGermlineCNVCalls \
    $MODELS \
    $CALLS \
    --sample-index ${i} \
    --output-genotyped-intervals ${WDIR}/vcfs/${SAMPLE_PREFIX}_intervals.cnv.vcf.gz \
    --output-genotyped-segments ${WDIR}/vcfs/${SAMPLE_PREFIX}_raw.cnv.vcf.gz \
    --output-denoised-copy-ratios ${WDIR}/counts/${COHORT}/${SAMPLE_PREFIX}_denoised_copy_ratios.tsv \
    --contig-ploidy-calls ${WDIR}/counts/${COHORT}/ploidy/ploidy-calls/ \
    --allosomal-contig chrX --allosomal-contig chrY \
    --sequence-dictionary $DICT &
done

wait

for SAMPLE in $SAMPLES
do
  gatk VariantFiltration \
    -V vcfs/${SAMPLE}_raw.cnv.vcf.gz \
    -filter "QUAL < 30.0" \
    --filter-name "CNVQUAL" \
    -O vcfs/${SAMPLE}_filtered.cnv.vcf.gz &
done

wait

for SAMPLE in $SAMPLES
do
zgrep -P -v "CNVQUAL|N\t\." vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
  | sed -e 's/##source=VariantFiltration/##source=VariantFiltration\n##reference=hg38.fasta/g' \
        -e 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o vcfs/${SAMPLE}.cnv.vcf.gz \
  && tabix vcfs/${SAMPLE}.cnv.vcf.gz
done

echo "All done!"
