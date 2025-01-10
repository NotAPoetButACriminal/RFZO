#!/bin/bash
#
#SBATCH -J copycat
#SBATCH --output /lustre/imgge/PharmGenHub/logs/%x_%A.out
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
WDIR="/lustre/imgge/PharmGenHub"
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

mkdir -p counts/${COHORT}/

gatk FilterIntervals \
  -L ${WDIR}/refs/read_counts_wes.interval_list \
  --annotated-intervals ${WDIR}/refs/read_counts_wes_annotated.interval_list \
  -imr OVERLAPPING_ONLY \
  $INPUTHDF5S \
  -O ${WDIR}/counts/${COHORT}/${COHORT}_filtered.interval_list \
  --low-count-filter-percentage-of-samples 65

gatk IntervalListTools \
  -I ${WDIR}/counts/${COHORT}/${COHORT}_filtered.interval_list \
  -O ${WDIR}/counts/${COHORT}/${COHORT}_scattered_intervals \
  --SUBDIVISION_MODE INTERVAL_COUNT \
  --SCATTER_CONTENT 13000

gatk DetermineGermlineContigPloidy \
  -L ${WDIR}/counts/${COHORT}/${COHORT}_filtered.interval_list \
  -imr OVERLAPPING_ONLY \
  $INPUTHDF5S \
  -O ${WDIR}/counts/${COHORT} \
  --output-prefix ${COHORT}_ploidy \
  --contig-ploidy-priors ${WDIR}/refs/contig_ploidy_priors.tsv

for SCATTER in {0001..0020}
do
  gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L ${WDIR}/counts/${COHORT}/${COHORT}_scattered_intervals/temp_${SCATTER}_of_20/scattered.interval_list \
    --annotated-intervals ${WDIR}/refs/read_counts_wes_annotated.interval_list \
    -imr OVERLAPPING_ONLY \
    $INPUTHDF5S \
    -O ${WDIR}/counts/${COHORT}/${COHORT}_scatters/ \
    --output-prefix ${COHORT}_scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/counts/${COHORT}/${COHORT}_ploidy-calls \
    --verbosity DEBUG &
done

wait

MODELS=$(ls -p ${WDIR}/counts/${COHORT}/${COHORT}_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/counts/${COHORT}/${COHORT}_scatters/#g")
CALLS=$(ls -p ${WDIR}/counts/${COHORT}/${COHORT}_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/counts/${COHORT}/${COHORT}_scatters/#g")

for i in $(seq 0 $(( $(echo $SAMPLES | wc -w ) - 1 )))
do 
  SAMPLE_PREFIX=$(cat ${WDIR}/counts/${COHORT}/${COHORT}_scatters/${COHORT}_scatter_0001-calls/SAMPLE_${i}/sample_name.txt)
  gatk PostprocessGermlineCNVCalls \
    $MODELS \
    $CALLS \
    --sample-index ${i} \
    --output-genotyped-intervals ${WDIR}/vcfs/${SAMPLE_PREFIX}_intervals.cnv.vcf.gz \
    --output-genotyped-segments ${WDIR}/vcfs/${SAMPLE_PREFIX}_raw.cnv.vcf.gz \
    --output-denoised-copy-ratios ${WDIR}/counts/${COHORT}/${SAMPLE_PREFIX}_denoised_copy_ratios.tsv \
    --contig-ploidy-calls ${WDIR}/counts/${COHORT}/${COHORT}_ploidy-calls/ \
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
