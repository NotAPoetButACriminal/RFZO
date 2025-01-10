#!/bin/bash
#
#SBATCH -J copycat
#SBATCH --output /lustre/imgge/RFZO/logs/%x_%A.out
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

module load samtools
module load gatk
module load miniconda3

set -ex

SAMPLE=$1
WDIR="/lustre/imgge/RFZO"
REF=$WDIR/refs/hg38.fasta

eval "$(conda shell.bash hook)"
conda activate gatk

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_ces.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/bams/${SAMPLE}.bam \
  -O ${WDIR}/counts/${SAMPLE}.hdf5

mkdir -p counts/${SAMPLE}/

gatk DetermineGermlineContigPloidy \
  --model ${WDIR}/counts/TSO20240108/ploidy/ploidy-model/ \
  -I ${WDIR}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/counts/${SAMPLE} \
  --output-prefix ${SAMPLE}_ploidy

for SCATTER in {0001..0019}
do
  gatk GermlineCNVCaller \
    --run-mode CASE  \
    --model ${WDIR}/counts/TSO20240108/gcnvcaller_scatters/scatter_${SCATTER}-model \
    -I ${WDIR}/counts/${SAMPLE}.hdf5 \
    -O ${WDIR}/counts/${SAMPLE}/${SAMPLE}_scatters/ \
    --output-prefix ${SAMPLE}_scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/counts/${SAMPLE}/${SAMPLE}_ploidy-calls \
    --verbosity DEBUG &
done

wait

MODELS=$(ls -p ${WDIR}/counts/${SAMPLE}/${SAMPLE}_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/counts/${SAMPLE}/${SAMPLE}_scatters/#g")
CALLS=$(ls -p ${WDIR}/counts/${SAMPLE}/${SAMPLE}_scatters/ \
  | grep calls \
  | sed "s#^#--calls-shard-path ${WDIR}/counts/${SAMPLE}/${SAMPLE}_scatters/#g")

gatk PostprocessGermlineCNVCalls \
  --sample-index 0 \
  --allosomal-contig chrX \
  --allosomal-contig chrY \
  --output-genotyped-intervals ${WDIR}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz \
  --output-genotyped-segments ${WDIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  --output-denoised-copy-ratios ${WDIR}/counts/${SAMPLE}/${SAMPLE}_denoised_copy_ratios.tsv \
  --contig-ploidy-calls ${WDIR}/counts/${SAMPLE}/${SAMPLE}_ploidy-calls \
  ${MODELS} \
  ${CALLS}

gatk VariantFiltration \
  -V ${WDIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz \
  -filter "QUAL < 30.0" \
  --filter-name "CNVQUAL" \
  -O ${WDIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz

zcat ${WDIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz \
  | sed -e 's/##source=VariantFiltration/##source=VariantFiltration\n##reference=hg38.fasta/g' \
  -e 's/\tEND/\tSVTYPE=CNV;END/g' \
  | bgzip -o vcfs/${SAMPLE}.cnv.vcf.gz \
  && tabix vcfs/${SAMPLE}.cnv.vcf.gz

echo "All done!"
