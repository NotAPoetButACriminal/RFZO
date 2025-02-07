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

SAMPLE=$1
WDIR="/lustre/imgge/RFZO"
REF=$WDIR/refs/hg38.fasta

eval "$(conda shell.bash hook)"
conda activate gatk

gatk CollectReadCounts \
  -R ${REF} \
  -L ${WDIR}/refs/read_counts_wes.interval_list \
  -imr OVERLAPPING_ONLY \
  -I ${WDIR}/output/${SAMPLE}/bams/${SAMPLE}.bam \
  -O ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5

echo "Finished counting reads!"

SCATTERS=$(basename ${WDIR}/output/WES250203/interval_scatters/temp_0001_of_* \
  | cut -d "_" -f 4)

gatk DetermineGermlineContigPloidy \
  --model ${WDIR}/output/WES250203/ploidy-model/ \
  -I ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5 \
  -O ${WDIR}/output/${SAMPLE}/ \
  --output-prefix ploidy

echo "Finished determining ploidy!"

for SCATTER in $(seq -w 0001 00${SCATTERS})
do
  gatk GermlineCNVCaller \
    --run-mode CASE \
    --model ${WDIR}/output/WES250203/gcnvcaller_scatters/scatter_${SCATTER}-model \
    -I ${WDIR}/output/${SAMPLE}/counts/${SAMPLE}.hdf5 \
    -O ${WDIR}/output/${SAMPLE}/gcnvcaller_scatters \
    --output-prefix scatter_${SCATTER} \
    --contig-ploidy-calls ${WDIR}/output/${SAMPLE}/ploidy-calls \
    --verbosity DEBUG
done

echo "Finished calling CNVs per scatter"

MODELS=$(ls -p ${WDIR}/output/WES250203/gcnvcaller_scatters/ \
  | grep model \
  | sed "s#^#--model-shard-path ${WDIR}/output/WES250203/gcnvcaller_scatters/#g")
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
