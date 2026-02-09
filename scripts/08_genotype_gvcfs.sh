#!/bin/bash
set -euo pipefail

ref="/home/ewf7555/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
indir="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_joint_hap1"
outdir="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_joint_hap1"

cohort_gvcf="$indir/gymno_hap1.cohort.g.vcf.gz"
outvcf="$outdir/gymno_hap1.raw.vcf.gz"
log="$outdir/logs/genotype_gvcfs.log"

mkdir -p "$outdir/logs"

if [[ ! -s "$cohort_gvcf" ]]; then
  echo "ERROR: Missing cohort gVCF: $cohort_gvcf" | tee "$log"
  exit 1
fi

echo "[$(date)] GenotypeGVCFs starting" | tee "$log"

gatk --java-options "-Xmx32g" GenotypeGVCFs \
  -R "$ref" \
  -V "$cohort_gvcf" \
  -O "$outvcf" \
  2>>"$log"

echo "[$(date)] GenotypeGVCFs done: $outvcf" | tee -a "$log"
