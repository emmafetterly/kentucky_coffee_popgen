#!/bin/bash
set -euo pipefail

ref="/home/ewf7555/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
in="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_preproc/md"
out="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_gvcf_hap1"

threads=8   # set number of threads
mkdir -p "$out"/{gvcf,logs}

for bam in "$in"/*.md.bam; do
  base=$(basename "$bam" .md.bam)

  gvcf="$out/gvcf/${base}.g.vcf.gz"
  log="$out/logs/${base}.haplotypecaller.log"

  echo "[$(date)] HaplotypeCaller: $base" | tee "$log"

  gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R "$ref" \
    -I "$bam" \
    -O "$gvcf" \
    -ERC GVCF \
    --native-pair-hmm-threads "$threads" \
    2>>"$log"
done
