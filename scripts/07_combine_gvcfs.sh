#!/bin/bash
set -euo pipefail

ref="/home/ewf7555/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
gdir="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_gvcf_hap1/gvcf"
outdir="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_joint_hap1"

mkdir -p "$outdir"/{logs,tmp}

cohort_gvcf="$outdir/gymno_hap1.cohort.g.vcf.gz"
log="$outdir/logs/combine_gvcfs.log"

echo "[$(date)] CombineGVCFs starting" | tee "$log"

# Build a temp list of -V arguments (one per gVCF)
vlist="$outdir/tmp/gvcf_inputs.list"
ls "$gdir"/*.g.vcf.gz > "$vlist"

# Convert list file into repeated "-V file" args
# shellcheck disable=SC2046
gatk --java-options "-Xmx32g" CombineGVCFs \
  -R "$ref" \
  $(awk '{print "-V", $0}' "$vlist") \
  -O "$cohort_gvcf" \
  2>>"$log"

echo "[$(date)] CombineGVCFs done: $cohort_gvcf" | tee -a "$log"
