#!/bin/bash
set -euo pipefail

in="/home/ewf7555/Gymnocladus_dioicus/output_files/bwa_mem2_hap1/bam"
out="/home/ewf7555/Gymnocladus_dioicus/output_files/gatk_preproc"
mkdir -p "$out"/{rg,md,metrics,logs}

for bam in "$in"/*.hap1.sorted.bam; do
  base=$(basename "$bam" .hap1.sorted.bam)

  # skip Undetermined 
  if [[ "$base" == Undetermined* ]]; then
    echo "Skipping Undetermined BAM: $base"
    continue
  fi

  # parse fields 
  sm=$(echo "$base" | cut -d_ -f1-3)          # Gymno_2013_109
  lb=$(echo "$base" | cut -d_ -f4)            # CKDL250009712-1A
  pu=$(echo "$base" | awk -F_ '{print $(NF-1)"_"$NF}')  # 22V7CCLT4_L2
  rgid="${sm}_${pu}"                          # Gymno_2013_109_22V7CCLT4_L2

  log="$out/logs/${base}.log"
  echo "[$(date)] RG + MarkDuplicates: $base (SM=$sm, LB=$lb, PU=$pu, RGID=$rgid)" | tee "$log"

  gatk AddOrReplaceReadGroups \
    -I "$bam" \
    -O "$out/rg/${base}.rg.bam" \
    -RGID "$rgid" \
    -RGLB "$lb" \
    -RGPL "ILLUMINA" \
    -RGPU "$pu" \
    -RGSM "$sm" \
    2>>"$log"

  gatk MarkDuplicates \
    -I "$out/rg/${base}.rg.bam" \
    -O "$out/md/${base}.md.bam" \
    -M "$out/metrics/${base}.markdup.metrics.txt" \
    --CREATE_INDEX true \
    2>>"$log"
done
