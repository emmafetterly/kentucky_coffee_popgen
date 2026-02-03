#!/bin/bash

ref="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp/trimmed"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/bwa_mem2_hap1"

t_bwa=20
t_sort=2

mkdir -p "$out"/{bam,stats,logs}

for r1 in "$in"/*_1.trim.fq.gz; do
  r2="${r1%_1.trim.fq.gz}_2.trim.fq.gz"
  sample=$(basename "$r1" _1.trim.fq.gz)

  echo "Mapping: $sample"

  bwa-mem2 mem -t "$t_bwa" "$ref" "$r1" "$r2" 2> "$out/logs/${sample}.bwa.log" \
    | samtools sort -@ "$t_sort" -o "$out/bam/${sample}.hap1.sorted.bam" -

  samtools index "$out/bam/${sample}.hap1.sorted.bam"
  samtools flagstat "$out/bam/${sample}.hap1.sorted.bam" > "$out/stats/${sample}.hap1.flagstat.txt"
done