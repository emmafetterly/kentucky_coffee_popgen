#!/bin/bash

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/01.RawData/all_fastq"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp"

mkdir -p "$out/trimmed" "$out/html" "$out/json"

for r1 in "$in"/*_1.fq.gz; do
  r2="${r1%_1.fq.gz}_2.fq.gz"

  if [ ! -f "$r2" ]; then
    echo "Skipping: missing R2 for $r1"
    continue
  fi

  base=$(basename "$r1" _1.fq.gz)

  fastp \
    -i "$r1" -I "$r2" \
    -o "$out/trimmed/${base}_1.trim.fq.gz" \
    -O "$out/trimmed/${base}_2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread 6 \
    --html "$out/html/${base}.fastp.html" \
    --json "$out/json/${base}.fastp.json"
done