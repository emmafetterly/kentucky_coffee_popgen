#!/bin/bash

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/01.RawData/all_fastq"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastqc_results"

mkdir -p "$out"

for f in "$in"/*.fq.gz; do
    fastqc -t 8 -o "$out" "$f"
done