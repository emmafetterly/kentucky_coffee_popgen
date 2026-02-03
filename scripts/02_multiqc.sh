#!/bin/bash

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastqc_results"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/multiqc"

mkdir -p "$out"
multiqc "$in" -o "$out"