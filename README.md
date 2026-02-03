# kentucky_coffee_popgen

**Process Overview**
1. Inspect read quality: fastqc and multiqc
2. Trim reads: fastp
3. Map reads to reference genome: bwa-mem2 + samtools
4. Evaluate the alignments: Picard tools and multiQC
5. Call SNPs: GATK


# Inspect read quality

## Step 1: Run fastqc for all samples
Fastqc is a way for us to check sample quality and identify any problematic samples before analysis. 

Guide for interpreation: https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html

### Conda enviornment for qc/trimming

```
#create the conda enviornment
conda create -n qc \
  -c conda-forge -c bioconda \
  fastqc \
  multiqc \
  fastp \
  python=3.9
  
 #activate 
 conda activate qc  
 
 #deactivate
 conda deactivate  
```

### 01_fastqc.sh
```
#!/bin/bash

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/01.RawData/all_fastq"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastqc_results"

mkdir -p "$out"

for f in "$in"/*.fq.gz; do
    fastqc -t 8 -o "$out" "$f"
done
```

Multiqc allows us to collate the fastqc reports into one report to view all samples together.

### 02_multiqc.sh

```
#!/bin/bash

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastqc_results"
out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/multiqc"

mkdir -p "$out"
multiqc "$in" -o "$out"
```

Viewing multiqc reports:
```
#how to download the multi qc file
#run this on your local machine

scp ewf7555@10.2.0.53:/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/multiqc/multiqc_report.html .
```


# Trim Reads

We will try using [fastp](https://github.com/OpenGene/fastp) to trim our adaptors and filter for low quality reads. This program automatically detects our adaptor type.


### 03_fastp.sh

```
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
```

After running 03_fastp.sh we will want to check our output files. We will use multiqc again to collate these reports.

Code for running the multiqc report on the fastp output:
```
multiqc /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp   -o /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp/multiqc
```

# Map to reference genome

We will now map our trimmed reads to the reference genome. 

### Indexing the genome
First, we need to prepare the reference file by indexing the genome.


1. Create a conda environment to use for mapping
```
#create a conda environment for mapping

#mapping conda environment
conda create -n mapping -c bioconda -c conda-forge bwa-mem2 samtools
conda activate mapping
```

Now activate the environment and index the reference genome
```
#creating index for mapping
conda activate mapping
cd /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/reference

#index using bwa-mem2
bwa-mem2 index Gymnocladus_dioicus_M_hap1.fa 

#index using samtools
samtools faidx Gymnocladus_dioicus_M_hap1.fa 
```

## Mapping script

### 04_bwamem.sh

```
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
```
What is this script doing?

This script maps paired-end, adapter-trimmed Illumina reads to the Gymnocladus dioicus male haplotype 1 reference genome using bwa-mem2, then generates sorted and indexed BAM files along with basic mapping statistics for each sample. It does this by:

1. Defines input and output locations

```
ref="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
#ref: the haplotype 1 reference genome FASTA file

in="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp/trimmed"
#in: directory containing trimmed paired-end FASTQ files produced by fastp

out="/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/bwa_mem2_hap1"
out: directory where mapping outputs (BAMs, logs, statistics) will be written
```
2. Sets computational resources

```
t_bwa=20 # Uses 20 threads for read alignment with bwa-mem2
t_sort=2 # Uses 2 threads for BAM sorting with samtools sort
```
3. Creates output directories

```
mkdir -p "$out"/{bam,stats,logs}

#bam/ → sorted BAM files
#stats/ → mapping summaries (samtools flagstat)
#logs/ → bwa-mem2 log files for each sample
```
4. Loops over all trimmed read pairs

```
for r1 in "$in"/*_1.trim.fq.gz; do
  r2="${r1%_1.trim.fq.gz}_2.trim.fq.gz"
  
#Identifies all read 1 files matching *_1.trim.fq.gz
#Automatically infers the corresponding read 2 file (*_2.trim.fq.gz)
```
5. Extracts a sample name from the filename
```
for r1 in "$in"/*_1.trim.fq.gz; do
  r2="${r1%_1.trim.fq.gz}_2.trim.fq.gz"
  sample=$(basename "$r1" _1.trim.fq.gz)
```
7. Maps reads to the reference genome by running bwa-mem2 to align paired-end reads to the haplotype 1 reference
```
bwa-mem2 mem -t "$t_bwa" "$ref" "$r1" "$r2" 2> "$out/logs/${sample}.bwa.log" \
```
10. Alignment output is sent directly into samtools sort to produce a coordinate-sorted BAM file

```
samtools index "$out/bam/${sample}.hap1.sorted.bam"
```

9. Generates basic alignment statistics (mapped reads, properly paired reads) using samtools flagstat

```
  samtools flagstat "$out/bam/${sample}.hap1.sorted.bam" > "$out/stats/${sample}.hap1.flagstat.txt"
```

#### Outputs produced

1. **A sorted BAM file:**
sample.hap1.sorted.bam
2. **A BAM index file**:
sample.hap1.sorted.bam.bai
3. **A mapping summary file:**
sample.hap1.flagstat.txt
4. **A bwa-mem2 log file** documenting the alignment process
