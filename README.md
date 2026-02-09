# kentucky_coffee_popgen

**Process Overview**
1. Inspect read quality: fastqc and multiqc
2. Trim reads: fastp
3. Map reads to reference genome: bwa-mem2 + samtools
4. Evaluate the alignments: Picard tools and multiQC
5. Call SNPs: GATK
6. Filter SNPs: vcftools
7. Pop gen analysis

![Screenshot 2026-02-03 at 12.53.01 PM](https://hackmd.io/_uploads/BJIQVTyDZe.png)

----


## File architecture

**Main directory in Curie:** `/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/`

![Screenshot 2026-02-03 at 3.58.46 PM](https://hackmd.io/_uploads/rkji1elPZl.png)

**This directory has three major divisions**
1. input_files: raw data 
2. output_files: outputs of scripts
3. scripts: script files and .txt note files

#### 1. Input files:
File path: `/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files`

![Screenshot 2026-01-30 at 3.51.30 PM](https://hackmd.io/_uploads/BJegdjc8Wl.png)


#### Input files contains: 


1. Gymnocladus_dioicus_male_assembly 
```
This directory contains the assemblies and associated data for 'Phased Chromosome-Scale Assembly of Male Gymnocladus dioicus'.

	- Gymnocladus_dioicus_M_hap1.fa - Haplotype 1 fasta assembly
	- Gymnocladus_dioicus_M_hap2.fa - Haplotype 2 fasta assembly
	- contam_hap1.fasta - Fragments removed from haplotype 1 during FCS-GX contaminant screening
	- contam_hap2.fasta - Fragments removed from haplotype 2 during FCS-GX contaminant screening
	- hap1_gymnocladusM_racon_25kb_fragments.fa - Fragments shorter than 25,000 bp removed from haplotype 1 via racon
	- hap2_gymnocladusM_racon_25kb_fragments.fa - Fragments shorter than 25,000 bp removed from hap2 via racon
	- mito_assembly - directory containing mitochondrial assembly and annotation
	- plastid_assembly - directory containing plastid assembly and annotation
    
```

2. NovaSeq raw reads: Each read was binned under a separate folder for each sample, I reorganized the reads to a new folder called "all_fastq"

Code used to move samples over:
```
#moving all fastq files to one folder in the raw data folder
# previously they were nested by sample name

find /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/01.RawData \
  -type f -name "*.fq.gz" \
  -exec cp -v {} /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/all_fastq \;
  
  
#### what is this code doing?
# search starting in this directory
find /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/ \

  # only return regular files (not directories)
  -type f \

  # only files ending in .fq.gz
  -name "*.fq.gz" \

  # for each file found, run the following command (-v verbose - print a message for each file copied)
  -exec cp -v {} \
  
  # copy the file to this destination directory
  /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/NovaSeq/all_fastq \

  # end the -exec command
  \;
```

Code to inspect the zipped files
```
zcat Gymno_2013_109_CKDL250009712-1A_22V7CCLT4_L2_1.fq.gz | head
```

3. reference: this is the directory that we will store our reference genome file and associated index files/sequence dictionary for our pipelines.

To move the hap1 assembly to this directory I used this code:
```
cp /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/Gymnocladus_dioicus_male_assembly/Gymnocladus_dioicus_M_hap1.fa \
   /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/input_files/reference/
```

#### 2. Output files
This will contain our outputs from scripts.

*File path:*
```
/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files
```

#### 3. Scripts
This folder contains all scripts and notes.

*File path:*
```
/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/scripts
```

**Moving files to home directory**

From Suzy - always use a local copy of files in home directory (not running from JBOD).

If you wanted to copy the directory over to your home directory:
`cp -a /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus ~/ewf7555`

or use rsync + screen to avoid interruptions in data transfer

Make sure to replace your own username when running this!
```
screen -L S rsynch_JBOD
rsync -avh /data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus ~/ewf7555
```


# Inspect read quality

## Step 1: Run fastqc for all samples
[Fastqc](https://github.com/s-andrews/FastQC) is a way for us to check sample quality and identify any problematic samples before analysis. 

Guides for interpreation: 
[Helpful guide](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html) 
[Another guide](https://rtsf.natsci.msu.edu/genomics/technical-documents/fastqc-tutorial-and-faq.aspx)

### Conda enviornment for qc/trimming
*Check coding quick start guide for more details on conda*

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
 
 # to deactivate (after you are done)
 conda deactivate  
```

*Process for running scripts:*

1. Write your script, save and then type `chmod +x script_name.sh` in your terminal - you must be in the directory where the script is located to do this and you only need to do this once per script.
3. In your terminal type`screen -L -S fastqc` This will open a new screen session named "fastqc" with a log file.
4. In the new screen type `conda activate qc` to activate the "qc" conda environment. You should see the (base) change to (qc) before your username.
5. Then in the screen type `bash script_name.sh` + enter
6. Then type Ctrl + A   then   D to detach the screen session. 

To check on your progress, check the log and output files. You can also type `top` or `htop` to check your memory/CPU usage as you script runs. 

#### Notes on running scripts + screens:

Screens let us run longer jobs by “detaching” them from the terminal, allowing processes to continue running after you log out.

Remember to exit (kill) your screen session when the job is finished.

```
#make a script executable/change file permissions so you can run it
chmod +x SCRIPT.sh

#Running a script in screen
screen -L bash SCRIPT.sh

#close a screen
crtl + A D or screen -d

#check how many screens you have open 
screen -ls

#iniate and name a screen 
screen -S screen_name

#how to kill a screen session
ctrl + A or type exit or screen -S screen_name -X quit
```


### 01_fastqc.sh

```
#!/bin/bash

BASE="/home/ewf7555/Gymnocladus_dioicus"
IN="$BASE/input_files/NovaSeq/01.RawData/all_fastq"
OUT="$BASE/output_files/fastqc_results"

THREADS=20

mkdir -p "$OUT"

for f in "$in"/*.fq.gz; do
    fastqc -t "$THREADS" -o "$out" "$f"
done
```

Multiqc allows us to collate the fastqc reports into one report to view all samples together.
*You don't really need a script for this, but it's nice to define the inputs/outputs.*

### 02_multiqc.sh

```
#!/bin/bash

BASE="/home/ewf7555/Gymnocladus_dioicus"
IN="$BASE/output_files/fastqc_results"
OUT="$BASE/output_files/multiqc"

mkdir -p "$OUT"

multiqc "$IN" -o "$OUT" --force
```

Viewing multiqc reports:
```
#how to download the multi qc file
#run this on your local machine

scp ewf7555@10.2.0.53:/home/ewf7555/Gymnocladus_dioicus/output_files/multiqc/multiqc_report.html /Users/emmafetterly/Documents/coffeetree/multiqc_pretrim.html
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
Download to view on your computer (change your username):

```
scp -r ewf7555@10.2.0.53:/data/labs/Love/Gymnocladus_AlphaHudson_data/Gymnocladus_dioicus/output_files/fastp/multiqc \
/Users/emmafetterly/Documents/coffeetree/multiqc_posttrim
```

# Map to reference genome

We will now map our trimmed reads to the reference genome. 

### Indexing the genome
First, we need to prepare the reference file by indexing the genome.


1. Create a conda environment to use for mapping
```
#create a conda environment for mapping
conda create -n mapping -c bioconda -c conda-forge bwa-mem2 samtools
conda activate mapping
```

Now activate the environment and index the reference genome
```
#creating index for mapping
conda activate mapping
cd /home/ewf7555/Gymnocladus_dioicus/input_files/reference

#index using bwa-mem2
bwa-mem2 index Gymnocladus_dioicus_M_hap1.fa 

#index using samtools
samtools faidx Gymnocladus_dioicus_M_hap1.fa 
```

## Mapping script

We will use bwa-mem2 and samtools to map our raw reads to the reference sequence. 

Software:
[bwa-mem2 Github](https://github.com/bwa-mem2/bwa-mem2)
[samtools Github](https://github.com/samtools)

The general process involves mapping samples to the reference genome which creates a .sam file, then the script will convert these sam files to .bam files and sort them by their coordinates to create a sorted .bam file (input for gatk).


### 04_bwamem.sh

```
#!/bin/bash

ref="/home/ewf7555/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
in="/home/ewf7555/Gymnocladus_dioicus/output_files/fastp/trimmed"
out="/home/ewf7555/Gymnocladus_dioicus/output_files/bwa_mem2_hap1"

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
ref="/home/ewf7555/Gymnocladus_dioicus/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"
#ref: the haplotype 1 reference genome FASTA file

in="/home/ewf7555/Gymnocladus_dioicus/output_files/fastp/trimmed"
#in: directory containing trimmed paired-end FASTQ files produced by fastp

out="/home/ewf7555/Gymnocladus_dioicus/output_files/bwa_mem2_hap1"
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


# Calling SNPs

## Step 1: Adding read groups

**What is a read group?**

A read group is a set of reads that share the same origin, (such as the same library preparation, sequencing lane, and flowcell), and that belong to the same biological sample. Read groups are stored in the .bam file header and are referenced by each read to record how and where that read was generated. They help the program to:

* Identify which reads come from the same biological sample (SM)
* Distinguish technical replicates, such as different lanes or libraries (ID, PU, LB)
* Correctly mark PCR duplicates, including duplicates that occur across sequencing lanes
* Model and account for technical variation during variant calling


Here is how our samples are named:
| Field | Value              | Meaning                                       |
| ----- | ------------------ | --------------------------------------------- |
| 1     | `Gymno`            | Species                   |
| 2     | `2013`             | Collection year                               |
| 3     | `109`              | Individual ID                                 |
| 4     | `CKDL250009712-1A` | Library / sample barcode (Illumina sample ID) |
| 5     | `22V7CCLT4`        | Flowcell ID                                   |
| 6     | `L2`               | Lane number                                   |
| 7     | `1`                | Read number (R1)                              |
| 8     | `fq.gz`            | FASTQ (compressed)                            |

The start of our script will define these groups:
```
for bam in "$in"/*.hap1.sorted.bam; do
  base=$(basename "$bam" .hap1.sorted.bam)
  
  # parse fields 
  sm=$(echo "$base" | cut -d_ -f1-3)          # Gymno_2013_109
  lb=$(echo "$base" | cut -d_ -f4)            # CKDL250009712-1A
  pu=$(echo "$base" | awk -F_ '{print $(NF-1)"_"$NF}')  # 22V7CCLT4_L2
  rgid="${sm}_${pu}"    
```

*This can also be done during mapping with bwa-mem2, but since we have already mapped our reads we will do this with gatk.*
```
  gatk AddOrReplaceReadGroups \
    -I "$bam" \
    -O "$out/rg/${base}.rg.bam" \
    -RGID "$rgid" \
    -RGLB "$lb" \
    -RGPL "ILLUMINA" \
    -RGPU "$pu" \
    -RGSM "$sm" \
    2>>"$log"
```

## Step 2: Marking duplicates

Marking duplicates helps us avoid errors to due PCR duplication (when the same piece of DNA gets copied multiple times). This over representation of sequences can result in biased allele frequecies, false confidence in variant ID and sequecing errors can propagate through these duplicates. 

Code for marking duplicates:

```
  gatk MarkDuplicates \
    -I "$out/rg/${base}.rg.bam" \
    -O "$out/md/${base}.md.bam" \
    -M "$out/metrics/${base}.markdup.metrics.txt" \
    --CREATE_INDEX true \
    2>>"$log"
```

### 05_gatk_mark_dups.sh

Putting it all together in one script.

:::spoiler 05_mark_dups.sh
```
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

```
:::

## Step 3: HaplotypeCaller (gVCF)

#### 06_haplotypecaller_gvcf.sh

```
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
```

### 07_combine_gvcfs.sh

```
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
```

### 08_genotype_gvcfs.sh


```
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
```
