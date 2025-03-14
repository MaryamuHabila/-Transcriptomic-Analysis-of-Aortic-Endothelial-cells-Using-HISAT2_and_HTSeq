# -Transcriptomic-Analysis-of-Aortic-Endothelial-cells-Using-HISAT2_and_HTSeq
 This guide provides a step-by-step explanation of an RNA-Seq analysis pipeline for beginners using HISAT2 for read alignment and HTSeq for gene quantification. It covers essential preprocessing steps, including quality control, adapter trimming, alignment, and read counting, explaining the purpose and output of each step.

In addition to analysing the samples using STAR_Salmon, HISAT_HTSeq was used. 

**Title: Step-by-Step Explanation of RNA-Seq Analysis Using HISAT2 and HTSeq**

## Introduction
RNA sequencing (RNA-seq) is a powerful technique used to analyze gene expression by sequencing RNA molecules in a sample. In this guide, we will walk through an RNA-seq analysis workflow using HISAT2 for read alignment and HTSeq for gene quantification. This guide is designed for beginners who are new to RNA-seq data analysis.

---

## **Step 1: Preparing Input Data**

Before aligning RNA-seq data, we need to organize the input files. In this analysis, raw sequencing data is stored in FASTQ format and is paired-end (i.e., for each sample, there are two corresponding files: one for forward reads and one for reverse reads). The paths to these files are recorded in a **samplesheet (Samplesheet_HAoECs.csv)**, which helps automate the analysis process.

### Example of the samplesheet:
```
Sample, Read 1, Read 2, Strandedness
C1_311,/path/to/C1_311_1.fq.gz,/path/to/C1_311_2.fq.gz,reverse
C2_212,/path/to/C2_212_1.fq.gz,/home/s2451842@bifx-core2/...
```
Each row in this file represents an RNA-seq sample, with columns for:
- **Sample name** (e.g., `C1_311`)
- **Path to read 1** (forward reads file)
- **Path to paired read file** (reverse reads file)
- **Strandedness** (important for correct gene expression quantification, here it is `reverse` strand-specific)

## Data Processing Workflow

### Step 1: Setting Up Directories and Preparing Adapter Sequences
Before processing the raw sequencing data, we need to organize output directories. We use the command:
```bash
mkdir -p "$OUTPUT_DIR" trimmed_data qc_reports adapters logs multiqc_reports
```
This ensures that outputs like trimmed reads, quality control (QC) reports, and logs are neatly organized.

**Trimmomatic Adapter Preparation:**
The script also checks for the presence of adapter sequence files used in **Trimmomatic**, which is a tool for trimming adapter sequences from raw reads to improve alignment accuracy. If the adapter sequences are not found in the expected directory (`adapters/TruSeq3-PE.fa`), they will be downloaded before proceeding.

## Step 1: Quality Control (FastQC)
**Why?** Raw sequencing data often contains low-quality reads, adapter sequences, and other artifacts. **FastQC** is a tool used to assess the quality of sequencing reads by evaluating factors such as:
- Read length distribution
- GC content (guanine-cytosine proportion)
- Presence of overrepresented sequences (which could indicate contamination or adapters)

Each sample is processed as follows:
```bash
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
    fastqc -t "$THREADS" "$fastq_1" "$fastq_2"
done
```
The `fastqc` tool generates an HTML report for each sequencing file, which can be later examined to assess the quality of the raw reads.

## Step 2: Read Trimming with Trimmomatic
**Why?**
Next, we remove low-quality reads and adapter sequences using **Trimmomatic**, which helps improve the accuracy of downstream analyses. This step filters out poor-quality sequences and trims adapters from the reads.

```bash
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
    trimmomatic PE -threads "$THREADS" \
        "$fastq_1" "$fastq_2" \
        trimmed_data/${sample}_1P.fastq.gz trimmed_data/${sample}_1U.fastq.gz \
        trimmed_data/${sample}_2P.fastq.gz trimmed_data/${sample}_2U.fastq.gz \
        ILLUMINACLIP:"$ADAPTERS":2:30:10
```

## Step 3: Read Alignment Using HISAT2
**Why?**
After trimming, we align the reads to a reference genome using **HISAT2**, a fast and efficient splice-aware aligner designed for RNA-seq data. This step helps map the short sequencing reads back to their genomic locations.

```bash
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
    hisat2 -x "$INDEX_PREFIX" \
        -1 trimmed_data/${sample}_1P.fastq.gz \
        -2 trimmed_data/${sample}_2P.fastq.gz \
        --rna-strandness RF \
        --summary-file "$OUTPUT_DIR/${sample}.hisat2.summary.log" \
        --threads "$THREADS" \
        | samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/${sample}_sorted.bam"
done
```
**What it produces:**
- **Sorted BAM files** for each sample.
- **Summary logs** containing alignment statistics.

## Step 4: Read Quantification Using HTSeq
**Why?**
HTSeq is used to count how many reads map to each gene in the genome. This output will be used in downstream differential expression analysis.

```bash
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness
do
    htseq-count -f bam -s reverse -i gene_id -r pos \n
   "$OUTPUT_DIR/${sample}.bam" "$GTF" > "$OUTPUT_DIR/${sample}_counts.txt"
done


