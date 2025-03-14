# Comprehensive Guide to RNA-Seq Data Processing with HISAT2 and HTSeq

## Introduction
RNA sequencing (RNA-Seq) is a powerful technique used in bioinformatics to analyze gene expression and understand the transcriptome. This guide provides an in-depth explanation of an RNA-Seq data processing pipeline, specifically using HISAT2 for alignment and HTSeq for read quantification. It is intended for bioinformatics students who are new to RNA sequencing data analysis.

## Overview of the Workflow
The RNA-Seq pipeline for this study includes the following major steps:
1. **Raw Data Quality Control** – Assessing the quality of the raw sequencing reads using FastQC.
2. **Preprocessing (Trimming Adapters and Low-Quality Bases)** – Trimming adapter sequences and filtering low-quality reads using Trimmomatic.
3. **Genome Indexing** – Preparing the reference genome index using HISAT2.
4. **Read Alignment** – Aligning the quality-filtered reads to the reference genome using HISAT2.
5. **Post-alignment Processing** – Converting, sorting, and indexing aligned reads.
6. **Gene Quantification** – Counting mapped reads per gene using HTSeq.
7. **Quality Assessment** – Aggregating quality reports using MultiQC.

## Input Data and Samplesheet
The input to the pipeline is a CSV file named `Samplesheet_HAoECs.csv`, which lists RNA-seq samples along with their corresponding paired-end read files and strand specificity:

```
sample,fastq_1,fastq_2,strandedness
C1_311,/home/s2451842/projects/HAoECs_RNA_seq/data/C1_311_1.fq.gz,/home/s2451842/projects/HAoECs_RNA_seq/data/C1_311_2.fq.gz,reverse
...
```

The dataset includes RNA-Seq samples from human aortic endothelial cells (HAoECs) that are either control (C) or treated (LD).

## Step 1: FastQC - Quality Control of Raw Reads
FastQC is used to assess the quality of raw sequencing reads before processing. This helps identify low-quality sequences, adapter contamination, and biases.

Command:
```
fastqc -t 8 -o fastqc_reports/ /home/s2451842/projects/HAoECs_RNA_seq/data/*.fq.gz
```

### Outputs:
- **HTML and text-based reports** on quality scores, GC content, sequence length distribution, and other quality metrics.
- A MultiQC report summarizing these results from all samples.

## Data Preprocessing: Adapter Trimming
Trimmomatic is used to remove sequencing adapters and low-quality reads to improve mapping accuracy.

Command:
```
trimmomatic PE -threads 8 \
    -phred33 \
    /home/s2451842/projects/HAoECs_RNA_seq/data/C1_311_1.fq.gz \
    /home/s2451842/projects/HAoECs_RNA_seq/data/C1_311_2.fq.gz \
    -baseout trimmed_data/C1_311 \
    ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```

### Outputs:
- Trimmed FASTQ files for each sample
- Log files indicating the number of reads removed

## Read Alignment using HISAT2
HISAT2 is a fast and memory-efficient aligner used to map RNA-Seq reads to the reference genome.

Command:
```
hisat2 -x Homo_sapiens.GRCh38_index \
    -1 trimmed_data/C1_311_1P.fq.gz \
    -2 trimmed_data/C1_311_2P.fq.gz \
    -S REVERSE_Hisat2_HAoECs/C1_311.sam \
    --rna-strandness RF \
    --dta -p 8
```

### Outputs:
- Aligned reads in **SAM format**
- Summary statistics from HISAT2, stored in logs

## BAM File Processing
After alignment, SAM files are converted into BAM format, sorted, and indexed using Samtools. These steps optimize data for further analysis.

Commands:
```
samtools view -bS REVERSE_Hisat2_HAoECs/C1_311.sam > REVERSE_Hisat2_HAoECs/C1_311.bam
samtools sort -@ 8 -o REVERSE_Hisat2_HAoECs/C1_311.sorted.bam REVERSE_Hisat2_HAoECs/C1_311.bam
samtools index REVERSE_Hisat2_HAoECs/C1_311.sorted.bam
```

### Outputs:
- **Unsorted BAM file**: Contains aligned sequence reads.
- **Sorted BAM file**: Required for downstream analysis.
- **Index file** for quick access to alignments.

## Gene Quantification using HTSeq
HTSeq is used to count the number of reads mapped to each gene. The `-s reverse` flag is used for libraries prepared with reverse-stranded protocols.

Command:
```
htseq-count -f bam -r pos -s reverse -t exon -i gene_id \
    REVERSE_Hisat2_HAoECs/C1_311.sorted.bam \
    /datastore/home/genomes/igenomes/Homo_sapiens/annotation.gtf \
    > REVERSE_Hisat2_HAoECs/C1_311.counts.txt
```

### Outputs:
- **Gene count files** (`.txt` files) for each sample, which can be used for differential expression analysis.

## Aggregating QC Reports with MultiQC
MultiQC is used to summarize various quality control reports into a single, easy-to-read HTML file.

Command:
```
multiqc -o multiqc_reports qc_reports REVERSE_Hisat2_HAoECs logs
```

### Outputs:
- **MultiQC Report**: A single interactive HTML report summarizing read quality, adapter content, alignment quality, and more.

## Final Output and Next Steps
1. **Aligned Reads**: BAM files in `REVERSE_Hisat2_HAoECs/`
2. **Sorted and Indexed BAM files**: `REVERSE_Hisat2_HAoECs/*.sorted.bam`
3. **Gene Count Matrix**: `REVERSE_Hisat2_HAoECs/*.counts.txt`
4. **Quality Control Reports**: `qc_reports/`, `multiqc_reports/`

### Next Steps:
With these results, the next steps include:
- **Quality Filtering and Analysis**: Assess the quality using MultiQC.
- **Differential Gene Expression Analysis**: Use DESeq2 or edgeR in R to compare gene expression levels between control and treatment groups.
- **Functional Enrichment Analysis**: Identify biological pathways affected by the treatment.

This guide provides a foundation for beginners to understand and implement an RNA-Seq pipeline using HISAT2 for alignment and HTSeq for quantification. For professional bioinformaticians, further optimizations, statistical modeling, and functional analysis can be conducted using additional tools like DESeq2, edgeR, and sleuth.




**Alternative Explanation: Title: Step-by-Step Explanation of RNA-Seq Analysis Using HISAT2 and HTSeq**

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







