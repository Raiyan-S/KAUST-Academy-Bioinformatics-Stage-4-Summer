8  Sequence Alignment (Mapping)



TipLearning Objectives

* Explain the process of mapping raw sequencing reads to a reference genome.
* Understand the difference between global and local alignment algorithms used by tools like `BWA`.
* Convert human-readable alignment files (SAM) to efficient binary formats (BAM).
* Perform coordinate sorting and indexing of alignment files for downstream analysis.

## 8.1 The DNA-seq Analysis Pipeline Workflow

We have now completed the **Reference Preparation** (the red box) part of our workflow and can move onto the next step which is **Alignment** (yellow).

[![The DNA-seq Analysis Pipeline Workflow](images/week_2_workflow.png)](images/week_2_workflow.png "The DNA-seq Analysis Pipeline Workflow")

The DNA-seq Analysis Pipeline Workflow

## 8.2 Sequence Alignment (Mapping)

Now that we have **clean reads** (FASTQ) and an **indexed reference** (FASTA), we can perform **read alignment** also referred to as **read mapping**. The goal is to find the best coordinates for each read in the genome, accounting for potential mismatches or small indels (insertions/deletions).

[![Diagram illustrating the steps involved in mapping sequencing reads to a reference genome. Mapping programs allow some differences between the reads and the reference genome (red mutation shown as an example). Before doing the mapping, reads are usually filtered for high-quality and to remove any sequencing adapters. The reference genome is also indexed before running the mapping step. The mapped file (BAM format) can be used in many downstream analyses. See text for more details.](images/ngs_mapping.svg)](images/ngs_mapping.svg "Diagram illustrating the steps involved in mapping sequencing reads to a reference genome. Mapping programs allow some differences between the reads and the reference genome (red mutation shown as an example). Before doing the mapping, reads are usually filtered for high-quality and to remove any sequencing adapters. The reference genome is also indexed before running the mapping step. The mapped file (BAM format) can be used in many downstream analyses. See text for more details.")

Diagram illustrating the steps involved in mapping sequencing reads to a reference genome. Mapping programs allow some differences between the reads and the reference genome (red mutation shown as an example). Before doing the mapping, reads are usually filtered for high-quality and to remove any sequencing adapters. The reference genome is also indexed before running the mapping step. The mapped file (BAM format) can be used in many downstream analyses. See text for more details.

### 8.2.1 How Aligners Work

Aligners like **`BWA` (Burrows-Wheeler Aligner)** or **`Bowtie2`** use the index files we created earlier to perform a rapid search. Instead of checking every base, they use “seeds” (short chunks of the read) to narrow down the possible locations before performing a more detailed alignment.

### 8.2.2 Performing Alignment with BWA-MEM

For most DNA-seq data (Whole Genome or Exome), **`BWA-MEM 2`** is one of the most popular tools. It is robust to sequencing errors and handles longer reads well.

A typical alignment command looks like this:

```
bwa-mem2 mem -t 4 genome.fasta R1_trimmed.fq.gz R2_trimmed.fq.gz > alignments.sam
```

Key Arguments:

* `mem`: The specific algorithm used (Maximal Exact Matches).
* `-t 4`: Tells the program to use 4 CPU cores to speed up the process.
* `genome.fasta`: The path to your reference (BWA will automatically look for the index files with the same name).
* `R1/R2`: Your forward and reverse trimmed reads.
* `> alignments.sam`: Redirects the output from the screen into a SAM file.

### 8.2.3 Post-Processing: From SAM to BAM

As we learned in the file formats section, SAM files are huge and slow. We almost always immediately convert the output to BAM, sort it by genomic coordinates, and index it. Most downstream tools (like variant callers) require the reads to be ordered by their position on the chromosome (e.g., all reads for Chromosome 1 first, then Chromosome 2) so they can process the genome linearly.

An index for the sorted BAM file is often required for downstream analysis and for visualising the alignment with programs such as the integrated genome viewer (IGV) or Artemis.

```
# 1. Convert SAM to BAM and Sort
samtools sort -@ 4 -o alignments.sorted.bam alignments.sam

# 2. Index the BAM file
samtools index alignments.sorted.bam
```

ExerciseExercise 1 - Exercise: Running an Alignment

1. Index the reference sequence with `bwa-mem2 index`
2. Align the trimmed yeast reads to the indexed yeast genome using `bwa-mem2 mem`.
3. Convert the resulting SAM file into a sorted BAM file using `samtools`.
4. Index your final BAM file.
5. **Question:** Use `ls -lh` to compare the size of your `.sam` file and your `.sorted.bam` file in the `results/mapping` directory. Which one is smaller, and by roughly how much?

AnswerAnswer

1. We ran the following commands:

```
# index the reference with bwa-mem2 index
bwa-mem2 index resources/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

# create an output directory
mkdir -p results/mapping

# map the reads
bwa-mem2 mem resources/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  results/qc/SRR36763845_1_trimmed.fastq.gz results/qc/SRR36763845_2_trimmed.fastq.gz > results/mapping/SRR36763845.sam 

# sort and convert sam to bam
samtools sort -o results/mapping/SRR36763845.sorted.bam results/mapping/SRR36763845.sam

# index the bam file
samtools index results/mapping/SRR36763845.sorted.bam
```

2. **Size Comparison**: You can see that the BAM file is roughly **six** times smaller than the SAM file. This is because BAM uses binary compression (BGZF), while SAM is plain text.

```
total 2.9G
-rw-r--r-- 2 ajv37 2.5G Feb 25 09:55 SRR36763845.sam
-rw-r--r-- 2 ajv37 388M Feb 25 09:59 SRR36763845.sorted.bam
-rw-r--r-- 2 ajv37  38K Feb 25 09:59 SRR36763845.sorted.bam.bai
```

ExerciseExercise 2 - Exercise: Calculate the depth of coverage (Bonus exercise)

1. Install the `qualimap` tool in the `intro_ngs` environment with `mamba install -c bioconda qualimap`
2. Run `qualimap` on the `SRR36763845.sorted.bam` file and save the results to the `results/mapping` directory.
3. Open `qualimapReport.html` in your web browser.
4. **Question:** Examine the ‘Coverage across Reference’ plot. What do you think is causing the very high peaks you see across the reference genome?

AnswerAnswer

1. We ran the following commands:

```
# install qualimap in the intro_ngs environment
mamba install -c bioconda qualimap

# run qualimap
qualimap bamqc -bam results/mapping/SRR36763845.sorted.bam -outdir results/qualimap
```

2. We opened the `qualimapReport.html` file in our web browser

[![qualimap Depth Coverage](images/qualimap_depth_coverage.png)](images/qualimap_depth_coverage.png "qualimap Depth Coverage")

qualimap Depth Coverage

3. The peaks we see are likely due to there being repeat regions in different parts of the reference genome. This means the aligner does not know which region to align the read to so it often ends up stacking the reads in one repeat instead of the correct repeat region which may be in another part of the genome.

### 8.2.4 Visualising BAM Files in IGV

One thing that can be useful is to visualise the BAM files we produce and we can use the program *IGV* (Integrative Genome Viewer) to do this:

* Open *IGV* and go to `File → Load from file…`.
* In the file browser that opens go to the folder `results/mapping/` and select the file `SRR36763845.sorted.bam` to open it.

There are several ways to search and browse through our alignments, exemplified in the figure below.

[![Screenshot IGV program. The search box at the top can be used to go to a specific region in the format “CHROM:START-END”.](images/igv_overview.svg)](images/igv_overview.svg "Screenshot IGV program. The search box at the top can be used to go to a specific region in the format “CHROM:START-END”.")

Screenshot IGV program. The search box at the top can be used to go to a specific region in the format “CHROM:START-END”.

## 8.3 Summary

TipKey Points

* Alignment is the process of finding where a short read fits best on a long reference sequence.
* `BWA-MEM 2` is a standard tool for aligning DNA-seq reads to a reference genome.
* SAM files are the text-based output of aligners; BAM files are their compressed, binary equivalents.
* Sorting and Indexing are mandatory post-alignment steps that allow software to jump to specific genomic regions instantly.