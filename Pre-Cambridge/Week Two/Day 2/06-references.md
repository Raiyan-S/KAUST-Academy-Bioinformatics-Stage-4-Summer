7  Reference genomes



TipLearning Objectives

* Identify the primary databases used to download reference genomes and annotations.
* Differentiate between a genome sequence (FASTA) and a genome annotation (GFF/GTF).
* Explain the “Matching Principle” for maintaining data integrity in bioinformatics pipelines.
* Understand the purpose of indexing and how it facilitates rapid data access for alignment tools.

## 7.1 The DNA-seq Analysis Pipeline Workflow

In the previous two chapters we did the routine work of checking the quality of our sequence data with `FastQC` and trimming the poor quality reads with `fastp`. We have now completed the **Data Pre-processing** (blue in the figure below) part of our workflow and can move onto the next step which is **Reference Preparation** (the red box).

[![The DNA-seq Analysis Pipeline Workflow](images/week_2_workflow.png)](images/week_2_workflow.png "The DNA-seq Analysis Pipeline Workflow")

The DNA-seq Analysis Pipeline Workflow

## 7.2 Reference Genomes and Annotations

To analyze NGS data, we usually need a **Reference Genome** to act as a map for our reads. Without this map, we are left with millions of individual fragments and no way to know which genes or regulatory regions they represent. A reference genome is a high-quality, representative sequence of the complete set of DNA for a specific species. It is typically stored in FASTA format and serves as a gold standard.

Think of it as the picture on the box of a 3-billion-piece jigsaw puzzle. Your sequencing reads are the individual puzzle pieces; the reference genome tells you exactly where each piece belongs.

A complete reference consists of two main parts: the **sequence** and the **annotation**.

### 7.2.1 Where to download Genomes (FASTA)

The genome sequence contains the actual nucleotides (\(A, C, G, T\)) for every chromosome. These are stored in **FASTA** format. Depending on your organism of interest, there are three primary “Gold Standard” databases:

* **[Ensembl](https://www.ensembl.org/info/data/ftp/index.html):** Preferred for vertebrate and human genetics. It provides well-curated versions of the genome.
* **[NCBI (GenBank)](https://www.ncbi.nlm.nih.gov/genome/):** The most comprehensive database, covering everything from viruses and bacteria to complex eukaryotes.
* **[UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html):** Widely used for visualization and providing specific “hg” versions of the human genome (e.g., `hg38`).

Note

When downloading a FASTA file, look for the **“Primary Assembly”** version. Avoid “Top-level” files if they contain “patches” or “haplotypes,” as these can cause issues during the alignment step.

### 7.2.2 Where to download Annotations (GFF/GTF)

The FASTA file tells us the sequence, but it doesn’t tell us where the genes are. For that, we need an **Annotation file**, usually in **GFF3** or **GTF** format.

These files provide the coordinates for: - **Genes** and **Transcripts** - **Exons** and **Introns** - **UTRs** (Untranslated Regions)

#### Which format should I choose?

* **GTF (Gene Transfer Format):** Highly standardized and preferred by most RNA-seq tools (like FeatureCounts or StringTie).
* **GFF3 (General Feature Format):** More flexible and hierarchical; often used for non-model organisms or plant genomes.

WarningCrucial Rule: The Matching Principle

Always download your FASTA and your GFF/GTF from the **same source** and the **same version**. For example, if you download the `GRCh38` FASTA from Ensembl, you **must** download the `GRCh38` GTF from Ensembl. If you mix NCBI sequences with Ensembl annotations, the chromosome names (e.g., `1` vs `chr1`) will not match, and your pipeline will fail.

ExerciseExercise 1 - Exercise: Finding a Reference

1. Go to the [Ensembl FTP Download](https://www.ensembl.org/info/data/ftp/index.html) site.
2. Search for *Saccharomyces cerevisiae* (Yeast).
3. Identify the link for the **DNA (FASTA)** and the **Gene annotation (GTF)**.
4. **Question:** What is the current assembly version name for the Yeast genome?

AnswerAnswer

1. In the Ensembl list, the assembly version for *Saccharomyces cerevisiae* is currently **R64-1-1**.

[![Ensembl page for Saccharomyces cerevisiae](images/saccharomyces_ensemble.png)](images/saccharomyces_ensemble.png "Ensembl page for Saccharomyces cerevisiae")

Ensembl page for *Saccharomyces cerevisiae*

2. The FASTA file name for the latest version: `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz`
3. The GTF file name for the latest version: `Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz`

## 7.3 Indexing Reference Files

Before we can use our FASTA and GFF/GTF files for alignment or analysis, they must be **indexed**.

Think of a human genome like a massive 3-billion-page book. If a sequencing tool (like an aligner) needs to find a specific sequence, it shouldn’t have to read the book from page 1 every single time. An **index** acts like the index at the back of a book, allowing the software to jump directly to the correct “page” (genomic coordinate).

### 7.3.1 Indexing the Sequence (FASTA)

There are two main types of indexing for FASTA files:

* **Simple Indexing (`.fai`):** Used by tools like `IGV` or `Samtools` to quickly jump to a specific chromosome or position.
* **Aligner Indexing:** Used by mappers like **`BWA`** or **`Bowtie2`**. These tools create a complex mathematical representation of the genome (often using a Burrows-Wheeler Transform) to make searching for millions of reads extremely fast.

How to create a simple index with Samtools:

```
samtools faidx genome.fasta
```

This creates a small file called `genome.fasta.fai`.

### 7.3.2 Indexing Annotations (GFF/GTF)

When working with large annotation files, we often compress them and create a `tabix` index. This allows tools to instantly extract features from a specific region (e.g., “Give me all exons on Chromosome 17”).

How to index a GFF3 file:

```
# 1. Block-compress the file
bgzip annotations.gff

# 2. Create the index
tabix -p gff annotations.gff.gz
```

This creates a `.tbi` file.

WarningDisk Space and Naming

Indexing often generates several extra files (e.g., `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`). **Do not rename or move these files** away from the original FASTA file. Most bioinformatics tools expect the index files to be in the same folder and have the exact same prefix as the reference.

ExerciseExercise 2 - Exercise: Indexing a Genome

1. Use `samtools` to index the yeast genome FASTA file you identified in the previous section. We’ve downloaded the file for you and it can be found in `resources/reference`.
2. Look at the directory contents of `resources/reference` using `ls -lh`.
3. **Question:** What is the file extension of the new file created? Is the index file larger or smaller than the original FASTA?

AnswerAnswer

1. We ran the following command:

```
samtools faidx resources/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
```

2. The command creates a `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai` file in `resources/reference`.
3. Running `ls -lh resources/reference` gives:

```
total 12M
-rw-r--r-- 1 ajv37 12M Feb 24 11:43 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
-rw-r--r-- 1 ajv37 416 Feb 24 11:46 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai
```

4. The `.fai` file is **much smaller** because it only contains the metadata about chromosome lengths and byte offsets, not the actual DNA sequence.

## 7.4 Summary

TipKey Points

* Reference genomes provide the “map” for NGS alignment and are typically stored in FASTA format.
* Annotations (GFF/GTF) describe the location and function of features like genes and exons.
* Matching sources (e.g., all from Ensembl) is critical to ensure chromosome naming conventions are consistent.
* Indexing is a mandatory step that allows software to navigate large genomic files efficiently.