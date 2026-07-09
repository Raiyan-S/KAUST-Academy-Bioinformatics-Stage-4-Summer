4  Common file formats



TipLearning Objectives

* Recognise the structure of common file formats in bioinformatics, in particular FASTQ, FASTA, GFF3 and CSV/TSV.
* To understand which file type corresponds to the outputs of a given analysis

## 4.1 Introduction

This page lists some common file formats used in Bioinformatics. The heading of each file links to a page with more details about each format.

Generally, files can be classified into two categories: text files and binary files.

* **Text files** can be opened with standard text editors, and manipulated using command-line tools (such as `head`, `less`, `grep`, `cat`, etc.). However, many of the standard files listed in this page can be opened with specific software that displays their content in a more user-friendly way.
* **Binary files** are often used to store data more efficiently. Typically, specific tools need to be used with those files. For example, the BAM format is used to store sequences aligned to a reference genome and can be manipulated with dedicated software such as `samtools`.

Before diving into the technical details of each format, it is helpful to have a “quick-glance” guide. In bioinformatics, we move through a specific hierarchy of files, from raw data to refined biological insights.

### 4.1.1 Common Bioinformatics File Formats: Quick Reference

The following table summarizes the files you will encounter most frequently in an NGS pipeline and the ones we will use in this course:

| Format | Full Name | Type | Primary Usage in Pipeline |
| --- | --- | --- | --- |
| **CSV/TSV** | Comma/Tab Separated Values | Text | Storing general tabular data, metadata, or software summaries. |
| **FASTQ** | Fast-Q (Quality) | Text (Compressed) | **Raw Data:** Contains the actual sequences from the sequencer plus quality scores. |
| **FASTA** | Fast-A | Text | **The Map:** Stores reference genomes or individual protein/nucleotide sequences. |
| **SAM/BAM** | Sequence Alignment Map | Text / Binary | **The Alignment:** Shows exactly where each read (from FASTQ) fits on the reference (FASTA). |
| **VCF** | Variant Call Format | Text / Binary | **The Result:** Lists the differences (mutations) between your sample and the reference. |
| **BED** | Browser Extensible Data | Text | **Coordinates:** Defines specific regions of interest (e.g., “look between base 100 and 500”). |
| **GFF/GTF** | General Feature Format | Text | **Annotations:** Highly detailed map of where genes, exons, and transcripts are located. |

### 4.1.2 [CSV](https://en.wikipedia.org/wiki/Comma-separated_values#Example) and [TSV](https://en.wikipedia.org/wiki/Tab-separated_values#Example)

**Comma-separated values** (CSV) and **tab-separated values** (TSV) files are text-based formats commonly used to store **tabular data**. While strictly not specific to bioinformatics, they are commonly used as the output of bioinformatic software. CSV files usually have `.csv` extension, while TSV files often have `.tsv` or the more generic `.txt` extension.

In both cases, the data is organized into rows and columns. Rows are represented across different lines of the file, while the columns are separated using a **delimiting character**: a command `,` in the case of CSV files and a tab space (`tab ↹`) for TSV files.

For example, for this table:

| sample | date | strain |
| --- | --- | --- |
| VCH001 | 2023-08-01 | O1 El Tor |
| VCH002 | 2023-08-02 | O1 Classical |
| VCH003 | 2023-08-03 | O139 |
| VCH004 | 2023-08-04 | Non-O1 Non-O139 |

This would be its representation as a CSV file:

```
sample,date,strain
VCH001,2023-08-01,O1 El Tor
VCH002,2023-08-02,O1 Classical
VCH003,2023-08-03,O139
VCH004,2023-08-04,Non-O1 Non-O139
```

And this is its representation as a TSV file (the space between columns is a `tab ↹`):

```
sample    date        strain
VCH001    2023-08-01  O1 El Tor
VCH002    2023-08-02  O1 Classical
VCH003    2023-08-03  O139
VCH004    2023-08-04  Non-O1 Non-O139
```

CSV and TSV files are human-readable and can be opened and edited using **basic text editors** or **spreadsheet software** like *Microsoft Excel*.

### 4.1.3 [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)

FASTQ files are used to store **nucleotide sequences along with a quality score** for each nucleotide of the sequence. These files are the typical format **obtained from NGS sequencing** platforms such as Illumina and Nanopore (after basecalling). Common file extensions used for this format include `.fastq` and `.fq`.

The file format is as follows:

```
@SEQ_ID                   <-- SEQUENCE NAME
AGCGTGTACTGTGCATGTCGATG   <-- SEQUENCE
+                         <-- SEPARATOR
%%).1***-+*''))**55CCFF   <-- QUALITY SCORES
```

In FASTQ files each sequence is always represented across 4 lines. The quality scores are encoded in a compact form, using a single character. They represent a score that can vary between 0 and 40 (see [Illumina’s Quality Score Encoding](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)). The reason single characters are used to encode the quality scores is that it saves space when storing these large files. Software that work on FASTQ files automatically convert these characters into their score, so we don’t have to worry about doing this conversion ourselves.

The quality value in common use is called a **Phred score** and it represents the **probability that the base is an error**. For example, a base with quality 20 has a probability \(10^{-2} = 0.01 = 1\%\) of being an error. A base with quality 30 has \(10^{-3} = 0.001 = 0.1\%\) chance of being an error. Typically, a Phred score threshold of >20 or >30 is used when applying quality filters to sequencing reads.

Because FASTQ files tend to be quite large, they are **often compressed** to save space. The most common compression format is called *gzip* and uses the extension `.gz`. To look at a *gzip* file, we can use the command `zcat`, which decompresses the file and prints the output as text.

For example, we can use the following command to count the number of lines in a compressed FASTQ file:

```
zcat sequences.fq.gz | wc -l
```

If we want to know how many sequences there are in the file, we can divide the result by 4 (since each sequence is always represented across four lines).

### 4.1.4 [FASTA](https://en.wikipedia.org/wiki/FASTA)

FASTA files are used to store **nucleotide or amino acid sequences**. Common file extensions used for this format include `.fasta`, `.fa`, `.fas` and `.fna`.

The general structure of a FASTA file is illustrated below:

```
>sample01                 <-- NAME OF THE SEQUENCE
AGCGTGTACTGTGCATGTCGATG   <-- SEQUENCE ITSELF
```

Each sequence is represented by a name, which always starts with the character `>`, followed by the actual sequence.

A FASTA file can contain several sequences, for example:

```
>sample01
AGCGTGTACTGTGCATGTCGATG
>sample02
AGCGTGTACTGTGCATGTCGATG
```

Each sequence can sometimes span multiple lines, and separate sequences can always be identified by the `>` character. For example, this contains the same sequences as above:

```
>sample01      <-- FIRST SEQUENCE STARTS HERE
AGCGTGTACTGT
GCATGTCGATG
>sample02      <-- SECOND SEQUENCE STARTS HERE
AGCGTGTACTGT
GCATGTCGATG
```

To count how many sequences there are in a FASTA file, we can use the following command:

```
grep ">" sequences.fa | wc -l
```

In two steps:

* find the lines containing the character “>”, and then
* count the number of lines of the result.

FASTA files are commonly used to **store genome sequences**, after they have been assembled. We will see FASTA files several times throughout these materials, so it’s important to be familiar with them.

#### Specific filename extensions

The generic form of FASTA file has the `.fas` extension. For more specific types, we can use the following:

| Extension | Meaning | Notes |
| --- | --- | --- |
| fna | FASTA nucleic acid | Used generically to specify nucleic acids |
| ffn | FASTA nucleotide coding regions | Contains coding regions for a genome |
| faa | FASTA amino acid | Contains amino acid sequences. A multiple protein fasta file can have the more specific extension mpfa |
| frn | FASTA non-coding RNA | Contains non-coding RNA regions for a genome, in DNA alphabet e.g. tRNA, rRNA |

### 4.1.5 [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) & [BAM](https://en.wikipedia.org/wiki/SAM_(file_format)#BAM)

SAM and BAM files are used to store **sequence reads aligned to a reference genome**. Common file extensions include `.sam` for the text version and `.bam` for the binary version.

The general structure of a SAM file is illustrated below:

```
@HD VN:1.6  SO:coordinate                              <-- HEADER SECTION
@SQ SN:chr1  LN:248956422
read01  99  chr1  10023  60  50M  =  10080  107  AGCT...  <-- ALIGNMENT LINE
```

A SAM file is a tab-delimited text file where each alignment line contains 11 mandatory fields, including:

* QNAME: The ID of the read.
* FLAG: A numeric code representing alignment information (e.g., is it paired?).
* RNAME: The name of the reference sequence (from your FASTA).
* POS: The left-most position where the read aligns.
* CIGAR: A string describing “edits” (matches, insertions, deletions).

#### The Difference Between SAM and BAM

SAM files are human-readable text files but can be incredibly large. BAM files are the binary, compressed version of SAM. They are much smaller and “indexable,” allowing software to jump to specific genomic coordinates instantly. Because BAM files are binary, you cannot read them with standard tools like cat or head. You must use a tool like `samtools`.

To view the content of a BAM file, we use:

```
samtools view alignments.bam | head
```

To count how many alignments (reads) are in a BAM file:

```
samtools view -c alignments.bam
```

In this command:

* `samtools view` opens the compressed file.
* `-c` tells the program to count the records instead of printing them.

SAM/BAM files are the standard format for mapping, variant calling, and RNA-seq analysis. You will almost always work with the BAM version to save disk space and processing time.

### 4.1.6 [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)

VCF (Variant Call Format) files are used to store **gene sequence variations** (like SNPs, insertions, and deletions) relative to a reference genome. Common file extensions for this format include `.vcf` and its compressed binary version `.bcf`.

The general structure of a VCF file consists of a header and data lines:

```
##fileformat=VCFv4.2                     <-- METADATA (Header)
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM POS ID REF ALT QUAL FILTER INFO   <-- COLUMN HEADERS
chr1 10177 . A C 100 PASS AC=1;AF=0.5    <-- VARIANT RECORD
```

A VCF file is a tab-delimited text file. After the header (lines starting with `##`), each row represents a single variant with mandatory columns:

* **CHROM**: The chromosome where the variant occurs.
* **POS**: The 1-based reference position.
* **REF**: The reference allele (the “normal” base).
* **ALT**: The alternate allele (the “mutated” base).
* **QUAL**: A Phred-scaled quality score for the assertion of a variant.
* **INFO**: Additional information like ancestral alleles or population frequencies.

#### Text vs. Compressed

* **VCF** is a plain text file that can be inspected with `less` or `grep`.
* **BCF** is the binary equivalent, designed for high-speed processing and efficiency.

Because VCF files can be massive, they are often compressed with `bgzip` and indexed with `tabix` to allow software to quickly find variants in a specific genomic region.

To view the variants without the long header section, we can use:

```
grep -v "^##" variants.vcf | head
```

To count how many variants (mutations) are stored in a VCF file:

```
grep -v "^#" variants.vcf | wc -l
```

In two steps:

* find the lines that do not start with the character “#” (excluding the header), and then
* count the number of lines of the result.

VCF files are the **final output of a variant calling pipeline**. They are used in clinical diagnostics and population genetics to identify the functional impact of specific mutations.

### 4.1.7 [BED](https://en.wikipedia.org/wiki/BED_(file_format))

BED (Browser Extensible Data) files are used to store **genomic regions or features**, such as genes, exons, or transcription factor binding sites. Common file extensions for this format include `.bed`.

The general structure of a BED file is illustrated below:

```
chr1    10000   10500   geneA    1000    +    <-- GENOMIC REGION
chr1    20000   21000   geneB    500     -    <-- GENOMIC REGION
```

A BED file is a tab-delimited text file. It is unique because it uses a **0-based, half-open coordinate system** (the start position is 0-indexed, but the end position is not included in the range).

There are three mandatory fields for every BED file: \* **chrom**: The name of the chromosome (e.g., `chr1`). \* **chromStart**: The starting position of the feature. \* **chromEnd**: The ending position of the feature.

Up to 9 additional optional fields can be included, such as **name**, **score**, and **strand** (`+` or `-`).

#### BED vs. VCF/SAM

* **BED** is much simpler than VCF or SAM; it is designed to define “where” something is, rather than “what” sequence or variation is there.
* It is the primary format used for **genome browsers** (like UCSC or IGV) to display tracks of data.

To count how many genomic regions are defined in a BED file, we can use:

```
wc -l regions.bed
```

Because BED files typically do not have a header, a simple line count usually gives the total number of features.

To find all regions on a specific chromosome (e.g., Chromosome 1):

```
grep "^chr1" regions.bed | wc -l
```

In two steps:

* find the lines starting with “chr1”, and then
* count the number of lines in that result.

BED files are essential for intersecting data—for example, finding which DNA mutations (from a VCF) fall inside specific gene boundaries (from a BED).

### 4.1.8 [GFF](https://en.wikipedia.org/wiki/General_feature_format) / [GTF](https://en.wikipedia.org/wiki/Gene_transfer_format)

GFF (General Feature Format) files are used to store **genomic annotations**, specifically describing the exact locations of genes, transcripts, and exons. Common file extensions include `.gff`, `.gff3`, and the closely related `.gtf`.

The general structure of a GFF file is illustrated below:

```
##gff-version 3                                  <-- HEADER
chr1  HAVANA  gene  11869  14409  .  +  .  ID=gene01 <-- GENOMIC FEATURE
chr1  HAVANA  exon  11869  12227  .  +  .  Parent=gene01
```

A GFF file is a tab-delimited text file consisting of 9 mandatory columns:

* **seqid**: The name of the chromosome or scaffold.
* **source**: The program or database that generated the feature (e.g., “GenBank”).
* **type**: The type of feature (e.g., `gene`, `exon`, `CDS`, `mRNA`).
* **start / end**: The 1-based coordinates of the feature.
* **score / strand**: Numerical confidence and the orientation (`+` or `-`).
* **phase**: Used for coding sequences to indicate where the next codon begins.
* **attributes**: A semicolon-separated list of tag-value pairs (e.g., `ID=gene123;Name=BRCA1`).

#### GFF vs. BED

* **GFF** uses **1-based** coordinates, whereas BED uses 0-based coordinates.
* **GFF** is hierarchical; it can show that an “exon” belongs to a specific “mRNA,” which belongs to a specific “gene.” BED files are usually “flat” lists of regions.

Because GFF files often contain a mix of different feature types, we frequently use `awk` or `grep` to isolate specific parts.

To count how many **genes** are defined in a GFF file, we can use:

```
grep -v "^#" annotations.gff | awk '$3=="gene"' | wc -l
```

In three steps:

* ignore the header lines (starting with #),
* find lines where the 3rd column is exactly “gene”, and then
* count those lines.

GFF files are the gold standard for genome annotation. When you download a reference genome in FASTA format, you almost always download the corresponding GFF file to know where the actual genes are located.

## 4.2 Summary

TipKey Points

* Common file formats in bioinformatics include FASTQ, FASTA and GFF. These are all text-based formats.
* FASTQ format (`.fastq` or `.fq`):
  + Designed to store sequences along with quality scores.
  + Contains a sequence identifier, sequence data, a separator line and quality scores.
  + Widely used for storing **sequence reads generated by NGS platforms**.
* FASTA format (`.fasta`, `.fa`, `.fas`, `.fna`):
  + Is used for storing biological sequences, including DNA, RNA, and protein.
  + It Comprises a sequence identifier (often preceded by “>”) and the sequence data.
  + Commonly used for sequence storage and exchange of **genome sequences**.
* GFF format (`.gff` or `.gff3`):
  + A structured, tab-delimited format for describing genomic features and annotations.
  + Consists of nine standard columns, including sequence identifier, feature type, start and end coordinates, strand information, and attributes.
  + Facilitates the representation of genes, transcripts, and other genomic elements, supporting hierarchical structures and metadata.
  + Commonly used for storing and sharing **genomic annotation** data in bioinformatics.
* CSV (`.csv`) and TSV (`.tsv`):
  + Plain text formats to store tables.
  + The columns in the CSV format are delimited by comma, whereas in the TSV format by a tab.
  + These files can be opened in standard spreadsheet software such as *Excel*.

#### References

https://github.com/cambiotraining/sars-cov-2-genomics

https://snipcademy.com/sequence-file-formats

Wellcome Genome Campus Advanced Courses and Scientific Conferences 2017 - WORKING WITH PATHOGEN GENOMES Course Manual http://www.wellcome.ac.uk/advancedcourses