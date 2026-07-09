9  Variant Calling



TipLearning Objectives

* Explain the statistical principles behind calling a variant from a pileup of reads.
* Use `BCFtools` to identify SNPs and Indels in a sample.
* Understand the structure of a VCF file and its key metadata fields.
* Apply basic filters to distinguish true biological variants from sequencing artifacts.

## 9.1 Variant Calling

**Variant Calling** is the process of identifying specific differences between your sequenced sample and the reference genome. These differences are usually **SNPs** (Single Nucleotide Polymorphisms), **Insertions**, or **Deletions** (collectively known as **Indels**).

### 9.1.1 How Variant Callers Work

A variant caller looks at the “pileup” of reads at every single position in the genome. It uses statistical models to decide if a mismatch is a **true biological variant** or just a **sequencing error**.

To make this decision, the software considers:

* **Read Depth**: How many reads cover that position?
* **Base Quality**: Are the mismatching bases high-quality or “noisy”?
* **Mapping Quality**: Do the reads belong at this location, or are they multi-mappers?
* **Allele Frequency**: In a diploid organism (like humans), we expect a heterozygous variant to be present in roughly 50% of the reads.

### 9.1.2 Common Tools

There are several standard tools for calling variants:

* **BCFtools**: Fast, lightweight, and great for standard germline variants.
* **GATK (Genome Analysis Toolkit)**: The “Gold Standard” for human clinical data, though more complex to run.
* **FreeBayes**: A Bayesian variant caller that handles complex polyploid organisms well.

WarningVariant calling can be complicated

For the purposes of this course we are showing you a very simple way to call and filter variants that is more likely to be appropriate when working with ‘simpler’ organisms such as viruses and bacteria. In reality, tools like `GATK` and `FreeBayes` are what researchers studying eukaryotes use for variant calling.

### 9.1.3 Workflow with BCFtools

`BCFtools` is a set of utilities for variant calling and that also manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed.

A common way to call variants is a two-step process using `BCFtools`:

```
# 1. Generate a pileup (mpileup) of all positions
bcftools mpileup -f genome.fasta alignments.sorted.bam > raw.bcf

# 2. Call the variants (determine the genotypes)
bcftools call -mv -Ov -o variants.vcf raw.bcf
```

Alternatively the two commands can be run together using the pipe `|` operator:

```
bcftools mpileup -f genome.fasta alignments.sorted.bam | bcftools call -mv -Ov > variants.vcf
```

Key Arguments:

* `-f`: The reference genome FASTA.
* `-m`: Allows for multi-allelic calling.
* `-v`: Output only the variant sites (skip the positions that match the reference).
* `-Ov`: Output the result in VCF format.

### 9.1.4 Filtering Variants

Not all variants in your initial VCF/BCF file are “real.” Many are false positives caused by alignment artifacts. We often filter our VCF files based on:

* **Depth (DP)**: Require at least 10x coverage.
* **Quality (QUAL)**: Require a Phred-scaled quality score > 30.

```
bcftools filter -s LOWQUAL -i 'QUAL>30 && DP>10' -Ov variants.vcf > filtered_variants.vcf
```

ExerciseExercise 1 - Exercise: Identifying Variants

1. Run the `bcftools mpileup` and `call` pipeline on your sorted yeast BAM file.
2. Filter the vcf file with `bcftools filter`.
3. Use `bcftools view`to look at the resulting VCF file.
4. **Question:** Look at the first few rows of your VCF. What is the most common type of mutation you found (SNP, Insertion, or Deletion)?

AnswerAnswer

1. The pileup and call command we ran:

```
# create an output directory
mkdir -p results/variants

# run the mpileup and call command
bcftools mpileup -f resources/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  results/mapping/SRR36763845.sorted.bam | bcftools call -mv -Ov > results/variants/SRR36763845_variants.vcf
```

2. Then we filtered our variants:

```
bcftools filter -s LOWQUAL -i 'QUAL>30 && DP>10' -Ov results/variants/SRR36763845_variants.vcf > results/variants/SRR36763845_filtered_variants.vcf
```

3. We opened the vcf file with `bcftools view` using the pipe `|` and `less -S` so we can scroll through the file:

```
bcftools view results/variants/SRR36763845_filtered_variants.vcf | less -S
```

2. **Observation:** We see many more SNPs (single base changes) than Indels.
3. **VCF check:** In the REF and ALT columns, a SNP looks like `A -> G`, while an Indel looks like `A -> ATCC` or `GTT -> G`.

## 9.2 Summary

TipKey Points

* Variant Calling identifies SNPs and Indels by comparing the pileup of reads against the reference sequence.
* `BCFtools` is an efficient and widely used suite for calling and manipulating variants.
* The VCF (Variant Call Format) file is the standard output, containing coordinates, alleles, and quality metrics.
* Filtering is a crucial post-calling step to remove low-confidence calls caused by mapping errors or low coverage.
* Variant calling requires a sorted and indexed BAM file as input.