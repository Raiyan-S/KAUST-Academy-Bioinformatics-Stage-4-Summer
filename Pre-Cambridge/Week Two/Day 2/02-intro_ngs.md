3  Introduction to NGS



TipLearning Objectives

* List the main high-throughput sequencing technologies in use.
* Describe the main differences between Illumina, Oxford Nanopore and PacBio platforms, including their advantages and disadvantages.

## 3.1 Next Generation Sequencing

The sequencing of genomes has become more routine due to the [rapid drop in DNA sequencing costs](https://www.genome.gov/about-genomics/fact-sheets/DNA-Sequencing-Costs-Data) seen since the development of Next Generation Sequencing (**NGS**) technologies in 2007. One main feature of these technologies is that they are *high-throughput*, allowing one to more fully characterise the genetic material in a sample of interest.

There are three main technologies in use nowadays, often referred to as 2nd and 3rd generation sequencing:

* Illumina’s sequencing by synthesis (2nd generation)
* Oxford Nanopore, shortened ONT (3rd generation)
* Pacific Biosciences, shortened PacBio (3rd generation)

The video below from the iBiology team gives a great overview of these technologies.

### 3.1.1 Illumina Sequencing

Illumina’s technology has become a widely popular method, with many applications to study transcriptomes (RNA-seq), epigenomes (ATAC-seq, BS-seq), DNA-protein interactions (ChIP-seq), chromatin conformation (Hi-C/3C-Seq), population and quantitative genetics (variant detection, GWAS), de-novo genome assembly, amongst [many others](https://moodle2.units.it/pluginfile.php/156101/mod_resource/content/0/sequencing-methods-review.pdf).

An overview of the sequencing procedure is shown in the animation video below. Generally, samples are processed to generate so-called *sequencing libraries*, where the genetic material (DNA or RNA) is processed to generate fragments of DNA with attached oligo adapters necessary for the sequencing procedure (if the starting material is RNA, it can be converted to DNA by a step of reverse transcription). Each of these DNA molecule is then sequenced from both ends, generating pairs of sequences from each molecule, i.e. *paired-end sequencing* (single-end sequencing, where the molecule is only sequenced from one end is also possible, although much less common nowadays).

This technology is a type of *short-read sequencing*, because we only obtain short sequences from the original DNA molecules. Typical protocols will generate 2x50bp to 2x250bp sequences (the 2x denotes that we sequence from each end of the molecule).

The main advantage of Illumina sequencing is that it produces very high-quality sequence reads (current protocols generate reads with an error rate of less than <1%) at a low cost. However, the fact that we only get relatively short sequences means that there are limitations when it comes to resolving particular problems such as long sequence repeats (e.g. around centromeres or transposon-rich areas of the genome), distinguishing gene isoforms (in RNA-seq), or resolving haplotypes (combinations of variants in each copy of an individual’s diploid genome).

**In summary, Illumina:**

* Utilizes sequencing-by-synthesis chemistry.
* Offers short read lengths.
* Known for high accuracy with low error rates (<1%).
* Well-suited for applications like DNA resequencing and variant detection.
* Scalable and cost-effective for large-scale projects.
* Limited in sequencing long DNA fragments.
* Expensive to set up.

### 3.1.2 Nanopore Sequencing

Nanopore sequencing is a type of *long-read sequencing* technology. The main advantage of this technology is that it can sequence very long DNA molecules (up to megabase-sized), thus overcoming the main shortcoming of short-read sequencing mentioned above. Another big advantage of this technology is its portability, with some of its devices designed to work via USB plugged to a standard laptop. This makes it an ideal technology to use in situations where it is not possible to equip a dedicated sequencing facility/laboratory (for example, when doing field work).

[![Overview of Nanopore sequencing showing the highly-portable MinION device. The device contains thousands of nanopores embedded in a membrane where current is applied. As individual DNA molecules pass through these nanopores they cause changes in this current, which is detected by sensors and read by a dedicated computer program. Each DNA base causes different changes in the current, allowing the software to convert this signal into base calls.](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41587-021-01108-x/MediaObjects/41587_2021_1108_Fig1_HTML.png?as=webp)](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41587-021-01108-x/MediaObjects/41587_2021_1108_Fig1_HTML.png?as=webp "Overview of Nanopore sequencing showing the highly-portable MinION device. The device contains thousands of nanopores embedded in a membrane where current is applied. As individual DNA molecules pass through these nanopores they cause changes in this current, which is detected by sensors and read by a dedicated computer program. Each DNA base causes different changes in the current, allowing the software to convert this signal into base calls.")

Overview of Nanopore sequencing showing the highly-portable MinION device. The device contains thousands of nanopores embedded in a membrane where current is applied. As individual DNA molecules pass through these nanopores they cause changes in this current, which is detected by sensors and read by a dedicated computer program. Each DNA base causes different changes in the current, allowing the software to convert this signal into base calls.

However, optimising this technology presents some challenges, notably in the production of sequencing libraries containing high molecular weight and intact DNA. It’s important to note that nanopore sequencing historically exhibited **higher error rates**, approximately 5% for older chemistries, compared to Illumina sequencing. However, significant advancements have emerged, enhancing the accuracy of nanopore sequencing technology, now achieving [accuracy rates exceeding 99%](https://nanoporetech.com/accuracy).

**In summary, ONT:**

* Operates on the principle of nanopore technology.
* Provides long read lengths, ranging from thousands to tens of thousands of base pairs.
* Ideal for applications requiring long-range information, such as *de novo* genome assembly and structural variant analysis.
* Portable, enabling fieldwork and real-time sequencing.
* Exhibits higher error rates (around 5%), with improvements in recent versions.
* Costs can be higher per base, compared to Illumina for certain projects.

### 3.1.3 PacBio Sequencing

Pacific Biosciences (PacBio) has pioneered a technology known as Single Molecule, Real-Time (SMRT) sequencing. Unlike Illumina, which relies on amplifying clusters of DNA and sequencing them in short bursts, PacBio sequences a single, long molecule of DNA in real-time.

The procedure takes place on a specialized chip called a SMRT Cell, which contains millions of microscopic wells known as Zero-Mode Waveguides (ZMWs). A single DNA polymerase enzyme is anchored at the bottom of each ZMW. As the enzyme incorporates fluorescently labeled nucleotides into a complementary DNA strand, the chip detects the light pulses emitted, “reading” the sequence as it is being built.

#### HiFi Reads and Circular Consensus Sequencing (CCS)

Historically, long-read sequencing was known for having higher error rates. However, PacBio solved this using Circular Consensus Sequencing (CCS). In this protocol, the DNA fragment is made circular by attaching “hairpin” adapters. The polymerase can then move around the circle multiple times, reading the same molecule over and over. By comparing these multiple passes, the software can correct random errors, producing HiFi (High Fidelity) reads that are both very long and highly accurate (>99.9%).

#### Advantages and Applications

This technology is a type of long-read sequencing. While Illumina reads are measured in hundreds of base pairs, PacBio HiFi reads are typically 10,000 to 25,000 base pairs long (and can reach much higher in “Continuous Long Read” mode).

Because these sequences are so long, they can easily bridge complex repetitive regions, resolve structural variants, and sequence entire transcripts (Iso-Seq) from start to finish without needing to break them apart. This makes PacBio the gold standard for *de-novo* genome assembly—building a genome from scratch without a reference.

**In summary, PacBio:**

* Utilizes Single Molecule, Real-Time (SMRT) sequencing.
* Offers long read lengths (typically 10kb - 25kb+).
* HiFi reads provide high accuracy (>99.9%), rivaling Illumina.
* Excellent for *de-novo* assembly, structural variant detection, and full-length isoform sequencing.
* Can detect epigenetic modifications (like DNA methylation) directly during the sequencing process.
* Higher cost per megabase compared to Illumina.
* Requires higher amounts of high-molecular-weight (HMW) input DNA.

### 3.1.4 Comparison of Sequencing Technologies

Choosing a sequencing technology involves balancing **read length**, **accuracy**, and **cost**. While Illumina remains the “gold standard” for high-throughput counting and small variant detection, Long-Read technologies (PacBio and Nanopore) are essential for resolving complex genomic structures and repetitive regions.

| Feature | Illumina | PacBio (HiFi) | Oxford Nanopore (ONT) |
| --- | --- | --- | --- |
| **Technology** | Sequencing-by-Synthesis | SMRT (Real-time) | Nanopore (Electrical current) |
| **Read Length** | Short (50–300 bp) | Long (10–25 kb) | Ultra-Long (up to 2 Mb) |
| **Accuracy** | Very High (>99.9%) | Very High (>99.9%) | High (~99%) |
| **Throughput** | Extremely High | Moderate to High | High |
| **Primary Error Type** | Substitutions | Random Indels | Systematic Indels |
| **Best Use Case** | Population SNPs, RNA-seq | *De-novo* assembly, Iso-Seq | Field work, Structural Variants |
| **Portability** | Benchtop/Production | Benchtop/Production | Handheld (MinION) to Production |

### 3.1.5 Summary: Which technology should I choose?

Selecting the right platform depends entirely on your biological question:

* **Illumina:** Best for projects requiring high sample multiplexing, such as differential gene expression (RNA-seq), ChIP-seq, or large-scale population studies (GWAS).
* **PacBio HiFi:** Ideal for “platinum-grade” *de-novo* genome assemblies, resolving full-length transcript isoforms, and high-accuracy structural variant discovery.
* **Oxford Nanopore:** Best for ultra-long read requirements, rapid real-time pathogen identification, and sequencing in remote or field-based environments.

## 3.2 Summary

TipKey Points

* High-throughput sequencing technologies, often called next-generation sequencing (NGS), enable rapid and cost-effective genome sequencing.
* Prominent NGS platforms include Illumina, Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio).
* Each platform employs distinct mechanisms for DNA sequencing, leading to variations in read length, error rates, and applications.
* Illumina sequencing:
  + Uses sequencing-by-synthesis chemistry, produces short read lengths and has high accuracy with low error rates (<1%).
  + While it is scalable and cost-effective for large-scale projects, it is expensive to set up and limited in sequencing long DNA fragments.
* Nanopore sequencing:
  + Uses nanopore technology, provides long read lengths, making it ideal for applications such as *de novo* genome assembly.
  + Although the costs can be higher per base, it is cheaper to set up.
  + Exhibits higher error rates (around 5%), but with significant improvements in recent versions (1%).
* PacBio:
  + Utilizes SMRT (Single Molecule, Real-Time) sequencing-by-synthesis.
  + Offers long read lengths (typically 10,000 to 25,000+ bp).
  + Features high accuracy with HiFi reads (>99.9%).
  + Excellent for **de novo** genome assembly, resolving structural variants, and full-length transcript sequencing (Iso-Seq).
  + Can detect DNA methylation and other base modifications directly.
  + Higher cost per base compared to Illumina.
  + Requires high-quality, high-molecular-weight (HMW) input DNA.