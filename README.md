# KAUST Academy Bioinformatics — Stage 4 (Summer)

Coursework, practicals, and mini-projects from the KAUST Academy Bioinformatics Stage 4 summer programme, run in partnership with the University of Cambridge. The programme runs in two phases: a two-week **Pre-Cambridge** foundations block (R, statistics, NGS, pipelines, HPC), followed by a four-week **Cambridge** block of specialised modules, each ending in a self-directed mini-project.

Every module folder contains a `Summary.md` — a concise write-up of that module in a consistent format: a short overview, one section per teaching day (prose + a representative code block), and a final section on the mini-project with its key figures.

## Repository structure

```
├── Pre-Cambridge/        Foundations (2 weeks): R, stats, NGS, pipelines, HPC
└── Cambridge/            Specialised modules (4 weeks), each with course/ + project/
    ├── Week1_Epidemiological_Modelling/
    ├── Week2_Expression_Proteomics/
    ├── Week3_Single-cell_RNAseq/
    └── Week4_Metagenomics/
```

The end-of-programme final project lives in a separate repository: [KASP-Cambridge_microbiome-meta-reanalysis](https://github.com/Raiyan-S/KASP-Cambridge_microbiome-meta-reanalysis).

Within each Cambridge week, `course/` holds the taught materials and the module `Summary.md`, and `project/` holds the mini-project data, scripts, and slides.

## Pre-Cambridge — Foundations

| Week | Topic | Summary |
|------|-------|---------|
| One | Data visualisation and statistical inference in R | [Summary](Pre-Cambridge/Week%20One/Summary.md) |
| Two | Bioinformatics, data visualisation, NGS, pipelines & HPC | [Summary](Pre-Cambridge/Week%20Two/Summary.md) |

Week One covers R programming, `tidyverse` data wrangling, `ggplot2`, and statistical inference (t-tests, ANOVA, regression, model selection, power analysis). Week Two moves into applied bioinformatics: data visualisation principles, NGS data processing (QC, trimming, mapping, variant calling), reproducible software and pipelines (Mamba, containers, Nextflow/nf-core), and working on an HPC cluster (SSH, SLURM, job arrays).

## Cambridge — Specialised Modules

Each module is three to five days of instruction plus a self-directed mini-project.

### Week 1 — Epidemiological Modelling
[Course summary](Cambridge/Week1_Epidemiological_Modelling/course_epi/Summary.md)

Compartmental (SIR/SEIR) models of infectious disease dynamics in R, solving ODEs with `deSolve`, extending models with demography and vaccination, comparing models to data, and fitting serocatalytic models. **Mini-project:** a SEIR-with-demography-and-vaccination model comparing measles re-introduction risk in the UK vs Saudi Arabia, framed around the herd-immunity threshold.

### Week 2 — Expression Proteomics
[Course summary](Cambridge/Week2_Expression_Proteomics/course_prot%20(DIA)/Summary.md)

Analysis of DIA mass-spectrometry proteomics in R/Bioconductor: importing DIA-NN output into `QFeatures`, cleaning and filtering, normalisation and aggregation to proteins, `limma` differential abundance, and GO functional enrichment. **Mini-project:** reanalysis of *Chlorella vulgaris* light-acclimation proteomes (Cecchin et al. 2023), showing that antenna, photosystems, and carbon-capture machinery down-regulate together under high light.

### Week 3 — Single-cell RNA-seq
[Course summary](Cambridge/Week3_Single-cell_RNAseq/course/Summary.md)

The Seurat v5 workflow: Cell Ranger, QC, normalisation and feature selection, dimensionality reduction, batch-correction/integration, clustering, marker genes, and differential expression/abundance. **Mini-project:** a reanalysis of the Bach et al. 2017 mouse mammary gland atlas showing that luminal progenitor cells retain a "parity memory" — keeping the milk-gene programme partly active after weaning.

### Week 4 — Metagenomics
[Course summary](Cambridge/Week4_Metagenomics/course/Summary.md)

Amplicon (16S, DADA2/phyloseq) and shotgun metagenomics: QC, k-mer/database profiling (MASH, Kraken2), alignment-based methods (bowtie2), de novo assembly (metaSPAdes), and metagenome-assembled genomes (MaxBin2, CheckM, GTDB-Tk). **Mini-project:** a *C. difficile* infection (CDI) vs faecal-transplant (FMT) gut-microbiome analysis showing that fibre-fermenting anaerobes lost during infection recover after transplant while opportunistic blooms recede.

## Final Project

The end-of-programme final project — a microbiome meta-analysis reanalysis — is maintained in its own repository: [Raiyan-S/KASP-Cambridge_microbiome-meta-reanalysis](https://github.com/Raiyan-S/KASP-Cambridge_microbiome-meta-reanalysis).

## Notes

The `Summary.md` write-ups are study notes: code blocks are condensed from the course practicals and personal project scripts and are meant as readable references rather than turnkey pipelines. Mini-project figures are the author's own, taken from the project scripts and presentation slides. Large data files and reference databases are kept inside the relevant `project/` folders.
