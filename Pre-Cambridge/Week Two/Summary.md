# Bioinformatics, Data Visualization, and HPC Analysis Pipeline.

This document covers Data Visualization, Next-Generation Sequencing (NGS) data processing, Software Management, and High-Performance Computing (HPC) workflows.

## Day 1: Data Visualization

### Overview
Data visualization serves two primary purposes: **exploration** to understand data features and discover patterns, and **communication** to clearly tell a story to an audience. Effective visualization relies on understanding how to construct figures using geometrical primitives called **marks** (e.g., points, lines, shapes) and **channels** that dictate their appearance (e.g., color, length, position, size). 

Choosing the right chart type is essential based on the data:
*   **Distribution:** Histograms, Boxplots, and Violins.
*   **Correlation:** Scatterplots, Heatmaps, and Bubble charts.
*   **Ranking:** Bar charts, Lollipops, and Spider plots.

When mapping data to visual channels, **expressiveness and effectiveness** are key. Quantitative data is best represented by length, position on a scale, or area, while qualitative/categorical data relies on spatial grouping or color hue. Relying purely on color saturation or luminance is considered poor for quantitative data, whereas spatial position is highly effective.

### Example: Data Exploration & Figure Design
This R code demonstrates how to apply the principles of marks and channels (position, color hue, and size) to effectively explore and communicate distributions and correlations using the built-in R tools and `ggplot2` logic.

```r
# Day 1: Data Visualization (Exploration and Communication)
# Exploring distributions and correlations matching marks and channels

# Load necessary library
library(ggplot2)

# 1. Simulate a dataset with categorical and quantitative variables
set.seed(42)
genomic_data <- data.frame(
  sample_id = paste0("Sample_", 1:200),
  expression_level = c(rnorm(100, mean = 15, sd = 3), rnorm(100, mean = 25, sd = 4)),
  mutation_count = rpois(200, lambda = 5) * c(runif(100, 1, 2), runif(100, 2, 4)),
  treatment_group = rep(c("Control", "Treated"), each = 100),
  gene_cluster = sample(c("Cluster_A", "Cluster_B", "Cluster_C"), 200, replace = TRUE)
)

# 2. DISTRIBUTION: Boxplot and Violin Plot
# Mark: Shapes/Lines. Channels: Position (Y-axis for quantitative), Color Hue (Categorical)
distribution_plot <- ggplot(genomic_data, aes(x = treatment_group, y = expression_level, fill = treatment_group)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = 16) +
  scale_fill_manual(values = c("Control" = "#3498db", "Treated" = "#e74c3c")) + # High discriminability
  labs(
    title = "Expression Level Distribution by Treatment",
    subtitle = "Combining Violin and Boxplots to explore data distributions",
    x = "Treatment Group",
    y = "Gene Expression Level"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 3. CORRELATION: Bubble Scatterplot
# Mark: Points. Channels: Position (X/Y), Color Hue (Categorical), Size (Area for quantitative)
correlation_plot <- ggplot(genomic_data, aes(x = expression_level, y = mutation_count, 
                                             color = gene_cluster, size = expression_level)) +
  geom_point(alpha = 0.7) +
  scale_color_brewer(palette = "Set2") + # Perceptually uniform qualitative palette
  labs(
    title = "Correlation between Expression and Mutation Count",
    subtitle = "Using X/Y position, Color Hue, and Area (Bubble size) as channels",
    x = "Expression Level",
    y = "Mutation Count",
    color = "Gene Cluster",
    size = "Expression Level"
  ) +
  theme_classic()

# Display plots
print(distribution_plot)
print(correlation_plot)
```
![Distribution_Plot](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/Day1_Distribution_Plot.png)
![Distribution_Plot](https://github.com/Raiyan-S/KAUST-Academy-Bioinformatics-Stage-4-Summer/blob/main/Pre-Cambridge/Images/Day1_Correlation_Plot.png)


---

## Day 2: NGS Data Processing

### Overview
Next-Generation Sequencing (NGS) allows for reading genetic code at unprecedented speeds, making it pivotal for clinical diagnostics, agriculture, and evolutionary biology. Sequencing technologies vary: **Illumina** is the gold standard for high-accuracy (<1% error rate) short reads, while **PacBio** and **Oxford Nanopore (ONT)** are long-read technologies critical for *de-novo* assembly and structural variants.

An NGS pipeline typically follows these core steps:
1.  **Common File Formats:** Raw data comes as **FASTQ** (sequences + quality scores). References use **FASTA** and annotations use **GFF/GTF**. Alignments are saved in **SAM/BAM** (binary), and the final mutations are stored in **VCF** format. 
2.  **QC & Pre-processing:** "Garbage in, garbage out" is a core principle. **FastQC** evaluates base quality (looking for dephasing and signal decay) and adapter contamination. **fastp** performs quality trimming and adapter removal. Outputs are summarized via **MultiQC**.
3.  **Reference Genomes:** You must download reference FASTA and GFF/GTF files from the *same source and version* (e.g., Ensembl) to ensure chromosome naming matches. The reference must be indexed using tools like `samtools` or mapper-specific indices.
4.  **Sequence Alignment:** **BWA-MEM 2** is a popular aligner that maps reads to the reference. Output SAM files are massive and must be converted to coordinate-sorted, indexed **BAM** files.
5.  **Variant Calling:** Tools like **BCFtools** calculate the probability of a variant (SNPs or Indels) existing at a specific coordinate based on read depth, base quality, and mapping quality.

### Shell Code Example:
This bash script demonstrates how to automate the shell commands required for QC, trimming, mapping, and variant calling.

```bash
# 1. Investigate common file formats
zcat data/reads/SRR36763845_1.fastq.gz | wc -l # Count lines in a compressed FASTQ
grep "^>" resources/reference/genome.fasta | wc -l # Count sequences in a FASTA

# 2. Pre-processing: Quality Control and Trimming
# Run FastQC on paired-end reads
fastqc data/reads/SRR36763845_1.fastq.gz data/reads/SRR36763845_2.fastq.gz -o results/qc/

# Clean reads with fastp (quality filtering and adapter trimming)
fastp -i data/reads/SRR36763845_1.fastq.gz -I data/reads/SRR36763845_2.fastq.gz \
      -o results/qc/SRR36763845_1_trimmed.fastq.gz -O results/qc/SRR36763845_2_trimmed.fastq.gz \
      --cut_tail --html results/qc/SRR36763845_report.html

# Aggregate reports with MultiQC
multiqc results/qc/ -o results/qc/

# 3. Reference Genomes & Mapping
# Create a simple FASTA index with samtools
samtools faidx resources/reference/genome.fasta

# Index the genome for the aligner
bwa-mem2 index resources/reference/genome.fasta

# Map reads to the reference using 4 CPU cores
bwa-mem2 mem -t 4 resources/reference/genome.fasta \
    results/qc/SRR36763845_1_trimmed.fastq.gz \
    results/qc/SRR36763845_2_trimmed.fastq.gz > results/mapping/alignments.sam

# Convert SAM to a coordinate-sorted BAM and index it
samtools sort results/mapping/alignments.sam -o results/mapping/alignments.sorted.bam
samtools index results/mapping/alignments.sorted.bam

# 4. Variant Calling
# Call variants and output to a VCF file
bcftools mpileup -f resources/reference/genome.fasta results/mapping/alignments.sorted.bam | \
    bcftools call -mv -Ov -o results/variants/variants.vcf

# Filter variants based on Quality (>30) and Depth (>10)
bcftools filter -i 'QUAL>30 && DP>10' results/variants/variants.vcf -o results/variants/variants_filtered.vcf

# View the filtered variants, ignoring the VCF header
grep -v "^#" results/variants/variants_filtered.vcf | wc -l
bcftools view results/variants/variants_filtered.vcf | less -S
```

---

## Day 3: Software & Pipelines

### Overview
Modern bioinformatics demands highly reproducible software management. **Package managers** like **Conda/Mamba** automate the installation of software and their dependencies inside isolated environments, eliminating dependency chaos and version conflicts across projects.

When Mamba environments aren't enough, **Container Virtualization** using **Singularity** (now Apptainer) or Docker is used. Containers package software, the operating system, and libraries into a single file. Singularity is heavily preferred on HPC systems because, unlike Docker, it does not require root/admin permissions to run, drastically improving security.

Managing multi-step pipelines requires **Automated Workflows** like **Nextflow** or Snakemake. 
*   The **nf-core** community maintains standardized, highly-curated pipelines (like `nf-core/rnaseq`).
*   Nextflow manages task execution, parallelization, and caches intermediate files in a `work/` directory, allowing you to `-resume` pipelines that fail without starting over.
*   You can utilize advanced configuration profiles (e.g., `-profile singularity`) to enforce container use. By altering a `nextflow.config` file, users dynamically alter resource allocations using selectors like `withLabel` (e.g., matching a high-memory task to a specific HPC queue).
*   Most nf-core pipelines require a **Samplesheet** (a CSV file) to pass experimental metadata and file paths into the pipeline.

**Shell Code Example:**
```bash
# 1. Package Managers (Mamba)
# Create an isolated environment installing fastqc and multiqc from bio-specific channels
mamba create -n qc -c conda-forge -c bioconda fastqc=0.12.1 multiqc=1.21
mamba activate qc

# Export the environment for reproducibility
mamba env export > env.yaml

# 2. Container Virtualisation (Singularity)
# Download a pre-built bioinformatics container (e.g., SeqKit)
singularity pull https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0

# Run a command securely inside the container, binding a custom filesystem path
singularity run --bind /scratch/robin/awesomeproject:/scratch/robin/awesomeproject \
    seqkit-2.8.2.sif seqkit stats reads/*.fastq.gz

# 3. Pipelines & Automated Workflows (Nextflow)
# Create a sample sheet for nf-core dynamically using bash
echo "sample,fastq_1,fastq_2" > samplesheet.csv
ls reads/*_1.fastq.gz | awk -F"/" '{print $2}' | awk -F"_1" '{print $1",""reads/"$0",reads/"$1"_2.fastq.gz"}' >> samplesheet.csv

# Set up a persistent terminal session for long-running workflows
tmux new -s nf_run

# Create an advanced custom HPC configuration file for Nextflow
cat << 'EOF' > custom_hpc.config
process {
    executor = 'slurm'
    queue = 'normal'
    withLabel: 'process_high_memory' {
        queue = 'highmem'
        memory = 200.GB
    }
}
executor {
    queueSize = 50
}
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/data/participant/.nextflow-singularity-cache'
}
EOF
#

# Launch an nf-core workflow using singularity, resume capabilities, and our custom HPC config
nextflow run nf-core/demo -revision "1.0.0" \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome/genome.fasta \
    -c custom_hpc.config \
    -resume

# Detach from tmux to leave the job running (Ctrl+B, then D)
# Later, reattach to the session:
tmux attach -t nf_run

# Clean up cached Nextflow intermediate files to free up disk space
nextflow clean -force
```

---

## Day 4 & 5: Working on an HPC

### Overview
High-Performance Computing (HPC) clusters consist of **Login nodes** (for editing scripts, moving files, submitting jobs) and **Compute nodes** (the heavy lifters with massive RAM/CPUs that actually run calculations). The filesystem is split into small, backed-up **home directories** (for config/software) and massive, high-speed **scratch spaces** (for data processing).

**File Transfer:** Files are moved onto the HPC via GUI tools (Filezilla) or command-line tools like `scp` (for basic transfers) and `rsync` (for advanced data synchronization).

**SLURM Scheduler:** You interact with compute nodes via a Job Scheduler like **SLURM**. 
*   Jobs are requested using `sbatch` wrapped in shell scripts that start with a shebang (`#!/bin/bash`) and `#SBATCH` directives to request CPUs (`-c`), memory (`--mem`), and partitions (`-p`).
*   Check the queue with `squeue` and track resource efficiency using `seff` or `sacct`.

**Job Parallelization & Dependencies:**
*   **Job Arrays:** Using `#SBATCH -a 1-10` submits 10 independent jobs. SLURM provides the `$SLURM_ARRAY_TASK_ID` variable, allowing a single script to process multiple files in parallel.
*   **Job Dependencies:** Advanced pipelines chain tasks together. Using `--dependency=afterok:JOBID` ensures a job only starts if a previous task finished successfully. Alternatively, `--dependency=afternotok:JOBID` can restart failed checkpoints.

**Shell Code Example:**
```bash
# 1. Remote Work & HPC Introduction
# Connect to the HPC using SSH
ssh username@login.hpc.cam.ac.uk

# Check available resources on the login node (Do not run heavy jobs here!)
free -h
nproc --all

# 2. Software Management via Modules
# Search for and load pre-installed cluster software
module avail bowtie2
module load bowtie2/2.5.0

# 3. Job Parallelisation (Arrays) and SLURM Scheduler
# Create a parallel job array script to parse parameters and map reads dynamically
cat << 'EOF' > parallel_mapping.sh
#!/bin/bash
#SBATCH -J parallel_mapping
#SBATCH -D /home/username/rds/hpc-work/
#SBATCH -o job_logs/mapping_%a.log
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -a 2-9

# Load software environment securely on compute nodes
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate bioinformatics

# Use the array ID to read the corresponding line from a CSV file
# Example: line 2 fetches parameters for the 1st sample 
SAMPLE_INFO=$(head -n $SLURM_ARRAY_TASK_ID data/samplesheet.csv | tail -n 1)
SAMPLE_NAME=$(echo $SAMPLE_INFO | cut -d ',' -f 1)
READ1=$(echo $SAMPLE_INFO | cut -d ',' -f 2)
READ2=$(echo $SAMPLE_INFO | cut -d ',' -f 3)

echo "Running mapping using $SLURM_CPUS_PER_TASK CPUs for $SAMPLE_NAME"
bowtie2 -p $SLURM_CPUS_PER_TASK -x genome_index -1 $READ1 -2 $READ2 -S results/${SAMPLE_NAME}.sam
EOF
#

# 4. Job Dependencies
# Submit the array and capture the Job ID cleanly using --parsable
ARRAY_JOB_ID=$(sbatch --parsable parallel_mapping.sh)

# Create a summary script that strictly depends on the array finishing successfully
cat << 'EOF' > summarize_mapping.sh
#!/bin/bash
#SBATCH -J summary
#SBATCH -o job_logs/summary.log
#SBATCH -c 1

echo "All mapping tasks have completed successfully. Summarizing BAM files."
multiqc results/
EOF
#

# Submit the dependent job
sbatch --dependency=afterok:$ARRAY_JOB_ID summarize_mapping.sh

# 5. Monitoring
# Check the queue for running and pending jobs
squeue -u username

# Check job efficiency (e.g., memory and CPU usage) after a job finishes
seff $ARRAY_JOB_ID

# Cancel jobs if an error was made
scancel $ARRAY_JOB_ID

# 6. File Transfer
# Synchronize output results from the HPC back to the local computer
rsync -auvh --progress username@login.hpc.cam.ac.uk:/home/username/rds/hpc-work/results/ /path/to/local/results/
```
# URLs
- https://cambiotraining.github.io/intro-ngs/setup.html
- https://cambiotraining.github.io/hpc-intro/
