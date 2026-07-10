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

### R-Only Code Example: Dynamic Configuration & Samplesheet Creation
This R script mimics the steps necessary to dynamically build a Nextflow samplesheet and configuration file directly from R, which is a highly recommended automation practice.

```r
# Day 3: Software Pipelines in R
# Automating the creation of a Nextflow samplesheet and a custom config file

# 1. Dynamically Create a Samplesheet for nf-core pipelines
generate_nextflow_samplesheet <- function(reads_directory, output_csv) {
  # Find all fastq.gz files in the target directory
  # Simulated file list
  fastq_files <- c("Sample1_R1.fastq.gz", "Sample1_R2.fastq.gz", 
                   "Sample2_R1.fastq.gz", "Sample2_R2.fastq.gz")
  
  # Extract unique sample names
  sample_names <- unique(gsub("_R.fastq.gz", "", fastq_files))
  
  # Initialize empty data frame conforming to nf-core standard
  samplesheet <- data.frame(
    sample = character(),
    fastq_1 = character(),
    fastq_2 = character(),
    strandedness = character(),
    stringsAsFactors = FALSE
  )
  
  # Populate the dataframe
  for (s in sample_names) {
    fq1 <- file.path(reads_directory, paste0(s, "_R1.fastq.gz"))
    fq2 <- file.path(reads_directory, paste0(s, "_R2.fastq.gz"))
    
    samplesheet <- rbind(samplesheet, data.frame(
      sample = s,
      fastq_1 = fq1,
      fastq_2 = fq2,
      strandedness = "auto"
    ))
  }
  
  # Write to CSV without row names
  write.csv(samplesheet, file = output_csv, row.names = FALSE, quote = FALSE)
  message("Samplesheet successfully written to ", output_csv)
}

# Execute function to build samplesheet
generate_nextflow_samplesheet(reads_directory = "/data/reads", output_csv = "samplesheet.csv")

# 2. Dynamically Generate Advanced Nextflow Configuration File
generate_nf_config <- function(config_file, max_cpus = 8, max_mem = "32.GB") {
  config_content <- paste0(
    "// Advanced Nextflow Configuration generated via R\n",
    "process {\n",
    "  executor = 'slurm'\n",
    "  queue = 'normal'\n",
    "  withLabel: process_high_memory {\n",
    "    queue = 'highmem'\n",
    "    memory = '", max_mem, "'\n",
    "  }\n",
    "}\n",
    "singularity {\n",
    "  enabled = true\n",
    "  autoMounts = true\n",
    "  cacheDir = '/data/.singularity_cache'\n",
    "}\n",
    "executor {\n",
    "  queueSize = 50\n",
    "  submitRateLimit = '10 sec'\n",
    "}\n"
  )
  
  writeLines(config_content, con = config_file)
  message("Nextflow config written to ", config_file)
}

# Execute function
generate_nf_config("custom_hpc.config", max_cpus = 16, max_mem = "64.GB")
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

### R-Only Code Example: Interacting with SLURM arrays and Dependencies
This R script demonstrates how to generate a SLURM job array script and parse parameters directly inside R, simulating a robust job orchestration mechanism commonly used when bridging statistical programming (R) with HPC administration (SLURM).

```r
# Day 4 & 5: HPC Management using R
# Generating SLURM array submission scripts and handling dependencies

# 1. Creating a parameter grid to run across an HPC Job Array
params_df <- expand.grid(
  learning_rate = c(0.01, 0.05, 0.1),
  iterations = c(1000, 5000)
)
params_df$task_id <- 1:nrow(params_df)

# Write parameters to a CSV file to be read by the array tasks
write.csv(params_df, "simulation_params.csv", row.names = FALSE)

# 2. Generate a SLURM array script from R
generate_slurm_array_script <- function(script_name, total_tasks) {
  slurm_header <- c(
    "#!/bin/bash",
    paste0("#SBATCH --job-name=R_Simulation_Array"),
    paste0("#SBATCH -o logs/sim_array_%a.out"),
    paste0("#SBATCH -e logs/sim_array_%a.err"),
    paste0("#SBATCH --time=02:00:00"),
    paste0("#SBATCH -p normal"),
    paste0("#SBATCH -c 2"),
    paste0("#SBATCH --mem=8G"),
    paste0("#SBATCH -a 1-", total_tasks),
    "",
    "module load R/4.2.0",
    "echo \"Starting array task ID: $SLURM_ARRAY_TASK_ID\"",
    "",
    "# In practice, Rscript would parse the SLURM_ARRAY_TASK_ID directly",
    "Rscript run_model.R $SLURM_ARRAY_TASK_ID simulation_params.csv"
  )
  
  writeLines(slurm_header, con = script_name)
  message("SLURM array script written to: ", script_name)
}

generate_slurm_array_script("submit_array.sh", total_tasks = nrow(params_df))

# 3. Simulate R parsing the SLURM_ARRAY_TASK_ID inside the compute node
# (This would be inside the 'run_model.R' script)
run_model_on_hpc <- function() {
  # Fetch task ID from environment (Simulated here as "3")
  task_id_str <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  
  if (task_id_str == "") {
    warning("Not running under SLURM array, defaulting to task 1")
    task_id_str <- "1"
  }
  
  task_id <- as.integer(task_id_str)
  
  # Load parameters
  all_params <- read.csv("simulation_params.csv")
  
  # Isolate specific parameters for THIS task
  my_params <- all_params[all_params$task_id == task_id, ]
  
  message("Running model with learning_rate: ", my_params$learning_rate, 
          " and iterations: ", my_params$iterations)
  
  # Dummy model execution
  Sys.sleep(2)
  message("Task ", task_id, " completed successfully.")
}

# Run the model logic
run_model_on_hpc()

# 4. Automating Job submission and dependency chaining via R
# We can use R to submit the array and hold a dependent job
submit_with_dependency <- function(array_script, dependent_script) {
  # Submit array and capture JOBID (using --parsable to only return the ID)
  submit_cmd <- paste("sbatch", "--parsable", array_script)
  
  # Execute the system command (simulated here)
  # job_id <- system(submit_cmd, intern = TRUE)
  job_id <- "99999" # Dummy ID for documentation
  
  # Chain the next script
  dep_cmd <- sprintf("sbatch --dependency=afterok:%s %s", job_id, dependent_script)
  message("Executing: ", dep_cmd)
  
  # system(dep_cmd)
}
```
# URLs
- https://cambiotraining.github.io/intro-ngs/setup.html
- https://cambiotraining.github.io/hpc-intro/
