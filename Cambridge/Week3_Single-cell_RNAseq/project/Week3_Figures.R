# =============================================================================
# Week 3 project - two figures for the presentation
#
# THE STORY
#   Premise  A mouse gets pregnant, nurses, weans. The gland builds a whole
#            milk-producing tissue and then takes it apart again. Does it go
#            back to how it was?
#   Slide 1  Csn2 (a milk protein) across the four stages. Off, on, on, off.
#            It LOOKS like the gland resets completely.
#   Slide 2  Measure the same gene in the luminal progenitor cells instead of
#            looking at it. After weaning it is NOT back to zero. The
#            progenitors keep the milk programme partly switched on.
#
# Slide 2 contradicts slide 1. That is the talk.
#
# CONFIRMED AGAINST Bach et al. 2017 (doi:10.1038/s41467-017-02001-5):
#   - 25,806 barcodes: 4376 NP, 6021 G, 9603 L, 5806 PI (Methods)
#   - "the post-parous luminal compartment differs from its nulliparous
#      counterpart ... luminal progenitor cells maintain memory of having
#      undergone gestation and lactation" (Results / Discussion)
#   - The genes in section 6 are the ones they use in Fig 5c
#
# Workflow follows course scripts 04-10. Figure 1 is script 10's
# FeaturePlot(split.by=) idea; figure 2 is script 10's per-pseudosample plot.
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

# tidyverse LAST. Bioconductor packages mask count(), filter(), rename().
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(patchwork)
library(tidyverse)

theme_set(theme_classic())
set.seed(42)

candidate_dirs <- c(
  "~/Desktop/Course_Materials/Week3_Single-cell_RNAseq/project_single_cell/week3_project_single_cell",
  "C:/Users/Raiyan Subedar/Downloads/Week 3/week3_project_single_cell"
)
project_dir <- candidate_dirs[dir.exists(path.expand(candidate_dirs))][1]
if (is.na(project_dir)) {
  stop("No candidate project directory exists here:\n  ",
       paste(candidate_dirs, collapse = "\n  "))
}
project_dir <- path.expand(project_dir)
message("project_dir: ", project_dir)

data_dir <- file.path(project_dir, "preprocessed", "cellranger")
out_dir  <- file.path(project_dir, "results")
pres_dir <- file.path(out_dir, "presentation")
dir.create(pres_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(data_dir))

N_PCS        <- 20
CLUSTER_RES  <- 0.4
MT_MAX       <- 5
stage_levels <- c("Nulliparous", "Gestation", "Lactation", "Post-involution")
# Plain English first so the audience follows without a key, the paper's own
# abbreviation in brackets so anyone who knows the dataset can map them.
# Swap "NP"/"G"/"L"/"PI" for the full words if you prefer - they are longer and
# will crowd the panel titles.
stage_plain  <- c("Nulliparous"     = "never pregnant (NP)",
                  "Gestation"       = "pregnant (G)",
                  "Lactation"       = "nursing (L)",
                  "Post-involution" = "after weaning (PI)")

# Milk genes the paper uses in Fig 5c to show the parity memory effect
memory_genes <- c("Csn2", "Csn1s1", "Lalba", "Lipa")

sct_cache <- file.path(out_dir, "figures_sct.rds")


# =============================================================================
# 1. LOAD, QC, NORMALISE   (cached - runs once, then reloads)
#
# If results/02_sct.rds already exists from the full analysis script and used
# the same QC settings, copy it to figures_sct.rds and skip ~20 minutes:
#   file.copy(file.path(out_dir, "02_sct.rds"), sct_cache)
# =============================================================================

if (file.exists(sct_cache)) {

  message("Loading cached object: ", sct_cache)
  seurat_obj <- readRDS(sct_cache)

} else {

  sampleinfo <- read_csv(file.path(project_dir, "sample_info.csv"),
                         show_col_types = FALSE) %>%
    mutate(stage = factor(stage, levels = stage_levels)) %>%
    arrange(stage, sample)

  # ReadMtx, not Read10X: CellRanger v2 output names the gene file
  # genes.tsv.gz, and Read10X looks for features.tsv.gz when it sees gzipped
  # barcodes, so it fails on this dataset.
  read_one <- function(gsm, sample_name) {
    d <- file.path(data_dir, gsm, "outs", "filtered_feature_bc_matrix")
    m <- ReadMtx(mtx      = file.path(d, "matrix.mtx.gz"),
                 features = file.path(d, "genes.tsv.gz"),
                 cells    = file.path(d, "barcodes.tsv.gz"),
                 feature.column = 2)
    colnames(m) <- paste(sample_name, colnames(m), sep = "_")
    m
  }

  mat_list <- map2(sampleinfo$sra_run, sampleinfo$sample, read_one)
  stopifnot(all(map_lgl(mat_list, ~ identical(rownames(.x),
                                              rownames(mat_list[[1]])))))
  counts <- do.call(cbind, mat_list)

  seurat_obj <- CreateSeuratObject(counts, project = "MammaryGland",
                                   min.cells = 1)

  meta <- seurat_obj[[]] %>%
    rownames_to_column("Cell") %>%
    mutate(SampleName = str_remove(Cell, "_[ACGT]+-1$")) %>%
    left_join(sampleinfo, by = c("SampleName" = "sample")) %>%
    mutate(SampleGroup = factor(stage, levels = stage_levels)) %>%
    select(-stage, -sra_run) %>%
    column_to_rownames("Cell")
  seurat_obj[[]] <- meta[colnames(seurat_obj), ]

  # CONFIRMATION AGAINST THE PAPER. Bach et al. Methods report 25,806 barcodes:
  # 4376 NP, 6021 G, 9603 L, 5806 PI. If this does not match, the matrices did
  # not load correctly and nothing downstream is worth running.
  paper_counts <- c(Nulliparous = 4376, Gestation = 6021,
                    Lactation = 9603, `Post-involution` = 5806)
  observed <- table(seurat_obj$SampleGroup)
  print(rbind(observed = as.integer(observed[names(paper_counts)]),
              paper    = paper_counts))
  stopifnot(all(as.integer(observed[names(paper_counts)]) == paper_counts))
  message("Cell counts match the paper exactly.")

  # Mouse mitochondrial prefix is lowercase mt-, not the human MT-
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  stopifnot(sum(seurat_obj$percent.mt) > 0)

  keep_cells <- seurat_obj[[]] %>%
    rownames_to_column("Cell") %>%
    group_by(SampleName) %>%
    filter(nFeature_RNA > (median(nFeature_RNA) - 2 * mad(nFeature_RNA)),
           nCount_RNA   > (median(nCount_RNA)   - 2 * mad(nCount_RNA)),
           percent.mt   < MT_MAX) %>%
    ungroup() %>%
    pull(Cell)
  message("Cells: ", ncol(seurat_obj), " -> ", length(keep_cells), " after QC")
  seurat_obj <- subset(seurat_obj, cells = keep_cells)

  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$SampleName)
  seurat_obj <- SCTransform(seurat_obj, assay = "RNA",
                            vars.to.regress = "percent.mt", verbose = FALSE)

  saveRDS(seurat_obj, sct_cache)
}


# =============================================================================
# 2. CLUSTERING
#
# No integration: each stage is its own 10x run, so "batch" and "stage" are the
# same variable and correcting for one removes the other.
# =============================================================================

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),
                     verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:N_PCS,
                      verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", k.param = 20,
                            dims = 1:N_PCS, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = CLUSTER_RES,
                           algorithm = 4, random.seed = 123,
                           cluster.name = "cluster", verbose = FALSE)
Idents(seurat_obj) <- "cluster"
message("Clusters: ", length(unique(seurat_obj$cluster)))

# Plain log-normalised RNA is what both figures plot. SCT residuals are for
# clustering, not for showing someone an expression level.
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Plain-English stage labels for the figure panels
seurat_obj$Stage <- factor(stage_plain[as.character(seurat_obj$SampleGroup)],
                           levels = unname(stage_plain))


# =============================================================================
# 3. FIND THE LUMINAL PROGENITORS   <<< THE ONLY THING YOU MUST LABEL <<<
#
# You do not need to annotate every cluster. Slide 2 needs one thing: which
# clusters are luminal progenitors. Everything else can stay a number.
#
# Identity is called from mean expression of canonical genes, NOT from
# FindAllMarkers. Milk transcripts leak into the ambient RNA of every droplet
# in a lactating sample, so they rank poorly as "markers" even in cells that
# are full of them - which is how an alveolar cluster ends up reported with
# myoepithelial marker genes.
# =============================================================================

ident_panel <- list(
  progenitor  = c("Aldh1a3", "Kit", "Elf5", "Cd14", "Plet1"),
  milk        = c("Csn2", "Wap", "Lalba", "Glycam1"),
  hormone     = c("Prlr", "Esr1", "Cited1", "Areg"),
  basal       = c("Krt5", "Krt14", "Acta2", "Myh11", "Oxtr"),
  immune      = c("Ptprc", "Lyz2"),
  fibroblast  = c("Col1a1", "Pdgfra"),
  endothelial = c("Cdh5", "Pecam1")
)
ident_panel <- map(ident_panel, ~ intersect(.x, rownames(seurat_obj)))
stopifnot(all(lengths(ident_panel) > 0))

p_ident <- DotPlot(seurat_obj, features = ident_panel, group.by = "cluster") +
  scale_colour_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste0("Canonical markers, resolution ", CLUSTER_RES))
p_ident
ggsave(file.path(pres_dir, "check_identity_dotplot.png"), p_ident,
       width = 14, height = 6, dpi = 150)

# Averages computed directly rather than with AverageExpression(): that
# function passes group names through make.names(), so cluster "1" comes back
# as "X1" and every downstream join by cluster name silently returns NA.
expr_mat    <- LayerData(seurat_obj, assay = "RNA", layer = "data")
panel_genes <- unique(unlist(ident_panel))
cl          <- as.character(seurat_obj$cluster)
cl_levels   <- sort(unique(cl))

avg <- vapply(cl_levels,
              function(k) Matrix::rowMeans(expr_mat[panel_genes, cl == k,
                                                    drop = FALSE]),
              numeric(length(panel_genes)))
avg <- as.matrix(avg)
rownames(avg) <- panel_genes
colnames(avg) <- cl_levels

# z-score each gene across clusters, leaving flat genes at 0 instead of NaN
gene_sd <- apply(avg, 1, sd)
avg_z   <- avg
avg_z[gene_sd > 0, ] <- t(scale(t(avg[gene_sd > 0, , drop = FALSE])))
avg_z[gene_sd == 0, ] <- 0
stopifnot(!any(is.na(avg_z)))

ident_summary <- map_dfr(names(ident_panel), function(k) {
  g <- intersect(ident_panel[[k]], rownames(avg_z))
  tibble(compartment = k, cluster = colnames(avg_z),
         score = colMeans(avg_z[g, , drop = FALSE]))
}) %>%
  pivot_wider(names_from = compartment, values_from = score) %>%
  mutate(across(-cluster, ~ round(.x, 2)))

# Slide 2 compares never-pregnant against after-weaning, so a usable progenitor
# cluster needs cells at BOTH. This table puts identity and availability
# side by side, sorted so the candidates are at the top.
cluster_context <- seurat_obj[[]] %>%
  group_by(cluster) %>%
  summarise(n = n(),
            n_NP = sum(SampleGroup == "Nulliparous"),
            n_PI = sum(SampleGroup == "Post-involution"),
            .groups = "drop") %>%
  mutate(cluster = as.character(cluster))

ident_summary <- left_join(cluster_context, ident_summary, by = "cluster") %>%
  arrange(desc(progenitor))
print(ident_summary, n = Inf)
write_csv(ident_summary, file.path(pres_dir, "check_identity_table.csv"))

# A cluster scoring high for both milk and basal is two cell types merged.
# Raise CLUSTER_RES and re-run section 2 if this fires.
merged <- ident_summary %>% filter(milk > 0.5, basal > 0.5)
if (nrow(merged) > 0) {
  print(merged)
  warning("Merged cluster(s) above - raise CLUSTER_RES and re-run section 2.")
}


# =============================================================================
# 4. YOUR CHOICE   <<< FILL THIS IN FROM SECTION 3 <<<
#
# Pick the cluster(s) with a high progenitor score AND enough cells in both
# n_NP and n_PI. Usually one cluster; occasionally two.
# =============================================================================

progenitor_clusters <- c()     # e.g. c("2") or c("2", "9")

stopifnot(length(progenitor_clusters) > 0)
stopifnot(all(progenitor_clusters %in% as.character(seurat_obj$cluster)))

seurat_obj$IsProgenitor <- as.character(seurat_obj$cluster) %in% progenitor_clusters
table(seurat_obj$IsProgenitor, seurat_obj$SampleGroup)

# --- Compartment labels for the reference panel in figure 1 -----------------
# Read off the section 3 dot plot. These are ORIENTATION ONLY - they tell the
# audience what part of the map they are looking at. The result on slide 2
# depends on progenitor_clusters above, not on these, so a debatable call
# between, say, milk-making and basal does not change any claim you make.
#
# Anything uncertain goes to "other" rather than being forced into a category.
compartment_map <- c(
  "1"  = "basal / myoepithelial",  # Krt5, Krt14, Acta2, Myh11, Oxtr
  "2"  = "luminal progenitor",     # Aldh1a3, Kit, Elf5, Cd14, Plet1
  "3"  = "hormone-sensing",        # Prlr, Esr1, Cited1, Areg
  "4"  = "hormone-sensing",        # Prlr, Esr1, Cited1, Areg
  "5"  = "milk-making",            # Csn2, Glycam1 + progenitor markers
  "6"  = "milk-making",            # Csn2 + Kit, Elf5, Cd14
  "7"  = "basal / myoepithelial",
  "8"  = "basal / myoepithelial",  # Oxtr highest here
  "9"  = "basal / myoepithelial",
  "10" = "other",                  # Ptprc, Lyz2 - immune
  "11" = "other",                  # Col1a1, Pdgfra - fibroblast
  "12" = "other",                  # Cdh5, Pecam1 - endothelial
  "13" = "basal / myoepithelial",
  "14" = "other",                  # milk + progenitor + basal all high: unclear
  "15" = "hormone-sensing"
)

compartment_levels <- c("luminal progenitor", "hormone-sensing",
                        "milk-making", "basal / myoepithelial", "other")
stopifnot(setequal(names(compartment_map),
                   as.character(unique(seurat_obj$cluster))))
stopifnot(all(compartment_map %in% compartment_levels))

seurat_obj$Compartment <- factor(
  unname(compartment_map[as.character(seurat_obj$cluster)]),
  levels = compartment_levels)
table(seurat_obj$Compartment, seurat_obj$SampleGroup)


# =============================================================================
# 5. GO / NO-GO   <<< CHECK THIS BEFORE BUILDING ANY SLIDES <<<
#
# The whole talk rests on one claim: in luminal progenitors, milk genes are
# higher after weaning than in never-pregnant mice. One value per mouse, so
# the two replicates are visible rather than hidden behind thousands of cells.
#
# PROCEED only if the after-weaning mice sit above the never-pregnant mice for
# most genes. If they do not, do not force it - pick different genes from the
# paper's Fig 5c list (Xdh, Cd36, Cidea) or change the story.
# =============================================================================

memory_genes <- intersect(memory_genes, rownames(seurat_obj))
stopifnot(length(memory_genes) > 0)

prog_cells <- colnames(seurat_obj)[
  seurat_obj$IsProgenitor &
    seurat_obj$SampleGroup %in% c("Nulliparous", "Post-involution")]
message("Progenitor cells used: ", length(prog_cells))

prog_expr <- FetchData(seurat_obj, cells = prog_cells, layer = "data",
                       vars = c(memory_genes, "SampleGroup", "SampleName")) %>%
  pivot_longer(all_of(memory_genes), names_to = "gene", values_to = "expr") %>%
  mutate(gene = factor(gene, levels = memory_genes))

per_mouse <- prog_expr %>%
  group_by(gene, SampleGroup, SampleName) %>%
  summarise(mean_expr = mean(expr),
            pct_cells = 100 * mean(expr > 0),
            n_cells   = n(), .groups = "drop")
print(per_mouse, n = Inf)
write_csv(per_mouse, file.path(pres_dir, "fig2_per_mouse.csv"))

# --- SPECIFICITY CONTROL ----------------------------------------------------
# The obvious objection to slide 2 is ambient RNA: milk transcripts float loose
# in the suspension and get packaged into every droplet, so maybe the
# progenitors are not expressing these genes at all, just sitting in a soup.
#
# If that were true the signal would be spread evenly across all cell types in
# the same sample. So compare progenitors against every other cell in the SAME
# after-weaning mice. If the progenitors are clearly higher, it is expression,
# not contamination.
pi_cells <- colnames(seurat_obj)[seurat_obj$SampleGroup == "Post-involution"]

specificity <- FetchData(seurat_obj, cells = pi_cells, layer = "data",
                         vars = c(memory_genes, "SampleName", "IsProgenitor")) %>%
  pivot_longer(all_of(memory_genes), names_to = "gene", values_to = "expr") %>%
  mutate(cell_group = ifelse(IsProgenitor, "progenitors", "all other cells")) %>%
  group_by(gene, SampleName, cell_group) %>%
  summarise(mean_expr = round(mean(expr), 3),
            pct_cells = round(100 * mean(expr > 0), 1),
            n_cells   = n(), .groups = "drop") %>%
  arrange(gene, SampleName, desc(cell_group))
print(specificity, n = Inf)
write_csv(specificity, file.path(pres_dir, "fig2_specificity_control.csv"))


# =============================================================================
# 6. FIGURE 1 - the milk gene across the cycle
#
# Four UMAP panels, one per stage, coloured by Csn2. Off, on, on, off.
# No cell-type labels involved, so nothing here depends on section 4.
# =============================================================================

# Equal cells per panel. Nursing has ~9600 cells and never-pregnant ~4400, so
# without this the nursing panel looks denser whatever the expression does.
set.seed(42)
n_per_stage <- min(table(seurat_obj$Stage))
balanced_cells <- seurat_obj[[]] %>%
  rownames_to_column("Cell") %>%
  group_by(Stage) %>%
  slice_sample(n = n_per_stage) %>%
  ungroup() %>%
  pull(Cell)
message("Cells per panel: ", n_per_stage)

obj_bal <- subset(seurat_obj, cells = balanced_cells)

# keep.scale = "all" is essential: without it each panel is rescaled to its own
# maximum and the "off" panels look as bright as the "on" ones.
# max.cutoff = "q99" matters as much as keep.scale. A handful of cells carry
# enormous Csn2 counts; without a ceiling they set the top of the colour scale
# and every ordinary expressing cell washes out to pale pink.
# coord_fixed() stops the panels being stretched vertically.
# Tight margins: with coord_fixed() every panel keeps a square aspect, so any
# spare room in its slot becomes a gap between panels. Small margins plus a
# canvas sized to the panel count keeps them close together.
blank_axes <- theme(axis.title  = element_blank(),
                    axis.text   = element_blank(),
                    axis.ticks  = element_blank(),
                    axis.line   = element_blank(),
                    plot.title  = element_text(size = 15, hjust = 0.5),
                    plot.margin = margin(2, 2, 2, 2))

# Reference panel: what part of the map is what. Without this the audience has
# to take on faith that the patch still lit up after weaning is the progenitors.
comp_cols <- c("luminal progenitor"    = "#E7298A",
               "hormone-sensing"       = "#7570B3",
               "milk-making"           = "#1B9E77",
               "basal / myoepithelial" = "#66A61E",
               "other"                 = "grey75")

# Legend to the RIGHT of this panel, not underneath it. Underneath, it creates
# an empty band running the full width of the figure and pushes every panel
# smaller, which is what left the big gaps.
p_ref <- DimPlot(obj_bal, reduction = "umap", group.by = "Compartment",
                 pt.size = 0.4) +
  scale_colour_manual(values = comp_cols) +
  coord_fixed() +
  ggtitle("cell types") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3.5))) +
  blank_axes +
  theme(legend.position = "right", legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key.height = unit(14, "pt"))

# max.cutoff = "q99" matters as much as keep.scale. A handful of cells carry
# enormous Csn2 counts; without a ceiling they set the top of the colour scale
# and every ordinary expressing cell washes out to pale pink.
# coord_fixed() stops the panels being stretched vertically.
p_csn2 <- FeaturePlot(obj_bal,
                      features   = "Csn2",
                      split.by   = "Stage",
                      reduction  = "umap",
                      keep.scale = "all",
                      order      = TRUE,
                      max.cutoff = "q99",
                      cols       = c("grey88", "#B22222"),
                      pt.size    = 0.4) &
  coord_fixed() &
  blank_axes &
  # Say what the colour bar is. Unlabelled numbers on a legend invite the
  # question "what is that scale?" during questions; this answers it on screen.
  labs(colour = "Csn2\nexpression\n(log-normalised\ncounts)") &
  theme(legend.position  = "right",
        legend.title     = element_text(size = 11, lineheight = 1.1),
        legend.text      = element_text(size = 10))
# --- Layout -----------------------------------------------------------------
# A single row of five panels is too wide for a 16:9 slide. This puts the
# reference map on the left and the four stages in a 2x2 block:
#
#     A A B C        A = cell types (spans 2x2, so coord_fixed keeps it square)
#     A A D E        B..E = the four stages
#
# FeaturePlot(split.by=) returns a patchwork, so its panels can be pulled out
# individually with [[ ]] while keeping the shared colour scale from
# keep.scale = "all".
#
# Alternatives if you want different proportions:
#   "ABC\nADE"  reference narrower, stages larger
#   "AB\nCD"    drop the reference entirely, stages only (use p_csn2[[1]]..[[4]])
fig1 <- p_ref + p_csn2[[1]] + p_csn2[[2]] + p_csn2[[3]] + p_csn2[[4]] +
  plot_layout(design = "AABC\nAADE", guides = "collect")
fig1

# coord_fixed() keeps every panel square, so height sets the panel size: each
# row is height/2, the reference is one full height, and each stage panel is
# half of it. Width then has to fit reference + 2 stages + both legends.
# If gaps open up between panels, reduce width - do not increase height.
ggsave(file.path(pres_dir, "fig1_csn2_across_cycle.png"), fig1,
       width = 16, height = 7, dpi = 300)


# =============================================================================
# 7. FIGURE 2 - the same gene, measured in the progenitors
#
# Violins are the cells; the black points are the two mice. Course script 10
# ends the same way - one point per pseudosample, so n = 2 stays visible.
# =============================================================================

fig2 <- ggplot(prog_expr, aes(x = SampleGroup, y = expr, fill = SampleGroup)) +
  geom_violin(scale = "width", colour = NA) +
  # small fixed jitter so the two mice do not land on top of each other
  geom_point(data = per_mouse, aes(y = mean_expr),
             colour = "black", size = 3,
             position = position_jitter(width = 0.07, height = 0, seed = 1)) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Nulliparous"     = "grey75",
                               "Post-involution" = "#B22222")) +
  # three short lines, not two long ones - at this facet width the labels
  # collide otherwise
  scale_x_discrete(labels = c("Nulliparous"     = "never\npregnant\n(NP)",
                              "Post-involution" = "after\nweaning\n(PI)")) +
  labs(x = NULL, y = "expression (log-normalised)",
       title = "Luminal progenitor cells only",
       subtitle = "violin = cells, black points = the two mice") +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.x = element_text(size = 11, lineheight = 0.95),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 15))
fig2
ggsave(file.path(pres_dir, "fig2_progenitor_memory.png"), fig2,
       width = 11, height = 5, dpi = 300)

saveRDS(seurat_obj, file.path(out_dir, "figures_final.rds"))
sessionInfo()
