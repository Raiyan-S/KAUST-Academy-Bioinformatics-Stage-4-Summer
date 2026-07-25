## =============================================================
## Lesson 07: Visualisation of differential abundance results
## =============================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("matrixStats")
library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")


## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------------------
## Heatmap of all replicated proteins
## ---------------------------------------------------------------------------------------------------------
quant_mtx <- assay(dia_qf[["norm_proteins_replicated"]])
colnames(quant_mtx)

sample_annotations <- colData(dia_qf) %>%
  data.frame() %>%
  select(group)

head(sample_annotations)



## Extract quantification matrix and rename columns with sample shortnames.
## We use correlation-based distance (1 - Pearson r) with pairwise complete
## observations to handle missing values, and Ward's linkage for clustering.


## Correlation-based distances (handles NAs via pairwise.complete.obs)
# Calculate correlation-based distance manually, handling NAs with pairwise complete obs
dist_cols <- as.dist(1 - cor(quant_mtx, use = "pairwise.complete.obs"))
dist_rows <- as.dist(1 - cor(t(quant_mtx), use = "pairwise.complete.obs"))

pheatmap(
  quant_mtx,
  scale = 'row', # standardise each row to mean=0 and sd=1 (Z-score)
  show_rownames = FALSE,
  annotation_col = sample_annotations,
  clustering_method = 'ward.D2', # Use Ward's method for hierarchical clustering
  clustering_distance_rows = dist_rows,
  clustering_distance_cols = dist_cols,
  color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100), # Use a diverging color palette
  annotation_names_col = FALSE, # Hide annotation title
  treeheight_row = 100, # default is 50 — increase to give row dendrogram more space
  treeheight_col = 100, # default is 50 — increase to give col dendrogram more space
  cellwidth = 10, # fix cell width in points
  cellheight = 0.25
)

## ---------------------------------------------------------------------------------------------------------
## Helper function for repeated heatmaps
## ---------------------------------------------------------------------------------------------------------
colnames(rowData(dia_qf[["norm_proteins_replicated"]]))
rowData(dia_qf[["norm_proteins_replicated"]]["Q96KN2",])

## Rather than repeating the same code for each protein subset, we define a
## function that takes a vector of UniProt IDs and plots the heatmap.
plot_heatmap_for_uids <- function(uids) {
  quant_mtx <- assay(dia_qf[["norm_proteins_replicated"]][uids, ])
  colnames(quant_mtx) <- dia_qf$runCol
  
  # Use concatenation of UniprotID and Gene name as row label
  rownames(quant_mtx) <- sprintf(
    '%s (%s)',
    rowData(dia_qf[["norm_proteins_replicated"]][uids, ])$Genes,
    rownames(quant_mtx)
  )
  # Calculate correlation-based distance manually, handling NAs with pairwise complete obs
  dist_cols <- as.dist(1 - cor(quant_mtx, use = "pairwise.complete.obs"))
  dist_rows <- as.dist(1 - cor(t(quant_mtx), use = "pairwise.complete.obs"))
  
  pheatmap(
    quant_mtx,
    scale = 'row',
    clustering_method = 'ward.D2',
    annotation_names_col = FALSE,
    color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
    annotation_col = sample_annotations,
    clustering_distance_rows = dist_rows,
    clustering_distance_cols = dist_cols,
    fontsize = 8
  )
}

## ---------------------------------------------------------------------------------------------------------
## Heatmap for the 50 most variable proteins
## ---------------------------------------------------------------------------------------------------------
var <- rowVars(assay(dia_qf[["norm_proteins_replicated"]]), na.rm = TRUE)

most_variable <- sort(var, decreasing = TRUE)[1:50] %>% names()

plot_heatmap_for_uids(most_variable)

## ---------------------------------------------------------------------------------------------------------
## Challenge: Heatmaps of significant proteins
##
## Use plot_heatmap_for_uids (defined above) and the F-test and pairwise
## contrast results to produce heatmaps for:
##
## 1. Proteins with a significant overall difference across groups according
##    to the F-test (adj.P.Val < 0.05)
## 2. Proteins significant (adj.P.Val < 0.05) in at least one pairwise contrast
##
## For each heatmap consider:
## - How do the proteins and samples cluster?
## - Do the sample groups separate cleanly?
## - How does restricting to significant proteins change the heatmap compared
##   to the most variable proteins plotted above?
##
## Hint: Filter limma_results_F and limma_results_all_contrasts to get UniProt
## IDs of significant proteins, then pass these to plot_heatmap_for_uids.
## ---------------------------------------------------------------------------------------------------------
## Format results

limma_results_F$UniprotID

limma_results_all_contrasts_list$MPXV_control$UniprotID

limma_results_all_contrasts_list$Covid19_control$UniprotID

limma_results_all_contrasts_list$MPXV_Covid19$UniprotID

plot_heatmap_for_uids(limma_results_F[limma_results_F$adj.P.Val <0.05,]$UniprotID)

plot_heatmap_for_uids(limma_results_all_contrasts_list$MPXV_control[limma_results_all_contrasts_list$MPXV_control$adj.P.Val <0.05,]$UniprotID)

plot_heatmap_for_uids(limma_results_all_contrasts_list$Covid19_control[limma_results_all_contrasts_list$Covid19_control$adj.P.Val <0.05,]$UniprotID)

plot_heatmap_for_uids(limma_results_all_contrasts_list$MPXV_Covid19[limma_results_all_contrasts_list$MPXV_Covid19$adj.P.Val <0.05,]$UniprotID)

sig_prot_contrasts <- limma_results_all_contrasts %>%
  filter(adj.P.Val < 0.05) %>%
  pull(UniprotID) %>%
  unique()

plot_heatmap_for_uids(sig_prot_contrasts)
## ---------------------------------------------------------------------------------------------------------
## Comparing log fold changes across contrasts
## ---------------------------------------------------------------------------------------------------------

## pivot_wider reshapes results so each protein occupies one row, with separate
## columns for the logFC and adj.P.Val from each contrast.
compare_limma_results <- limma_results_all_contrasts %>%
  select(UniprotID, Genes, adj.P.Val, logFC, contrast) %>%
  pivot_wider(names_from = contrast, values_from = c(adj.P.Val, logFC))


## Initial scatter plot of logFC values across two contrasts
compare_limma_results %>%
  ggplot(aes(x = logFC_MPXV_control, y = logFC_Covid19_control)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Log Fold Change (MPXV vs Control)",
    y = "Log Fold Change (COVID-19 vs Control)"
  )

## Categorise each protein by significance status across the two contrasts
compare_limma_results <- compare_limma_results %>%
  mutate(
    sig = case_when(
      adj.P.Val_MPXV_control < 0.05 & adj.P.Val_Covid19_control < 0.05 ~ "Both",
      adj.P.Val_MPXV_control < 0.05 &
        adj.P.Val_Covid19_control >= 0.05 ~ "MPXV only",
      adj.P.Val_MPXV_control >= 0.05 &
        adj.P.Val_Covid19_control < 0.05 ~ "COVID-19 only",
      TRUE ~ "None"
    )
  )

## Scatter plot coloured by significance category, labelling proteins significant
## in both contrasts
p <- compare_limma_results %>%
  arrange(desc(sig)) %>%
  ggplot(aes(x = logFC_MPXV_control, y = logFC_Covid19_control)) +
  geom_point(aes(colour = sig), pch = 20, size = 3) +
  geom_text_repel(
    data = filter(compare_limma_results, sig == "Both"),
    aes(label = Genes),
    size = 3,
    colour = 'grey20',
    max.overlaps = 100
  ) +
  theme_minimal() +
  labs(
    x = "Log Fold Change (MPXV vs Control)",
    y = "Log Fold Change (COVID-19 vs Control)"
  ) +
  scale_colour_manual(
    values = c(
      "Both" = "purple",
      "MPXV only" = "orangered4",
      "COVID-19 only" = "steelblue",
      "None" = "grey70"
    )
  )

print(p)

## Add transparency and a linear regression line
p2 <- p +
  aes(alpha = sig) +
  scale_alpha_manual(
    values = c(
      "Both" = 0.8,
      "MPXV only" = 0.8,
      "COVID-19 only" = 0.8,
      "None" = 0.25
    )
  ) +
  geom_smooth(
    aes(x = logFC_MPXV_control, y = logFC_Covid19_control),
    inherit.aes = FALSE, # Don't inherit aesthetics to avoid inheriting alpha aesthetics
    method = 'lm',
    se = FALSE,
    linetype = 2,
    colour = 'grey'
  ) +
  theme(legend.title = element_blank())

print(p2)
