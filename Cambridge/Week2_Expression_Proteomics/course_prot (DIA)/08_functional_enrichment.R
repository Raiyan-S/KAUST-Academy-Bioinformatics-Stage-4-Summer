## =============================================================
## Lesson 08: Functional enrichment analysis
## =============================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("ggplot2")
library("clusterProfiler")
library("tidyr")
library("tibble")
library("org.Hs.eg.db")
library("enrichplot")
library("GOSemSim")
library("pheatmap")
library("dplyr")

limma_results_all_contrasts_list$MPXV_control[limma_results_all_contrasts_list$MPXV_control$adj.P.Val <0.05,]$Genes
## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------
## Subset results to the contrast of interest
## ---------------------------------------------------------------------------------------------------------

## We focus on COVID-19 vs Control as the comparison of primary biological interest.
limma_results_covid_control <- limma_results_all_contrasts %>%
  filter(contrast == 'Covid19_control')

## ---------------------------------------------------------------------------------------------------------
## GO over-representation analysis
## ---------------------------------------------------------------------------------------------------------

## We test proteins with significantly decreased abundance in COVID-19 vs Control.
## The universe is all proteins that were tested (all quantified proteins).
sig_down <- limma_results_covid_control %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  pull(UniprotID)

go_down <- enrichGO(
  gene = sig_down, # list of down proteins
  universe = limma_results_covid_control$UniprotID, # all proteins
  OrgDb = org.Hs.eg.db, # database to query
  keyType = "UNIPROT", # protein ID encoding
  pvalueCutoff = 0.05,
  ont = "ALL", # can be CC, MF, BP, or ALL
  readable = TRUE
)

## ---------------------------------------------------------------------------------------------------------
## Visualising significant GO terms
## ---------------------------------------------------------------------------------------------------------
dotplot(
  go_down,
  x = "Count",
  font.size = 10,
  showCategory = Inf,
  label_format = 100,
  split = "ONTOLOGY",
  color = "p.adjust"
) +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y')


## ---------------------------------------------------------------------------------------------------------
## Reducing GO term redundancy using semantic similarity
## ---------------------------------------------------------------------------------------------------------

## GOSemSim computes pairwise semantic similarity between GO terms.
## simplify removes highly similar terms, keeping the most significant
## representative from each cluster.
gd <- godata('org.Hs.eg.db', ont = "BP")

go_down_similarity <- pairwise_termsim(go_down, method = "Wang", semData = gd)

go_down_simplified <- simplify(
  go_down_similarity,
  cutoff = 0.7, # similarity threshold for clustering
  by = "p.adjust",
  select_fun = min
)

## Re-visualise with reduced redundancy
dotplot(
  go_down_simplified,
  x = "Count",
  font.size = 10,
  showCategory = Inf,
  label_format = 100,
  split = "ONTOLOGY",
  color = "p.adjust"
) +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y')


## ---------------------------------------------------------------------------------------------------------
## Visualising relationships between GO terms with a treeplot
## ---------------------------------------------------------------------------------------------------------
treeplot(
  go_down_similarity,
  showCategory = Inf,
  nCluster = 3,
  fontsize_tiplab = 2,
  fontsize_cladelab = 3,
  cladelab_offset = 6,
  tiplab_offset = 0.5,
  hexpand = 0.2
)


## ---------------------------------------------------------------------------------------------------------
## Linking GO results back to protein abundances
## ---------------------------------------------------------------------------------------------------------

## Extract gene symbols for a GO term of interest: "plasma lipoprotein particle clearance"
## Because we set readable = TRUE in enrichGO, geneID contains gene symbols.
pm_lipo_part_clear <- go_down@result %>%
  filter(Description == "plasma lipoprotein particle clearance") %>%
  pull(geneID) %>%
  strsplit("/") %>%
  unlist()

go_down@result
## Map gene symbols to UniProt accessions
# Gene symbol -> UniProt map
g2p <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = pm_lipo_part_clear,
  columns = "UNIPROT",
  keytype = "SYMBOL"
) %>%
  filter(UNIPROT %in% limma_results_covid_control$UniprotID)

print(g2p)
## ---------------------------------------------------------------------------------------------------------
## Challenge: Visualising abundances of plasma lipoprotein particle clearance proteins
##
## Use the UniProt IDs retrieved above to assess whether the protein abundances
## for the 'plasma lipoprotein particle clearance' proteins are reliable, then
## examine their abundance patterns across samples.
##
## 1. Use subsetByFeature to subset dia_qf to these proteins. Then plot the
##    precursor-level and protein-level abundances side by side across samples,
##    adapting the approach from lesson 05.
##
## 2. Plot a heatmap of the protein-level abundances for these proteins, with
##    samples annotated by group. Do the samples separate as you would expect
##    from the enrichment result?
##
## Hint: Pass rowvars = 'Genes' to longForm to carry gene names through to
## the plot, then facet by gene name and assay.
## ---------------------------------------------------------------------------------------------------------
# extract UniProt IDs for the proteins driving the enrichment of the 'plasma lipoprotein particle clearance' term
uids <- g2p$UNIPROT


# 1. Subset to the proteins of interest and
# plot precursor and protein level abundances across samples
dia_qf_lipids <- subsetByFeature(dia_qf, uids)

dia_qf_lipids[,, c("precursors_filtered_missing", "proteins")] %>%
  longForm(rowvars = c('Genes')) %>%
  as_tibble() %>%
  mutate(
    assay_order = factor(
      assay,
      levels = c("precursors_filtered_missing", "proteins")
    )
  ) %>%
  mutate(
    value = ifelse(assay == 'precursors_filtered_missing', log2(value), value)
  ) %>%
  ggplot(aes(x = colname, y = value, colour = assay)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = rowname)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7)
  ) +
  facet_grid(Genes ~ assay_order, scales = 'free') +
  labs(x = 'Sample', y = 'Abundance (log2)') +
  theme(aspect.ratio = 1 / 2)


# 2. Protein-level heatmap for plasma lipoprotein particle clearance proteins
sample_annotations <- colData(dia_qf) %>%
  data.frame() %>%
  select(group)

quant_mtx <- assay(dia_qf[["norm_proteins_replicated"]][uids, ])
colnames(quant_mtx) <- dia_qf$runCol

rownames(quant_mtx) <- rowData(dia_qf[["norm_proteins_replicated"]][
  uids,
])$Genes

dist_cols <- as.dist(1 - cor(quant_mtx, use = "pairwise.complete.obs"))
dist_rows <- as.dist(1 - cor(t(quant_mtx), use = "pairwise.complete.obs"))

pheatmap(
  quant_mtx,
  scale = "row",
  annotation_col = sample_annotations,
  clustering_distance_rows = dist_rows,
  clustering_distance_cols = dist_cols,
  clustering_method = 'ward.D2',
  fontsize = 8
)
