## =============================================================
## Lesson 06: Statistical analysis
## =============================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("limma")
library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("patchwork")


## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------
## Handling missing values prior to testing
## ---------------------------------------------------------------------------------------------------------

## We only test proteins quantified in >= 3 replicates of each condition.
## getWithColData extracts a SummarizedExperiment from the QFeatures object,
## carrying the colData across so we can subset directly by group.

dia_qf[["norm_proteins"]] 
se <- getWithColData(dia_qf, "norm_proteins")
print(se) # same but with colData names

## Subset by condition
control <- se[, se$group == "control"]
mpxv    <- se[, se$group == "MPXV"]
cov     <- se[, se$group == "Covid19"]

## Keep only proteins quantified in >= 3 replicates per condition
control_replicated <- rowSums(!is.na(assay(control))) >= 3
mpxv_replicated    <- rowSums(!is.na(assay(mpxv))) >= 3
cov_replicated     <- rowSums(!is.na(assay(cov))) >= 3

tokeep <- control_replicated & mpxv_replicated & cov_replicated
print(table(tokeep))

## Add new set with only the replicated proteins
dia_qf <- addAssay(dia_qf, y = se[tokeep, ], name = "norm_proteins_replicated")

dia_qf <- addAssayLink(
  dia_qf,
  from = "norm_proteins",
  to = "norm_proteins_replicated"
)

dia_qf
plot(dia_qf)
## ---------------------------------------------------------------------------------------------------------
## Defining the statistical model
## ---------------------------------------------------------------------------------------------------------

## No-intercept design estimates the mean for each group directly.
## This makes it straightforward to specify all pairwise contrasts.
dia_qf$group

group <- factor(dia_qf$group, levels = c("control", "MPXV", "Covid19"))
group

limma_design <- model.matrix(formula(~ 0 + group))
limma_design

# Rename design matrix columns to make them easier to refer to
colnames(limma_design) <- levels(group)

## Verify the design matrix
limma_design

## ---------------------------------------------------------------------------------------------------------
## Specifying contrasts
## ---------------------------------------------------------------------------------------------------------
contrasts_matrix <- makeContrasts(
  MPXV_control = MPXV - control,
  Covid19_control = Covid19 - control,
  MPXV_Covid19 = MPXV - Covid19,
  levels = limma_design
)

## Verify the contrasts matrix
contrasts_matrix


## ---------------------------------------------------------------------------------------------------------
## Running the empirical Bayes-moderated test using limma
## ---------------------------------------------------------------------------------------------------------
data <- assay(dia_qf[['norm_proteins_replicated']])

limma_fit <- lmFit(data, limma_design)

limma_fit_contrasts <- contrasts.fit(
  fit = limma_fit,
  contrasts = contrasts_matrix
)

limma_fit_contrasts <- eBayes(limma_fit_contrasts, trend = TRUE, robust = TRUE)

## ---------------------------------------------------------------------------------------------------------
## Diagnostic plots
## ---------------------------------------------------------------------------------------------------------

## SA plot: residual SD vs log abundance — trend line should be flat or gently
## decreasing; a strong increasing trend suggests trend = TRUE was appropriate.
plotSA(fit = limma_fit_contrasts, cex = 0.5, xlab = "Average log2 abundance")


## ---------------------------------------------------------------------------------------------------------
## F-test for overall significance across groups
## ---------------------------------------------------------------------------------------------------------

## coef = NULL returns the F-statistic for overall significance across all groups.
limma_results_F <- topTable(
  fit = limma_fit_contrasts,
  coef = NULL,
  adjust.method = "BH", # Method for multiple hypothesis testing
  number = Inf
) %>% # Print results for all proteins
  rownames_to_column("UniprotID")

head(limma_results_F)

# Check distribution of p-values
limma_results_F %>%
  ggplot(aes(x = P.Value)) +
  geom_histogram() +
  theme_classic() +
  ggtitle("Distribution of p-values for F-test")
## ---------------------------------------------------------------------------------------------------------
## Pairwise contrasts: MPXV vs Control
## ---------------------------------------------------------------------------------------------------------

## Extract gene name annotations from rowData to merge into results
limma_results_mpxv_control <- topTable(
  limma_fit_contrasts,
  coef = "MPXV_control",
  number = Inf,
  sort.by = "p"
) %>%
  data.frame() %>%
  rownames_to_column('UniprotID')
view(limma_results_mpxv_control)

uniprot_annotations <- rowData(dia_qf[['norm_proteins_replicated']]) %>%
  data.frame() %>%
  select(Genes, Protein.Names)
view(uniprot_annotations)

limma_results_mpxv_control <- limma_results_mpxv_control %>%
  merge(uniprot_annotations, by.x = 'UniprotID', by.y = 'row.names') %>%
  arrange(P.Value)
view(limma_results_mpxv_control)

# plot
limma_results_mpxv_control %>%
  ggplot(aes(P.Value)) +
  geom_histogram() +
  theme_classic() +
  ggtitle("P-value distribution for MPXV vs Control contrast")

## How many significant changes?
sig <- limma_results_mpxv_control$adj.P.Val < 0.05
limma_results_mpxv_control$adj.P.Val[sig]
nrow(limma_results_mpxv_control[sig,])
table(sig, ifelse(limma_results_mpxv_control$logFC > 0, 'Increased', 'Decreased'))


## ---------------------------------------------------------------------------------------------------------
## Volcano plot: MPXV vs Control
## ---------------------------------------------------------------------------------------------------------

# extract data for highlighting
data_for_highlighting <- limma_results_mpxv_control %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(P.Value) %>%
  head(30)

# plot
ggplot(
  limma_results_mpxv_control,
  aes(
    logFC,
    -log10(P.Value),
    colour = adj.P.Val < 0.05,
    alpha = adj.P.Val < 0.05
  )
) +
  geom_point(pch = 20, size = 3) +
  ggrepel::geom_text_repel(
    data = data_for_highlighting,
    aes(label = Genes),
    vjust = 1.5,
    size = 2,
    colour = 'grey20'
  ) +
  theme_classic() +
  xlab("Log2 Fold Change (MPXV vs Control)") +
  ylab("-Log10 P-value") +
  scale_colour_manual(values = c('grey', 'orangered3')) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(legend.position = "none")

## ---------------------------------------------------------------------------------------------------------
## Extracting results for all contrasts
## ---------------------------------------------------------------------------------------------------------
contrasts <- colnames(contrasts_matrix)

limma_results_all_contrasts_list <- vector('list', length(contrasts))
names(limma_results_all_contrasts_list) <- contrasts

for (contrast in contrasts) {
  limma_results_all_contrasts_list[[contrast]] <- topTable(
    limma_fit_contrasts,
    coef = contrast,
    number = Inf,
    sort.by = "p"
  ) %>%
    data.frame() %>%
    rownames_to_column('UniprotID') %>%
    merge(uniprot_annotations, by.x = 'UniprotID', by.y = 'row.names') %>%
    mutate(contrast = contrast) %>%
    arrange(P.Value)
}

## ---------------------------------------------------------------------------------------------------------
## P-value distributions across all contrasts
## ---------------------------------------------------------------------------------------------------------
# bind into one data.frame
limma_results_all_contrasts <- bind_rows(limma_results_all_contrasts_list)

limma_results_all_contrasts %>%
  ggplot(aes(P.Value)) +
  geom_histogram() +
  theme_classic() +
  labs(x = "P-value", y = "Frequency") +
  facet_wrap(~contrast)

## ---------------------------------------------------------------------------------------------------------
## Global multiple testing correction across all contrasts
## ---------------------------------------------------------------------------------------------------------

## Adjusting globally across all proteins and contrasts is more rigorous than
## correcting each contrast separately. BH correction is valid under the positive
## correlation structure typical of multi-contrast proteomics experiments.
limma_results_all_contrasts <- limma_results_all_contrasts %>%
  mutate(adj.P.Val = p.adjust(P.Value, method = 'BH'))


## How many significant changes in each contrast?
table(
  sig = limma_results_all_contrasts$adj.P.Val < 0.05,
  ifelse(limma_results_all_contrasts$logFC > 0, 'Increased', 'Decreased'),
  limma_results_all_contrasts$contrast
)

## ---------------------------------------------------------------------------------------------------------
## Volcano plots for all contrasts
## ---------------------------------------------------------------------------------------------------------

# extract data for highlighting
data_for_highlighting_all_contrasts <- limma_results_all_contrasts %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(P.Value) %>%
  group_by(contrast) %>% # group by contrast to get top 30 for each contrast
  slice_min(P.Value, n = 30)

# plot
ggplot(
  limma_results_all_contrasts,
  aes(
    logFC,
    -log10(P.Value),
    colour = adj.P.Val < 0.05,
    alpha = adj.P.Val < 0.05
  )
) +
  geom_point(pch = 20, size = 3) +
  ggrepel::geom_text_repel(
    data = data_for_highlighting_all_contrasts,
    aes(label = Genes),
    vjust = 1.5,
    size = 2,
    colour = 'grey20'
  ) +
  theme_classic() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  scale_colour_manual(values = c('grey', 'orangered3')) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(legend.position = "none") +
  facet_wrap(~contrast) # facet by contrast to get separate volcano plots for each contrast

## ---------------------------------------------------------------------------------------------------------
## Challenge: Using treat to test against a fold-change threshold
##
## The eBayes function tests whether log fold change equals zero. The treat
## function instead tests whether the absolute log fold change is less than a
## specified threshold, providing a statistically rigorous alternative to
## post-hoc fold-change filtering.
##
## 1. Replace the eBayes call with treat(fc = 1.2, trend = TRUE, robust = TRUE)
##    and extract results for all contrasts using topTreat instead of topTable.
## 2. Apply global Benjamini-Hochberg correction to the p-values as before.
## 3. How many proteins are significant (adj.P.Val < 0.05) in each contrast
##    compared to the eBayes results? Why is the number lower?
## ---------------------------------------------------------------------------------------------------------
limma_fit_treat <- treat(
  limma_fit_contrasts,
  fc = 1.2,
  trend = TRUE,
  robust = TRUE
)

limma_results_treat_list <- lapply(
  colnames(contrasts_matrix),
  function(contrast) {
    topTreat(limma_fit_treat, coef = contrast, n = Inf) %>%
      rownames_to_column("UniprotID") %>%
      mutate(contrast = contrast)
  }
)

limma_results_treat_all <- bind_rows(limma_results_treat_list) %>%
  mutate(adj.P.Val = p.adjust(P.Value, method = 'BH'))

table(
  sig = limma_results_treat_all$adj.P.Val < 0.05,
  limma_results_treat_all$contrast
)
