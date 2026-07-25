## ================================================
## Lesson 04: Normalisation and aggregation
## ================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("naniar")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")


## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------



## ---------------------------------------------------------------------------------------------------------
## Examine the distribution of raw precursor intensities
## ---------------------------------------------------------------------------------------------------------
dia_qf[["precursors_filtered_missing"]] %>%
  assay() %>%
  longForm() %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  theme_bw() +
  xlab("Abundance (raw)")

## ---------------------------------------------------------------------------------------------------------
## Log2 transformation
## ---------------------------------------------------------------------------------------------------------

## Log transformation is applied before aggregation because robustSummary
## fits a linear model that assumes approximately Gaussian residuals.
dia_qf <- logTransform(
  dia_qf,
  base = 2,
  i = "precursors_filtered_missing",
  name = "precursors_filtered_missing_log"
)

dia_qf[["precursors_filtered_missing_log"]] %>%
  assay() %>%
  longForm() %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  theme_bw() +
  xlab("Abundance (raw)")

## ---------------------------------------------------------------------------------------------------------
## Summarising to protein-level abundance
## ---------------------------------------------------------------------------------------------------------

## robustSummary handles missing values at the precursor level by modelling
## log-transformed precursor quantities as protein abundance + precursor effect.
## A precursor must be quantified in at least 2 samples.
dia_qf <- aggregateFeatures(
  dia_qf,
  i = "precursors_filtered_missing_log",
  fcol = "Protein.Ids",
  name = "proteins",
  fun = MsCoreUtils::robustSummary,
  maxit = 10000
)
head(rowData(dia_qf[["proteins"]])) #rowData is now at Protein level
names(rowData(dia_qf[["proteins"]]))

## ---------------------------------------------------------------------------------------------------------
## Challenge 1: Quantification completeness at the protein level
##
## Although robustSummary handles precursor-level missing values, a protein
## will still be missing in any sample where ALL of its precursors are absent.
##
## Use longForm() to convert the "proteins" set to long format (one row per
## protein-sample combination). The colvars argument carries colData columns
## across. Run:
##
##   longForm(dia_qf[,,'proteins'], colvars = 'group') %>%
##     data.frame() %>%
##     head()
##
## Then create a plot showing how many samples each protein is quantified in,
## broken down by group.
##
## Hint: Count finite values per protein per group using sum(is.finite(value)).
## ---------------------------------------------------------------------------------------------------------
# dia_qf[,,'proteins'] subsets the QFeatures object to the "proteins" set
longForm(dia_qf[,, 'proteins'], colvars = 'group') %>%
  data.frame() %>%
  head()

longForm(dia_qf[,, 'proteins'], colvars = 'group') %>%
  data.frame() %>%
  group_by(rowname, group) %>%
  summarise(n_quant = sum(is.finite(value))) %>%
  ggplot(aes(n_quant)) +
  geom_histogram() +
  facet_wrap(~group) +
  theme_classic() +
  labs(x = 'Quantified samples', y = 'Proteins')

## ---------------------------------------------------------------------------------------------------------
## Challenge 2: Patterns of missingness
##
## Use naniar::gg_miss_upset to visualise patterns of missingness in the
## protein-level data. An upset plot shows which combinations of samples
## tend to have missing values for the same proteins.
##
## - Are missing values random across samples, or do they cluster within
##   particular groups?
## - Is the missingness more consistent with MCAR or MNAR?
##
## Hint: Extract the protein assay with assay(), convert to a data frame,
## and pass to naniar::gg_miss_upset. Use nintersects to limit the number
## of intersections shown.
## ---------------------------------------------------------------------------------------------------------
missing_data <- dia_qf[['proteins']] %>% assay() %>% data.frame()

# The sets argument is required to make the keep_order work as intended,
# otherwise the order of the sets is determined by the total number of missing
# values in each set, which is not what we want here. By providing the sets
# argument with the column names of the missing data, we can ensure that the order
# of the sets in the upset plot matches the order of the columns in our data frame.
naniar::gg_miss_upset(
  missing_data,
  sets = paste0(colnames(missing_data), "_NA"),
  keep.order = TRUE,
  nintersects = 50
)



## ---------------------------------------------------------------------------------------------------------
## Normalisation
## ---------------------------------------------------------------------------------------------------------

## diff.median shifts each sample's intensity distribution so that all sample
## medians match the grand median across all samples.
dia_qf <- normalize(
  dia_qf,
  i = "proteins",
  name = "norm_proteins",
  method = "diff.median"
)

pre_norm <- longForm(dia_qf[,, 'proteins'], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Pre-normalisation')

post_norm <- longForm(dia_qf[,, 'norm_proteins'], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Post-normalisation')

pre_norm + post_norm


## Visualise the effect of normalisation with density plots

longForm(dia_qf[,,'proteins'], colvars = 'group')
rowData(dia_qf[[5]])

