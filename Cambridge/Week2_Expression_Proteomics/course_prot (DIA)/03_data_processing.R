## ================================================
## Lesson 03: Data cleaning
## ================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")


## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------
## Filtering by Q-value
## ---------------------------------------------------------------------------------------------------------

## DIA-NN Q-value columns represent the local FDR for each identification.
## We apply a 1% FDR threshold (Q-value <= 0.01) to all four columns.
## Omitting the i argument applies the filter to all per-sample sets simultaneously.

dia_qf

dia_qf <- dia_qf %>%
  filterFeatures(~ Q.Value <= 0.01) %>%       # Run-level precursor Q-value
  filterFeatures(~ PG.Q.Value <= 0.01) %>%    # Run-level protein group Q-value
  filterFeatures(~ Lib.Q.Value <= 0.01) %>%   # Library-level precursor Q-value
  filterFeatures(~ Lib.PG.Q.Value <= 0.01)    # Library-level protein group Q-value

dia_qf # filtered

hist(nrows(dia_qf),
     xlab = "Precursors",
     main = "Precursor counts per sample\nafter Q-value filtering")


## ---------------------------------------------------------------------------------------------------------
## Joining per-sample sets into a single precursor set
## ---------------------------------------------------------------------------------------------------------

## joinAssays merges all per-sample sets into one combined set, aligning by Precursor.Id.
## Where a precursor was not detected in a sample, the entry will be NA.
max(nrows(dia_qf))
dia_qf <- joinAssays(
  x = dia_qf,
  i = names(dia_qf),
  name = "precursors",
  fcol = "Precursor.Id"
)
head(rowData(dia_qf[[32]]))

## Compare rowData columns before and after joining

sample_level_rdata_names <- rowDataNames(dia_qf)[[1]]
sample_level_rdata_names

joined_rdata_names <- rowDataNames(dia_qf)[["precursors"]]
joined_rdata_names

setdiff(sample_level_rdata_names, joined_rdata_names) #metadata rowData that are getting removed

sampleMap(dia_qf)

## ---------------------------------------------------------------------------------------------------------
## Removing individual sample-level assays
## ---------------------------------------------------------------------------------------------------------
sets_to_rm <- which(names(dia_qf) != "precursors")
dia_qf <- removeAssay(dia_qf, i = sets_to_rm)
dia_qf

sampleMap(dia_qf)
## ---------------------------------------------------------------------------------------------------------
## Identifying contaminant proteins
## ---------------------------------------------------------------------------------------------------------

## Contaminant entries from the database search are named with a "Cont_" prefix
## in the Protein.Ids column. We check how many precursors are flagged.
rd <- rowData(dia_qf[['precursors']]) %>%
  data.frame()

head(rd$Protein.Ids)
table(grepl('Cont_', rd$Protein.Ids))
## ---------------------------------------------------------------------------------------------------------
## Filtering contaminants
## ---------------------------------------------------------------------------------------------------------

## We create a copy of the "precursors" set, add an assay link, then filter.
## The original "precursors" set is retained for reference.
dia_qf <- addAssay(dia_qf, dia_qf[["precursors"]], name = "precursors_no_cont")
plot(dia_qf)

dia_qf <- addAssayLink(dia_qf, from = "precursors", to = "precursors_no_cont")
plot(dia_qf)

dia_qf <- filterFeatures(
  dia_qf,
  ~ !grepl("Cont_", Protein.Ids),
  i = "precursors_no_cont"
)
dia_qf
## ---------------------------------------------------------------------------------------------------------
## Exploring missing values
## ---------------------------------------------------------------------------------------------------------
precursors_missing <- nNA(dia_qf, i = "precursors_no_cont")

print(precursors_missing)

## ---------------------------------------------------------------------------------------------------------
## Challenge: Visualising missing values
##
## Using the precursors_missing object, create plots to explore the
## distribution of missing values in the data.
##
## 1. Plot a histogram for the proportion of missing values across precursors.
## 2. Plot the proportion of missing values per sample, coloured by group.
##
## Consider:
## - Do most precursors have missing values?
## - Are there any precursors or samples with an unusually high proportion?
## - Does the proportion of missing values differ between groups?
##
## Hint: precursors_missing contains two data frames: nNArows (one row per
## precursor) and nNAcols (one row per sample), each with a pNA column.
## ---------------------------------------------------------------------------------------------------------
hist(
  precursors_missing$nNArows$pNA,
  xlab = "Proportion missing values",
  main = "Missing values per precursor"
)
precursors_missing$nNAcols %>%
  as_tibble() %>%
  merge(data.frame(colData(dia_qf)), by.x = "name", by.y = "runCol") %>%
  ggplot(aes(x = name, y = pNA, group = group, fill = group)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion missing values") +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic()


## ---------------------------------------------------------------------------------------------------------
## Filtering sparse precursors
## ---------------------------------------------------------------------------------------------------------

## pNA = 0.75 means a precursor must be observed in at least 25% of samples
## (approximately 8 of 31 samples).
dia_qf <- addAssay(
  dia_qf,
  dia_qf[["precursors_no_cont"]],
  name = "precursors_filtered_missing"
)

dia_qf <- addAssayLink(
  dia_qf,
  from = "precursors_no_cont",
  to = "precursors_filtered_missing"
)
plot(dia_qf)

dia_qf <- dia_qf %>%
  filterNA(pNA = 0.75, i = "precursors_filtered_missing")
dia_qf
