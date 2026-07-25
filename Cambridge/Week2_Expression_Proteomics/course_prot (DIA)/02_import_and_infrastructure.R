## ================================================
## Lesson 02: Import and infrastructure
## ================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------
library("arrow")
library("readxl")
library("QFeatures")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")


## ---------------------------------------------------------------------------------------------------------
## Read the DIA-NN parquet report
## ---------------------------------------------------------------------------------------------------------
diann_df <- read_parquet("C:/Users/Raiyan Subedar/Downloads/Kaust BioSummer/Cambridge/New/Week2_Expression_Proteomics/dia/data/monkeypox_plasma_proteomes.parquet")

## Inspect the dimensions and column names
dim(diann_df)
colnames(diann_df)
head(diann_df)

## ---------------------------------------------------------------------------------------------------------
## Import sample metadata
## ---------------------------------------------------------------------------------------------------------
sample_metadata <- read_excel("C:/Users/Raiyan Subedar/Downloads/Kaust BioSummer/Cambridge/New/Week2_Expression_Proteomics/dia/data/monkeypox_metadata.xlsx")
print(sample_metadata)

## ---------------------------------------------------------------------------------------------------------
## Clean Run names
## ---------------------------------------------------------------------------------------------------------
# Clean Run names
diann_df[, "Run"] %>%
  unique() %>%
  head(2)


## DIA-NN pastes the .wiff2 folder and .wiff.scan filenames together in the Run column.
## We strip everything from the first "." onwards to match the metadata sample_id field.
diann_df <- diann_df %>%
  mutate(Run = sub('\\..*', '', Run))

print(head(unique(diann_df$Run), 2))

## Merge the runCol label from the metadata into the DIA-NN data frame
diann_df <- diann_df %>%
  merge(
    sample_metadata[, c("sample_id", "runCol")],
    by.x = "Run",
    by.y = "sample_id"
  )


dim(diann_df)
names(diann_df)
head(diann_df$Precursor.Id)


## ---------------------------------------------------------------------------------------------------------
## Read into a QFeatures object
## ---------------------------------------------------------------------------------------------------------
dia_qf <-
  readQFeaturesFromDIANN(
    diann_df,
    quantCols = "Precursor.Quantity",
    fnames = "Precursor.Id",
    runCol = "runCol",
    colData = sample_metadata
  )


## ---------------------------------------------------------------------------------------------------------
## Exploring the QFeatures object
## ---------------------------------------------------------------------------------------------------------
experiments(dia_qf)
dia_qf[[1]]
head(rowData(dia_qf[[1]]))
head(assay(dia_qf[[1]]))
colData(dia_qf)

## ---------------------------------------------------------------------------------------------------------
## Challenge 1: Precursors per sample
##
## 1. Find a function to determine how many precursors are quantified in each
##    sample. Use it to plot a histogram of precursor counts across samples.
##
## 2. Create a bar plot showing the number of precursors identified in each
##    sample, coloured by group. Which group tends to have the most precursors
##    identified? Is there any sample that looks like an outlier?
##
## Hint: Browse the QFeatures documentation with ?QFeatures to find a suitable
## function for obtaining per-sample feature counts.
## ---------------------------------------------------------------------------------------------------------
?QFeatures
nrows(dia_qf)
hist(nrows(dia_qf), xlab = "Precursors", main = "Precursor counts per sample")

# Extract the numeric vector of quantities for Sample 1
# sample1_quantities <- assay(dia_qf[[1]])[, 1]
# Plot the histogram of precursor intensities
# hist(
#   log2(sample1_quantities),
#   xlab = "log2(Precursor Quantity)",
#   main = "Log2 Precursor Intensity Distribution (Sample 1)",
#   col = "steelblue",
#   breaks = 50
# )
colData(dia_qf)

data.frame(
  n_precursors = nrows(dia_qf),
  sample = names(nrows(dia_qf))
) %>%
  merge(data.frame(colData(dia_qf)), by.x = "sample", by.y = "runCol") %>%
  ggplot(aes(x = sample, y = n_precursors, fill = group)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Number of precursors") +
  coord_flip() +
  theme_classic()

## ---------------------------------------------------------------------------------------------------------
## Exploring rowData quality metrics
## ---------------------------------------------------------------------------------------------------------

## rbindRowData collects rowData from all per-sample sets into a single data frame.
## Merge with colData to get group information for plotting.
rdata <- rbindRowData(dia_qf, i = names(dia_qf))
head(rdata)
cdata <- data.frame(colData(dia_qf))
head(cdata)
cdata$group

adata <- rdata %>%
  data.frame() %>%
  merge(cdata, by.x = "assay", by.y = "runCol")


## Plot FWHM distribution per sample, coloured by group

adata %>%
  ggplot(aes(FWHM, group = Run, colour = group)) +
  geom_density() +
  theme_classic()


## ---------------------------------------------------------------------------------------------------------
## Challenge 2: Exploring rowData quality metrics
##
## The rdata object contains many columns from the DIA-NN rowData beyond FWHM.
##
## 1. Use str(rdata) to identify other numerical quality metric columns.
## 2. Choose one metric and create a plot to visualise its distribution
##    across samples, coloured by group.
## 3. Do any metrics appear to differ systematically between groups, or
##    between individual samples?
## ---------------------------------------------------------------------------------------------------------
str(adata)
adata |> 
  ggplot(aes(log2(Q.Value), group = Run, colour = group)) +
  geom_density() +
  theme_classic()

