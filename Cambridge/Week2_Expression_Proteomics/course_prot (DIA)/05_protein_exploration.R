## =============================================================
## Lesson 05: Exploration and visualisation of protein data
## =============================================================


## ---------------------------------------------------------------------------------------------------------
## Load R/Bioconductor libraries
## ---------------------------------------------------------------------------------------------------------

library("QFeatures")
library("factoextra")
library("tidyr")
library("dplyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("ggbeeswarm")
library("patchwork")


## ---------------------------------------------------------------------------------------------------------
## Load data from previous lesson
## ---------------------------------------------------------------------------------------------------------

load("~/Course_Materials/Week2_Expression_Proteomics/course_proteomics/dia/preprocessed/lesson04.rda")
dia_qf
dia_qf[,,"norm_proteins"]
head(rowData(dia_qf[["norm_proteins"]]))
sampleMap(dia_qf[,,"norm_proteins"])
colData(dia_qf)

head(longForm(dia_qf[,,"norm_proteins"]))
## ---------------------------------------------------------------------------------------------------------
## Plotting quantification values for proteins of interest
## ---------------------------------------------------------------------------------------------------------

## Convert normalised protein-level data to long format, carrying across
## both sample-level and feature-level metadata for plotting.
long_form_protein <- longForm(
  dia_qf[,, "norm_proteins"],
  colvars = c('runCol', 'group'),
  rowvars = c('Genes', 'Protein.Names')
) %>%
  data.frame() %>%
  dplyr::rename(UniprotID = rowname)


## Plot abundance of three proteins of interest across groups
proteins_of_interest <- c('P02766', # TTR (Transthyretin)
                          'P18428',  # LBP (Lipopolysaccharide binding protein)
                          'P13796')  # LCP1 (Lymphocyte cytosolic protein 1)

# plot
long_form_protein %>%
  filter(UniprotID %in% proteins_of_interest) %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  geom_quasirandom(pch=21)+
  #geom_point(pch=21) +
  facet_wrap(~Genes, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))


## ---------------------------------------------------------------------------------------------------------
## Principal Component Analysis (PCA)
## ---------------------------------------------------------------------------------------------------------

## PCA requires a complete data matrix, so we first remove proteins with
## any missing values using filterNA before calling prcomp.
head(assay(dia_qf[["norm_proteins"]])) #Protein rows, sample columns
head(t(assay(dia_qf[["norm_proteins"]]))) #Sample rows, protein columns

protein_pca <- dia_qf[["norm_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp(scale = TRUE, center = TRUE)

summary(protein_pca)

## Scree plot
fviz_screeplot(protein_pca)

## PCA plot coloured by group and shaped by age group
str(protein_pca)
colData(dia_qf)

protein_pca$x %>%
  merge(colData(dia_qf), by = 'row.names') %>%
  data.frame() %>%
  head()

p <- protein_pca$x %>%
  merge(colData(dia_qf), by = 'row.names') %>%
  data.frame() %>%
  ggplot(aes(x = PC1, y = PC2, colour = group, shape = age.group)) +
  geom_point(size = 3) +
  theme_bw()

print(p)
print(p + aes(x = PC3, y = PC4))

## ---------------------------------------------------------------------------------------------------------
## The subsetByFeature function
## ---------------------------------------------------------------------------------------------------------

## subsetByFeature returns a new QFeatures object containing only the data
## for a specified feature across all assay levels.
TTR <- subsetByFeature(dia_qf, "P02766")

experiments(TTR)


## Visualise how precursor data was aggregated to protein level for TTR
longForm(TTR[,, c("precursors_filtered_missing_log", "proteins")]) %>%
  data.frame() %>%
  head()

TTR[,, c("precursors_filtered_missing_log", "proteins")] %>%
  longForm() %>%
  as_tibble() %>%
  ggplot(aes(x = colname, y = value, colour = assay)) +
  geom_point() +
  geom_line(aes(group = rowname)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)
  ) +
  facet_wrap(~assay) +
  labs(x = 'Sample', y = 'Intensity (log2)')


## Centre the values per precursor to make profile similarity clearer
TTR[,, c("precursors_filtered_missing_log", "proteins")] %>%
  longForm() %>%
  as_tibble() %>%
  group_by(rowname) %>%
  mutate(value = value - mean(value, na.rm = TRUE)) %>%
  ggplot(aes(x = colname, y = value, colour = assay)) +
  geom_point() +
  geom_line(aes(group = rowname)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)
  ) +
  facet_wrap(~assay) +
  labs(x = 'Sample', y = 'Intensity (log2)')
