# Load the libraries we will need for this practical
library(Seurat)
library(tidyverse)

# set default ggplot theme
theme_set(theme_classic())

# load the preprocessed Seurat object with 500 cells per sample
seurat_object <- readRDS("RObjects/SCT.500.rds")
seurat_object

#### PCA ####

# Run PCA
seurat_object <- RunPCA(seurat_object, 
                        features = VariableFeatures(seurat_object))
seurat_object@reductions

# check the new reductions slot
Reductions(seurat_object)

# check the variance explained by each PC
Stdev(seurat_object, reduction = "pca")

# plot the variance explained by each PC
ElbowPlot(seurat_object, ndims = 50)

unique(seurat_object$SampleName)
# plot the PCA
DimPlot(seurat_object,
        reduction = "pca")

# plot PC2 and PC3
DimPlot(seurat_object,
        reduction = "pca",
        dims = c(2, 3))

# plot PC1 and PC2 split by sample group
DimPlot(seurat_object,
        reduction = "pca",
        split.by = "SampleGroup")


#### t-SNE ####

# run t-SNE using default options
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca")

# confirm a new reduction was added to the object
Reductions(seurat_object)
seurat_object

# plot the t-SNE
DimPlot(seurat_object,
        reduction = "tsne")


#### Exercise: t-SNE ####

# the following code:
# adds the t-SNE result to the Seurat object
# names this reduction "TSNE_perplex30"
# sets perplexity = 30 (which is the default if we don't specify it)
# uses the first 10 principal components
# sets a seed for reproducibility

seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 30,
                         seed = 123,
                         reduction.name = "TSNE_perplex30",
                         reduction.key = "TSNE30_")
DimPlot(seurat_object,
        reduction = "TSNE_perplex30")

# Part A
# Re-run the algorithm but change the random seed number.
# Do the results change dramatically between runs?
# YOUR CODE HERE
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 30,
                         seed = 42,
                         reduction.name = "TSNE_perplex42",
                         reduction.key = "TSNE42_")
DimPlot(seurat_object,
        reduction = "TSNE_perplex42")


# Part B
# Facet these plots by SampleName to better understand where each marker is mostly expressed
# YOUR CODE HERE
DimPlot(object = seurat_object, reduction = "TSNE_perplex42", split.by = "SampleName", ncol = 4)

# Explore different perplexity values (for example 5 and 500)
# Do you get tighter or looser clusters?
# YOUR CODE HERE
seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 5,
                         seed = 42,
                         reduction.name = "TSNE_perplex5",
                         reduction.key = "TSNE5_")
DimPlot(seurat_object,
        reduction = "TSNE_perplex5")

seurat_object <- RunTSNE(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         perplexity = 500,
                         seed = 42,
                         reduction.name = "TSNE_perplex500",
                         reduction.key = "TSNE500_")
DimPlot(seurat_object,
        reduction = "TSNE_perplex500")


# Instead of colouring by SampleName, colour by expression of known cell markers
# CD79A (B cells)
# CST3 (monocytes)
# CD3D (T cells)
# HBA1 (erythrocytes)
# YOUR CODE HERE
FeaturePlot(
  seurat_object,
  c("CD79A","CST3","CD3D","HBA1"),
  reduction = "TSNE_perplex30"
)

FeaturePlot(
  seurat_object,
  features = c("PC_1","PC_2","PC_3","PC_40"),
  reduction = "TSNE_perplex30"
)


#### UMAP ####

# run UMAP using the first 10 PCs
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10)
seurat_object

# confirm a new reduction was added to the object
Reductions(seurat_object)

# visualise the UMAP
DimPlot(seurat_object,
        reduction = "umap")


#### Exercise: UMAP ####

# Part A
# run the UMAP with 50 neighbours
# YOUR CODE HERE
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         n.neighbors = 50,
                         seed = 42,
                         reduction.name = "UMAP_n50",
                         reduction.key = "UMAPn50_")
seurat_object

# Part B
# visualise the resulting UMAP projection
# YOUR CODE HERE
DimPlot(seurat_object,
        reduction = "UMAP_n50")


# Part C
# run the UMAP with 5 and 500 neighbours and compare the results
# YOUR CODE HERE
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         n.neighbors = 5,
                         seed = 42,
                         reduction.name = "UMAP_n5",
                         reduction.key = "UMAPn5_")

seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         n.neighbors = 500,
                         seed = 42,
                         reduction.name = "UMAP_n500",
                         reduction.key = "UMAPn500_")
seurat_object

DimPlot(seurat_object,
        reduction = "UMAP_n5")
DimPlot(seurat_object,
        reduction = "UMAP_n500")

# Part D
# compare the UMAP projection with the t-SNE projections
# would you prefer one over the other?
# YOUR CODE HERE
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:10,
                         n.neighbors = 30,
                         seed = 123,
                         reduction.name = "UMAP_n30",
                         reduction.key = "UMAP_n30")

DimPlot(seurat_object,
        reduction = "UMAP_n30")

DimPlot(seurat_object,
        reduction = "TSNE_perplex30")
