library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/HarderLab/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(x = pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(object = pbmc)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)

pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(x = pbmc[["pca"]])

VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca")

DimPlot(object = pbmc, reduction = "pca")

DimHeatmap(object = pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(object = pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)

JackStrawPlot(object = pbmc, dims = 1:15)

ElbowPlot(object = pbmc)

pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(x = Idents(object = pbmc), 5)
#library(reticulate)

# # create a new environment 
# conda_create("r-reticulate")
# 
# # install SciPy
# conda_install("r-reticulate", "umap-learn")
# 
# # import SciPy (it will be automatically discovered in "r-reticulate")
# umap <- import("umap-learn")
# # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# # 'umap-learn')
# pbmc <- RunUMAP(object = pbmc, dims = 1:10)
# # note that you can set `label = TRUE` or use the LabelClusters function to help label
# # individual clusters

pbmc <- RunTSNE(object=pbmc)
TSNEPlot(object = pbmc)
saveRDS(pbmc, file = "~/HarderLab/pbmc_tutorial.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(x = cluster5.markers, n = 5)

VlnPlot(object = pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(object = pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", 
                                        "PPBP", "CD8A"), cols = c("lightgrey", "black"),  blend.threshold = 0.3)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()



new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Mk")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids)
TSNEPlot(object = pbmc,  label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "~/HarderLab/pbmc_tutorial_final.rds")
