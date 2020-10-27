# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(dplyr)
library(ggplot2)

l6.immune <- connect(filename = "/home/nicola/Scrivania/DATASET_MOUSE/l1_olfactory.loom", mode = "r+")
l6.immune
l6.seurat <- as.Seurat(l6.immune)
# VlnPlot(l6.seurat, features = c("Ace2", "Tmprss2", "Snap25", "Ctrb", "Ctsl", "Bsg", "Hspa5", "Dpp4",
#                                 "Furin", "Anpep", "Tmprss11d", "St6gal1", "St3gal4",
#                                 "Ceacam1", "Th", "Mlc1", "Itgam", "Top2a", "Olig2", "Foxc1", "Mog",
#                                 "Kcnj8", "Frzb"), ncol = 2, pt.size = 0.1)

VlnPlot(l6.seurat, features = c("Ace2", "Tmprss2", "Snap25", "Th"), ncol = 2, pt.size = 0.1)

plot1 <- FeatureScatter(object = l6.seurat, feature1 = "Tmprss2", feature2 = "Snap25")
plot1
plot2 <- FeatureScatter(object = l6.seurat, feature1 = "Tmprss2", feature2 = "Ace2")
plot2

# Normalize
mouse <- NormalizeData(object = l6.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
mouse

# Identification of highly variable features 
mouse <- FindVariableFeatures(object = mouse, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = mouse), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = mouse)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# # Scaling Data
all.genes <- rownames(x = mouse)
mouse <- ScaleData(object = mouse)


# Perform linear dimensional reduction
mouse <- RunPCA(object = mouse, features = VariableFeatures(object = mouse))
mouse


mouse <- FindNeighbors(object = mouse, dims = 1:18)
mouse <- FindClusters(object = mouse, resolution = 0.5)

mouse <- RunUMAP(object = mouse, dims = 1:18)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(object = mouse, reduction = "umap", label = TRUE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mouse.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(mouse, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(mouse, features = c("Ace2", "Tmprss2", "Snap25", "Th"))

# you can plot raw counts as well
#VlnPlot(patient, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(mouse, features = c("Ace2", "Tmprss2", "Snap25", "Ctrb", "Ctsl", "Bsg", "Hspa5", "Dpp4", 
                                  "Furin", "Anpep", "Tmprss11d", "St6gal1", "St3gal4",
                                  "Ceacam1", "Th", "Mlc1", "Itgam", "Top2a", "Olig2", "Foxc1", "Mog", 
                                  "Kcnj8", "Frzb"))
                                 


top10 <- mouse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mouse, features = top10$gene) 


new.cluster.ids <- c("VIP", "Calretinin", "Oligodendrocytes", "Mitral/Tufted", "ETCs", "Interneurons", "Neurons", "Granule",
                     "Immature", "Dopamine", "Immune", "Pericytes", "")
names(new.cluster.ids) <- levels(mouse)[1:12]
mouse <- RenameIdents(mouse, new.cluster.ids)
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5) 

saveRDS(mouse, file = "/home/nicola/Scrivania/DATASET_MOUSE/final.rds")













