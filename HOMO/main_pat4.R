library(dplyr)
library(Seurat)
library(ggplot2)

# Load dataset
patient.data <- Read10X(data.dir = "/home/nicola/Scrivania/DATASET_HOMO/PATIENT_4")

# Initialize the Seurat object with the raw (non-normalized data).
patient <- CreateSeuratObject(counts = patient.data, project = "patientdata", min.cells = 3, min.features = 200)
patient

# Lets examine a few genes in the first thirty cells, 
#"ACE2", "TMPRSS2", "FOXJ1", "MUC5AC" these are the genes which are considered for the patient 1 and 4 
#pat 2 and 3 -> "ACE2", "TMPRSS2", "GNG8", "GNG13", "KRT5", "CYP2A13"
patient.data[c("ACE2", "TMPRSS2", "FOXJ1", "MUC5AC"), 1:30]

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
patient[["percent.mt"]] <- PercentageFeatureSet(object = patient, pattern = "^MT-")
patient[["percent.mt"]]

# Show QC metrics for the first 5 cells
head(x = patient@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(object = patient, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by
# the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(object = patient, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(object = patient, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
CombinePlots(plots = list(plot1, plot2))

patient <- subset(x = patient, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
patient

# Normalize
patient <- NormalizeData(object = patient, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features 
patient <- FindVariableFeatures(object = patient, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = patient), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = patient)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
#CombinePlots(plots = list(plot1, plot2))

# Scaling Data
all.genes <- rownames(x = patient)
patient <- ScaleData(object = patient, features = all.genes)

# Perform linear dimensional reduction
patient <- RunPCA(object = patient, features = VariableFeatures(object = patient))
patient

# Examine and visualize PCA results a few different ways
print(x = patient[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = patient, dims = 1:2, reduction = "pca")
DimPlot(object = patient, reduction = "pca")
DimHeatmap(object = patient, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = patient, dims = 1:12, cells = 500, balanced = TRUE)

# Each square shows the correlation between the variables on each axis. 
# Correlation ranges from -1 to +1. Values closer to zero means there is no 
# linear trend between the two variables. The close to 1 the correlation 
# is the more positively correlated they are; that is as one increases so 
# does the other and the closer to 1 the stronger this relationship is. 
# A correlation closer to -1 is similar, but instead of both increasing one 
# variable will decrease as the other increases. The diagonals are all 
# 1/dark green because those squares are correlating each variable to itself 
# (so it's a perfect correlation). For the rest the larger the number 
# and darker the color the higher the correlation between the two variables. 
# The plot is also symmetrical about the diagonal since the same two variables
# are being paired together in those squares.

# NOTE: This process can take a long time for big datasets, comment out for expediency. More approximate techniques such
# as those implemented in ElbowPlot() can be used to reduce computation time
patient <- JackStraw(object = patient, num.replicate = 100)
patient <- ScoreJackStraw(object = patient, dims = 1:20)

JackStrawPlot(object = patient, dims = 1:15)


patient <- FindNeighbors(object = patient, dims = 1:10)
patient <- FindClusters(object = patient, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(x = Idents(object = patient), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
patient <- RunUMAP(object = patient, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(object = patient, reduction = "umap", label = TRUE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
patient.markers <- FindAllMarkers(patient, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
patient.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(patient, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(patient, features = c("ACE2", "TMPRSS2", "FOXJ1", "MUC5AC"))

# you can plot raw counts as well
VlnPlot(patient, features = c("ACE2", "TMPRSS2", "FOXJ1", "MUC5AC"), slot = "counts", log = TRUE)

#important UMAP!!!!!!!!!
FeaturePlot(patient, features = c("ACE2", "TMPRSS2", "FOXJ1", "MUC5AC"))


top10 <- patient.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(patient, features = top10$gene) 


new.cluster.ids <- c("Resp. ciliated","Resp. secretory", "Resp. HBC")
names(new.cluster.ids) <- levels(patient)[1:3]
patient <- RenameIdents(patient, new.cluster.ids)
DimPlot(patient, reduction = "umap", label = TRUE, pt.size = 0.5) 

saveRDS(patient, file = "/home/nicola/Scrivania/DATASET_HOMO/PATIENT_4/final.rds")
