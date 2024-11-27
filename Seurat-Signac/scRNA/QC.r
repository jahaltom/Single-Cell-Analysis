library(Signac)
library(Seurat)
library(dplyr)


juvenile=readRDS("06122024_combined_juvenile_integrated.allen_brain_projection.filtered.rds")


DefaultAssay(juvenile) <- 'RNA'

#Low-quality / dying cells often exhibit extensive mitochondrial contamination
juvenile[["percent.mt"]] <- PercentageFeatureSet(juvenile, pattern = "^MT-")

Idents(juvenile)="seuratProject"
png("VlnPlot_QC.png",width=15,height=15,units="in",res=300)
VlnPlot(juvenile, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


juvenile <- subset(juvenile, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Normalizing the data
juvenile <- NormalizeData(juvenile)
# Identification of highly variable features (feature selection)

juvenile <- FindVariableFeatures(juvenile, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(juvenile), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(juvenile)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


png("VarFeat.png",width=15,height=15,units="in",res=300)
plot1 + plot2
dev.off()

#Scaling the data
all.genes <- rownames(juvenile)
juvenile <- ScaleData(juvenile, features = all.genes)



# Perform linear dimensional reduction

juvenile <- RunPCA(juvenile)

# Examine and visualize PCA results a few different ways
print(juvenile[["pca"]], dims = 1:5, nfeatures = 5)


png("Loadings.png",width=15,height=15,units="in",res=300)
VizDimLoadings(juvenile, dims = 1:2, reduction = "pca")
dev.off()


png("Dimplot.png",width=15,height=15,units="in",res=300)
DimPlot(juvenile, reduction = "pca") + NoLegend()
dev.off()


png("Dimheatmap.png",width=15,height=15,units="in",res=300)
DimHeatmap(juvenile, dims = 1, cells = 500, balanced = TRUE)
dev.off()

png("Dimheatmaps.png",width=15,height=15,units="in",res=300)
DimHeatmap(juvenile, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine the ‘dimensionality’ of the dataset 
png("Elbow.png",width=15,height=15,units="in",res=300)
ElbowPlot(juvenile,ndims=50)
dev.off()

# Cluster the cells
juvenile <- FindNeighbors(juvenile, dims = 1:10)
juvenile <- FindClusters(juvenile, resolution = 0.5)


# Run non-linear dimensional reduction (UMAP/tSNE)

juvenile <- RunUMAP(juvenile, dims = 1:10)

png("UMAP.png",width=15,height=15,units="in",res=300)
DimPlot(juvenile, reduction = "umap")
dev.off()


saveRDS(juvenile, file = "06122024_combined_juvenile_integrated.allen_brain_projection.filtered.rds")


