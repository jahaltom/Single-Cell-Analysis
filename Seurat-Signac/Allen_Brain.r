library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(readr)



cell_metadata <- read_csv("cell_metadata_with_cluster_annotation.csv")
cell_metadata=cell_metadata[cell_metadata$region_of_interest_acronym=="ACA",]


df1 <- LoadH5Seurat("WMB-10Xv2-Isocortex-1-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df1[["cell_label"]]=rownames(df1@meta.data)
df1 <- subset(df1, cells = cell_metadata$cell_label)

df2 <- LoadH5Seurat("WMB-10Xv2-Isocortex-2-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df2[["cell_label"]]=rownames(df2@meta.data)
df2 <- subset(df2, cells = cell_metadata$cell_label)

df3 <- LoadH5Seurat("WMB-10Xv2-Isocortex-3-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df3[["cell_label"]]=rownames(df3@meta.data)
df3 <- subset(df3, cells = cell_metadata$cell_label)

df4 <- LoadH5Seurat("WMB-10Xv2-Isocortex-4-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df4[["cell_label"]]=rownames(df4@meta.data)
df4 <- subset(df4, cells = cell_metadata$cell_label)




WMB_10xV2_ACA  <- merge(df1, y = c(df2, df3, df4), project = "ACA")

#Assign gene symbol to metadata
geneSymbol=df1@assays[["RNA"]]@meta.features[["gene_symbol"]]
WMB_10xV2_ACA@assays[["RNA"]]@meta.features[["gene_symbol"]]=geneSymbol

#Add metadata
rownames(cell_metadata)=cell_metadata$cell_label
cell_metadata <- cell_metadata[colnames(WMB_10xV2_ACA), ]
WMB_10xV2_ACA <- AddMetaData(WMB_10xV2_ACA, metadata = cell_metadata)




saveRDS(WMB_10xV2_ACA, "ACA_10Xv2_AllenBrain.rds")
