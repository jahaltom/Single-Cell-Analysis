library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(readr)



cell_metadata <- read_csv("cell_metadata_with_cluster_annotation.csv")
cell_metadata=cell_metadata[cell_metadata$region_of_interest_acronym=="ACA",]


df1 <- LoadH5Seurat("WMB-10Xv2-Isocortex-1-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df1 <- subset(df1, cells = cell_metadata$cell_label)

df2 <- LoadH5Seurat("WMB-10Xv2-Isocortex-2-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df2 <- subset(df2, cells = cell_metadata$cell_label)

df3 <- LoadH5Seurat("WMB-10Xv2-Isocortex-3-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df3 <- subset(df3, cells = cell_metadata$cell_label)

df4 <- LoadH5Seurat("WMB-10Xv2-Isocortex-4-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df4 <- subset(df4, cells = cell_metadata$cell_label)




####Gene names
gene_symbols<-df1@assays[["RNA"]]@meta.features[["gene_symbol"]]

df<-data.frame(gene_symbols, duplicate= duplicated(gene_symbols))

head(df)

df$duplicate<-replace(df$duplicate, df$duplicate =="FALSE", "")

df$concatenate<-paste0(df$gene_symbols, df$duplicate)

df$duplicate_2<- duplicated(df$concatenate)

df$duplicate_2<-replace(df$duplicate_2, df$duplicate_2 =="FALSE", "")

df$gene_symbol<-paste0(df$concatenate, df$duplicate_2)


##Recreate
df1<-LayerData(df1, assay = "RNA", layer = "counts")
df1@Dimnames[[1]]<-df$gene_symbol
df1<-CreateSeuratObject(counts = df1, project = "Allen Brain 10Xv2 2023", meta.data = HPF_metadata1)

df2<-LayerData(df2, assay = "RNA", layer = "counts")
df2@Dimnames[[1]]<-df$gene_symbol
df2<-CreateSeuratObject(counts = df2, project = "Allen Brain 10Xv2 2023", meta.data = HPF_metadata2)

df3<-LayerData(df3, assay = "RNA", layer = "counts")
df3@Dimnames[[1]]<-df$gene_symbol
df3<-CreateSeuratObject(counts = df3, project = "Allen Brain 10Xv2 2023", meta.data = HPF_metadata3)

df4<-LayerData(df4, assay = "RNA", layer = "counts")
df4@Dimnames[[1]]<-df$gene_symbol
df4<-CreateSeuratObject(counts = df4, project = "Allen Brain 10Xv2 2023", meta.data = HPF_metadata4)




######################################################################################################################I used df2 to get the gene symbol, its not in WMB_10xV2_HPF. WMB_10xV2_HPF and df2 ssame dim
WMB_10xV2_ACA  <- merge(df1, y = c(df2, df3, df4), add.cell.ids = c("1", "2","3","4"), project = "ACA")




saveRDS(WMB_10x_seurat, "ACA_10Xv2_seurat.rds")
