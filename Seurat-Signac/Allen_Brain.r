library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(readr)



cell_metadata <- read_csv("cell_metadata_with_cluster_annotation.csv")
cell_metadata=cell_metadata[cell_metadata$region_of_interest_acronym=="ACA",]



df1 <- LoadH5Seurat("WMB-10Xv2-Isocortex-1-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df1[["cellnames"]]<-colnames(df1)
cellnames<-colnames(df1)
HPF_metadata1<-cell_metadata[cell_metadata$cell_label %in% cellnames, ]
cellnames1<-HPF_metadata1$cell_label
df1<- subset(df1, subset = cellnames %in% cellnames1 )
df1<-AddMetaData(df1, HPF_metadata1)



df2 <- LoadH5Seurat("WMB-10Xv2-Isocortex-2-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df2[["cellnames"]]<-colnames(df2)
cellnames<-colnames(df2)
HPF_metadata2<-cell_metadata[cell_metadata$cell_label %in% cellnames, ]
cellnames1<-HPF_metadata2$cell_label
df2<- subset(df2, subset = cellnames %in% cellnames1 )
df2<-AddMetaData(df2, HPF_metadata2)



df3 <- LoadH5Seurat("WMB-10Xv2-Isocortex-3-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df3[["cellnames"]]<-colnames(df3)
cellnames<-colnames(df3)
HPF_metadata3<-cell_metadata[cell_metadata$cell_label %in% cellnames, ]
cellnames1<-HPF_metadata3$cell_label
df3<- subset(df3, subset = cellnames %in% cellnames1 )
df3<-AddMetaData(df3, HPF_metadata3)

df4 <- LoadH5Seurat("WMB-10Xv2-Isocortex-4-raw.h5seurat",  meta.data = FALSE, misc = FALSE)
df4[["cellnames"]]<-colnames(df4)
cellnames<-colnames(df4)
HPF_metadata4<-cell_metadata[cell_metadata$cell_label %in% cellnames, ]
cellnames1<-HPF_metadata4$cell_label
df4<- subset(df4, subset = cellnames %in% cellnames1 )
df4<-AddMetaData(df4, HPF_metadata4)

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

WMB_10x_seurat <- NormalizeData(WMB_10xV2_ACA)

WMB_10x_seurat <- FindVariableFeatures(WMB_10x_seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(WMB_10x_seurat)
WMB_10x_seurat <- ScaleData(WMB_10x_seurat, features = all.genes)
WMB_10x_seurat <- RunPCA(WMB_10x_seurat, features = VariableFeatures(object = WMB_10x_seurat))

##
pdf("ElbowPlot.pdf",width=12,height=13)
ElbowPlot(WMB_10x_seurat)
dev.off()

WMB_10x_seurat <- RunUMAP(WMB_10x_seurat, dims = 1:10)

pdf("DimPlot.allen_brain_ACA.class.pdf",width=12,height=13)
DimPlot(WMB_10x_seurat, group.by = "class", label = TRUE, repel = TRUE)+NoLegend()
dev.off()

Idents(WMB_10x_seurat)<-"class"

markers<-FindAllMarkers(WMB_10x_seurat, only.pos = TRUE)


write.csv(markers, "ACA_10Xv2_seurat_class_markers.csv" )



Idents(WMB_10x_seurat)<-"subclass"

pdf("DimPlot.allen_brain_ACA.subclass.pdf",width=12,height=13)
DimPlot(WMB_10x_seurat, group.by = "subclass", label = TRUE, repel=TRUE)+NoLegend()
dev.off()

markers<-FindAllMarkers(WMB_10x_seurat, only.pos = TRUE)

write.csv(markers, "ACA_10Xv2_seurat_subclass_markers.csv" )


saveRDS(WMB_10x_seurat, "ACA_10Xv2_seurat.rds")
