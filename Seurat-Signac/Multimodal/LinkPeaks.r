library(Signac)
library(ggplot2)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(readr)
#library(DESeq2)
library(dplyr)
library(tibble)
library(EnsDb.Hsapiens.v86)
library("biovizBase")
set.seed(1234)



srtRNA=readRDS("RNA_AllAssays_Final_Subset.rds")
srtATAC <-readRDS("ATAC_AllAssays_Final.rds")


# Merge the RNA and ATAC Seurat objects
srt <- merge(
    x = srtRNA,
    y = srtATAC,
    add.cell.ids = c("RNA", "ATAC"),  # Prefix to distinguish cells from each dataset
    project = "MultiModal"
)




DefaultAssay(srt) <- "RNA"
srt <- NormalizeData(srt)
srt <- FindVariableFeatures(srt)
srt <- ScaleData(srt)

DefaultAssay(srt)<-"ATAC"
#Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
srt <- RunTFIDF(srt, assay = "ATAC")
#Run partial singular value decomposition using irlba. Dimensionality reduction
srt <- RunSVD(srt, assay = "ATAC")



# Add gene annotations to the ATAC assay
Annotation(srt) <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(Annotation(srt)) <- "UCSC"

#Create metadata that has cell type and condition
srt$celltype.stim <- paste(srt$MarkerAnnotations, srt$Status, sep = "_")
#Lable Idents with this metadata
Idents(srt) <- "celltype.stim"





# ##Motif Analysis
# #to avoid an error: https://github.com/timoast/signac/issues/486#issuecomment-788964712
# main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
# keep.peaks <- as.logical(seqnames(granges(srt)) %in% main.chroms)############################################
# srt <- srt[keep.peaks, ]




rna_cells <- colnames(srt[["RNA"]])
atac_cells <- colnames(srt[["ATAC"]])
# Find common cells
common_cells <- intersect(rna_cells, atac_cells)



# Subset the Seurat object to include only the common cells
srt <- subset(srt, cells = common_cells)

#Link Peaks to Genes
srt <- LinkPeaks(
object = srt,
peak.assay = "ATAC",
expression.assay = "RNA",
distance = 100000
    )

linked_peaks <- Links(srt)

write.csv(linked_peaks, file = "Linkpeaks.csv")



## Centering and scaling data matrix
saveRDS(srt, "06122024_combined_juvenile_integrated.allen_brain_projection.LinkPeaks.rds")





