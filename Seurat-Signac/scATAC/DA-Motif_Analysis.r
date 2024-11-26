library(Signac)
library(ggplot2)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
set.seed(1234)






srt <-readRDS("ATAC_AllAssays_Final.rds")

#Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
srt <- RunTFIDF(srt, assay = "ATAC")
#Run partial singular value decomposition using irlba. Dimensionality reduction
srt <- RunSVD(srt, assay = "ATAC")

DefaultAssay(srt)<-"ATAC"




#Create metadata that has cell type and condition
srt$celltype.stim <- paste(srt$MarkerAnnotations, srt$Status, sep = "_")
#Lable Idents with this metadata
Idents(srt) <- "celltype.stim"





##Motif Analysis
#to avoid an error: https://github.com/timoast/signac/issues/486#issuecomment-788964712
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(srt)) %in% main.chroms)
srt <- srt[keep.peaks, ]

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
srt <- AddMotifs(
  object = srt,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)




#Run DA and motif analysis
getDaMotif <- function(cell,condition,condition2) {


        Vcell <- FindMarkers(
          object = srt,
          ident.1 = paste(cell,condition,sep="_"),
          ident.2 = paste(cell,condition2,sep="_"),
          min.pct = 0.1,
          test.use = 'LR',
          assay="ATAC"
        )

        # Get peak coordinates
        peak_ranges <- StringToGRanges(rownames(Vcell), sep = c("-", "-"))

        # Annotate peaks with the closest gene using the genome annotations
        closest_genes <- ClosestFeature(srt, regions = peak_ranges)

        # Merge the gene annotations with the markers data frame

        Vcell$query_region=rownames(Vcell)

        Vcell <- merge(Vcell, closest_genes,by="query_region")

        write.csv(Vcell, file=paste("DA/",cell,"_",condition,"VS",condition2,".da.csv",sep=""))



        ##Motif Analysis
        Vcell <- read_csv(paste("DA/",cell,"_",condition,"VS",condition2,".da.csv",sep=""))


        Vcell<-column_to_rownames(Vcell, var = "query_region")

        top.da.peak <- rownames(Vcell[Vcell$p_val < 0.005, ])


        enriched.motifs <- FindMotifs(
          object = srt,
          features = top.da.peak
        )

        write.csv(enriched.motifs, file=paste("motifs/",cell,"_",condition,"VS",condition2,".motif.csv",sep=""))



        MotifPlot(
          object = srt,
          motifs = head(rownames(enriched.motifs))
        )
        ggsave(paste("motifs/",cell,"_",condition,"VS",condition2,".motif.png",sep=""),scale=1,width=20)



}






getDaMotif("CELL","Severe_Followup","Healthy")
getDaMotif("CELL","Severe","Healthy")
getDaMotif("CELL","Mild_Followup","Healthy")
getDaMotif("CELL","Mild","Healthy")
















