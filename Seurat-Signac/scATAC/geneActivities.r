library(Signac)
library(Seurat)
library(dplyr)

juvenile <- readRDS("06122024_combined_juvenile_integrated.allen_brain_projection.rds")


DefaultAssay(juvenile)<-"ATAC"



# Compute GeneActivity scores (this uses peaks linked to genes)
activities <- GeneActivity(juvenile, assay = "ATAC")
# Ensure GeneActivity scores are added as a new assay
juvenile[["GeneActivity"]] <- CreateAssayObject(counts = activities)



juvenile <- NormalizeData(
  object = juvenile,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(juvenile$nCount_RNA)
)


saveRDS(juvenile, "06122024_combined_juvenile_integrated.allen_brain_projection.gene-activities.rds")


###
DefaultAssay(juvenile) <- 'GeneActivity'


png("MarkerGeneCLuster.png",width=15,height=15,units="in",res=300)
FeaturePlot( 
  object = juvenile,
  features = c('Gfap', 'Mbp', 'Cx3cr1', 'Cldn5', 'Sox2', 'Foxj1','Prox1','Slc1a3'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()



Idents(juvenile)<-"HPF_class"
# Run differential analysis on GeneActivity (for example, between two clusters)
Markers_ATAC <- FindMarkers(juvenile, assay = "GeneActivity")


write.csv(Markers_ATAC, file = "GeneActivityATAC.csv")

top10 <- Markers_ATAC %>% group_by(HPF_class) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "GeneActivityATAC_top10.csv")



## Centering and scaling data matrix
saveRDS(juvenile, "06122024_combined_juvenile_integrated.allen_brain_projection.gene-activities.rds")













###Find differentially accessible peaks between cell types

