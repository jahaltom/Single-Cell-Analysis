library(Signac)
library(Seurat)
library(dplyr)

juvenile <- readRDS("06122024_combined_juvenile_integrated.allen_brain_projection.rds")


DefaultAssay(juvenile)<-"ATAC"
#mononucleosomal / nucleosome-free ratio.  low NucleosomeSignal indicates a higher proportion of fragments in nucleosome-free regions (open chromatin, typically active regulatory regions).
juvenile <- NucleosomeSignal(juvenile)
# compute TSS enrichment score per cell. ratio of fragments centered at the TSS to fragments in TSS-flanking regions. Higher the better. 
juvenile <- TSSEnrichment(object = juvenile)

#Total number of fragments in peaks: A measure of cellular sequencing depth / complexity. Cells with very few reads may need to be excluded due to low sequencing depth. Cells with extremely high levels may represent doublets, nuclei clumps, or other artefacts.

# add fraction of reads in peaks. Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with low values (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed. Note that this value can be sensitive to the set of peaks used.
#juvenile$pct_reads_in_peaks <- juvenile$peak_region_fragments / juvenile$passed_filters * 100

# add blacklist ratio
juvenile$blacklist_ratio <- FractionCountsInRegion(
  object = juvenile, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)
#Note that the last three metrics can be obtained from the output of CellRanger

png("DensityScatter_QC_TSS_vs_nCount_peaks.png",width=15,height=15,units="in",res=300)
DensityScatter(juvenile, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()



#'pct_reads_in_peaks'
png("VlnPlot_QC.png",width=15,height=15,units="in",res=300)
VlnPlot(
  object = juvenile,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


juvenile
juvenile <- subset(
  x = juvenile,
  subset = nCount_peaks >= 1000 &
    nCount_peaks <= 25000 &
   # pct_reads_in_peaks > 40 &
    blacklist_ratio <= 0.05 &
    nucleosome_signal < 10 &
    TSS.enrichment >= 2
)
juvenile


saveRDS(juvenile, "06122024_combined_juvenile_integrated.allen_brain_projection.filtered.rds")











#Signac performs term frequency-inverse document frequency (TF-IDF) normalization. This is a two-step normalization procedure, that both normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.
juvenile <- RunTFIDF(juvenile)

juvenile <- FindTopFeatures(juvenile, min.cutoff = 'q0')
#ingular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above
juvenile <- RunSVD(juvenile)

# correlation between each LSI component and sequencing depth
png(""LSIvsSeqDepth.png",width=15,height=15,units="in",res=300)
DepthCor(juvenile)
dev.off()







#Non-linear dimension reduction and clustering

juvenile <- RunUMAP(object = juvenile, reduction = 'lsi', dims = 2:30)
juvenile <- FindNeighbors(object = juvenile, reduction = 'lsi', dims = 2:30)
juvenile <- FindClusters(object = juvenile, verbose = FALSE, algorithm = 3)


png("Cluster.png",width=15,height=15,units="in",res=300)
DimPlot(object = juvenile, label = TRUE) + NoLegend()
dev.off()










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

